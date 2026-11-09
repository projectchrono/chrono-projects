// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2026 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Author: Radu Serban
// =============================================================================
//
// Steady shear of glass beads in an annular cell with Chrono::Multicore (SMC), for the mu(I) rheology. The beads and
// their contact settings are those of the direct shear programs (6 mm glass beads, after Hartl and Ooi 2008).
//
// The beads fill the ring between two coaxial frictionless walls (a solid inner cylinder and an outer ring of box
// plates) and are sheared between two rough plates: flat disks covered with a random layer of glued spheres of the bead
// size. The plates turn about the axis in opposite directions, each at half the relative angular speed, so that the
// mean azimuthal velocity is zero at mid-height. A stress servo on the height of the top plate holds the normal
// pressure. Gravity acts only while the beads settle; without gravity during the test, the top plate alone sets the
// pressure.
//
// The run has three phases:
// - fill        : the beads, generated at random positions in a tall annular column, settle under gravity
// - consolidate : at t_place, gravity is switched off and the top plate is placed just above the beads; it moves down at
//                 approach_speed until it carries a tenth of the pressure, then its servo holds the pressure
// - shear       : from t_shear, the plates turn (the relative angular speed ramps up over a nominal strain of
//                 strain_ramp) until the nominal strain reaches strain_end
//
// The relative angular speed follows from a nominal inertial number at the mid radius r_m = (r_in + r_out) / 2:
//   Omega = I_nominal * sqrt(pressure / density) * H / (d * r_m),
// with d the mean bead diameter and H the distance between the centers of the two glued layers at the start of the
// shear. The nominal shear rate at radius r is Omega * r / H, and the nominal strain is the relative rotation times
// r_m / H. The local shear rate, stresses, solid fraction and inertial number come from the snapshots, by
// coarse-graining.
//
// The global reference frame has Y up, along the axis of the cell; the face of the bottom plate is at y = 0. All units
// SI. Forces and torques on a plate are those of the beads on it; torques are about the +y axis.
//
// Outputs, in the output directory:
//   meta.json           all parameters, run information and summary
//   timeseries.csv      one row per output interval (written as the run goes)
//   particles/*.csv     particle snapshots, at the nominal strains snap_from, snap_from + snap_every, ...
//   contacts/*.csv      contact snapshots, at the same strains
//   states/*.bin        the state at each snapshot (with --snap_states)
// A run can start from a settled sample (--save_fill, --load_fill: one filling for several runs). The shear can start
// from a consolidated sample (--save_state saves it at t_shear, --prepare_only then stops the run; --load_state starts
// the shear from it, with this run's --I_nominal) or resume from a snapshot state (--load_state).
// =============================================================================

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "chrono/ChConfig.h"
#include "chrono/ChVersion.h"
#include "chrono/assets/ChVisualShapeSphere.h"
#include "chrono/core/ChDataPath.h"
#include "chrono/physics/ChContactMaterial.h"
#include "chrono/utils/ChUtilsCreators.h"

#include "chrono_multicore/physics/ChSystemMulticore.h"

#include "chrono_thirdparty/cxxopts/ChCLI.h"

#include "dem_utils.h"

using namespace chrono;
using namespace dem_validation;

using std::cout;
using std::endl;

// -----------------------------------------------------------------------------
// Parameters
// -----------------------------------------------------------------------------

struct Params {
    // Beads
    double radius;         // mean bead radius, m
    double radius_spread;  // half-width of a uniform radius distribution, relative to the mean (0: monodisperse)
    int num_beads;         // beads in the cell
    double density;        // kg/m^3
    double E;              // Young's modulus (beads, plates and walls), Pa
    double nu;             // Poisson ratio
    double cor;            // coefficient of restitution
    double mu;             // friction of the beads and of the plates (glued spheres and disks)
    double mu_fill;        // bead friction while the beads settle

    // Cell
    double r_in;          // radius of the inner wall, m
    double r_out;         // radius of the outer wall, m
    double h_wall;        // height of the side walls, m
    double wall;          // thickness of the plates of the outer wall, m
    int plates;           // plates in the outer ring
    double h_plate;       // thickness of the plate disks, m
    double glued_gap;     // smallest gap between the glued spheres of a plate, m
    double gen_fraction;  // solid fraction of the column in which the beads are generated

    // Loading
    double pressure;        // normal pressure held by the top plate, Pa
    double I_nominal;       // nominal inertial number at the mid radius
    double strain_end;      // nominal strain at the end of the shear
    double strain_ramp;     // nominal strain over which the angular speed ramps up
    double gravity;         // while the beads settle, m/s^2
    double t_place;         // switch gravity off and place the top plate, s
    double t_shear;         // start of the shear, s
    double approach_speed;  // speed of the top plate until it carries a tenth of the pressure, m/s
    double servo_gain;      // speed of the top plate per unit stress error, m/s/Pa
    double servo_vmax;      // largest speed of the top plate, m/s
    std::string history;    // tangential displacement model: multistep, onestep, none

    // Numerics
    double dt;
    int threads;
    int seed;
    std::vector<double> bins;  // collision bins per axis (x, y, z)

    // Output
    std::string out_dir;
    double out_dt;           // time series interval before the shear, s
    double out_strain;       // time series interval during the shear, nominal strain
    double snap_from;        // nominal strain of the first snapshot
    double snap_every;       // nominal strain between snapshots
    std::string save_fill;   // file for the settled sample
    std::string load_fill;   // start from a settled sample
    bool fill_only;          // stop after saving the settled sample
    std::string save_state;  // file for the consolidated sample, at the start of the shear
    std::string load_state;  // start the shear from a consolidated sample, or resume it from a snapshot state
    bool prepare_only;       // stop after saving the consolidated sample
    bool snap_states;        // also save the state at each snapshot
    std::string provenance;  // JSON text added to meta.json (git state, run command)
};

ChSystemSMC::TangentialDisplacementModel HistoryModel(const std::string& name) {
    if (name == "multistep")
        return ChSystemSMC::TangentialDisplacementModel::MultiStep;
    if (name == "onestep")
        return ChSystemSMC::TangentialDisplacementModel::OneStep;
    if (name == "none")
        return ChSystemSMC::TangentialDisplacementModel::None;
    std::cerr << "Unknown tangential displacement model: " << name << endl;
    std::exit(1);
}

bool GetParams(int argc, char* argv[], Params& p) {
    ChCLI cli(argv[0], "Steady shear of glass beads in an annular cell, for mu(I) (Chrono::Multicore, SMC)");

    cli.AddOption<double>("Beads", "radius", "Mean bead radius (m)", "0.003");
    cli.AddOption<double>("Beads", "radius_spread", "Half-width of a uniform radius distribution, relative to the mean (0: monodisperse)", "0");
    cli.AddOption<int>("Beads", "num_beads", "Beads in the cell", "14400");
    cli.AddOption<double>("Beads", "density", "Bead density (kg/m3)", "2550");
    cli.AddOption<double>("Beads", "E", "Young's modulus of beads, plates and walls (Pa)", "4e7");
    cli.AddOption<double>("Beads", "nu", "Poisson ratio", "0.22");
    cli.AddOption<double>("Beads", "cor", "Coefficient of restitution", "0.87");
    cli.AddOption<double>("Beads", "mu", "Friction coefficient of the beads and of the rough plates (glass on glass)", "0.18");
    cli.AddOption<double>("Beads", "mu_fill", "Bead friction while the beads settle (default: mu)");

    cli.AddOption<double>("Cell", "r_in", "Radius of the inner wall (m)", "0.09");
    cli.AddOption<double>("Cell", "r_out", "Radius of the outer wall (m)", "0.15");
    cli.AddOption<double>("Cell", "h_wall", "Height of the side walls (m)", "0.2");
    cli.AddOption<double>("Cell", "wall", "Thickness of the plates of the outer wall (m)", "0.01");
    cli.AddOption<int>("Cell", "plates", "Plates in the outer ring", "180");
    cli.AddOption<double>("Cell", "h_plate", "Thickness of the plate disks (m)", "0.01");
    cli.AddOption<double>("Cell", "glued_gap", "Smallest gap between the glued spheres of a plate (m)", "1e-4");
    cli.AddOption<double>("Cell", "gen_fraction", "Solid fraction of the column in which the beads are generated", "0.3");

    cli.AddOption<double>("Loading", "pressure", "Normal pressure held by the top plate (Pa)", "3100");
    cli.AddOption<double>("Loading", "I_nominal", "Nominal inertial number at the mid radius", "0.01");
    cli.AddOption<double>("Loading", "strain_end", "Nominal strain at the end of the shear", "5");
    cli.AddOption<double>("Loading", "strain_ramp", "Nominal strain over which the angular speed ramps up", "0.1");
    cli.AddOption<double>("Loading", "gravity", "Gravitational acceleration while the beads settle (m/s2)", "9.81");
    cli.AddOption<double>("Loading", "t_place", "Switch gravity off and place the top plate (s)", "0.8");
    cli.AddOption<double>("Loading", "t_shear", "Start of the shear (s)", "1.6");
    cli.AddOption<double>("Loading", "approach_speed", "Speed of the top plate until it carries a tenth of the pressure (m/s)", "0.02");
    cli.AddOption<double>("Loading", "servo_gain", "Speed of the top plate per unit stress error (m/s/Pa)", "3e-6");
    cli.AddOption<double>("Loading", "servo_vmax", "Largest speed of the top plate (m/s)", "0.1");
    cli.AddOption<std::string>("Loading", "history", "Tangential displacement model: multistep, onestep or none", "multistep");

    cli.AddOption<double>("Numerics", "dt", "Time step (s)", "1e-5");
    cli.AddOption<int>("Numerics", "threads", "Number of OpenMP threads", "4");
    cli.AddOption<int>("Numerics", "seed", "Seed of the random bead positions and sizes and of the glued layers", "1");
    cli.AddOption<std::vector<double>>("Numerics", "bins", "Collision bins per axis (x, y, z)", "40,24,40");

    cli.AddOption<std::string>("Output", "out", "Output directory", "annularShear");
    cli.AddOption<double>("Output", "out_dt", "Time series interval before the shear (s)", "0.01");
    cli.AddOption<double>("Output", "out_strain", "Time series interval during the shear (nominal strain)", "0.01");
    cli.AddOption<double>("Output", "snap_from", "Nominal strain of the first snapshot", "1");
    cli.AddOption<double>("Output", "snap_every", "Nominal strain between snapshots", "0.25");
    cli.AddOption<std::string>("Output", "save_fill", "Save the settled sample (before the top plate is placed) to this file", "");
    cli.AddOption<std::string>("Output", "load_fill", "Start from a settled sample saved with --save_fill", "");
    cli.AddOption<bool>("Output", "fill_only", "Stop after the settled sample is saved");
    cli.AddOption<std::string>("Output", "save_state", "Save the consolidated sample, at the start of the shear, to this file", "");
    cli.AddOption<std::string>("Output", "load_state", "Start the shear from a consolidated sample, or resume it from a snapshot state", "");
    cli.AddOption<bool>("Output", "prepare_only", "Stop after the consolidated sample is saved");
    cli.AddOption<bool>("Output", "snap_states", "Also save the state at each snapshot (states/*.bin)");
    cli.AddOption<std::string>("Output", "provenance", "JSON text recorded in meta.json", "");

    if (!cli.Parse(argc, argv, true))
        return false;

    p.radius = cli.GetAsType<double>("radius");
    p.radius_spread = cli.GetAsType<double>("radius_spread");
    p.num_beads = cli.GetAsType<int>("num_beads");
    p.density = cli.GetAsType<double>("density");
    p.E = cli.GetAsType<double>("E");
    p.nu = cli.GetAsType<double>("nu");
    p.cor = cli.GetAsType<double>("cor");
    p.mu = cli.GetAsType<double>("mu");
    p.mu_fill = cli.CheckOption("mu_fill") ? cli.GetAsType<double>("mu_fill") : p.mu;

    p.r_in = cli.GetAsType<double>("r_in");
    p.r_out = cli.GetAsType<double>("r_out");
    p.h_wall = cli.GetAsType<double>("h_wall");
    p.wall = cli.GetAsType<double>("wall");
    p.plates = cli.GetAsType<int>("plates");
    p.h_plate = cli.GetAsType<double>("h_plate");
    p.glued_gap = cli.GetAsType<double>("glued_gap");
    p.gen_fraction = cli.GetAsType<double>("gen_fraction");

    p.pressure = cli.GetAsType<double>("pressure");
    p.I_nominal = cli.GetAsType<double>("I_nominal");
    p.strain_end = cli.GetAsType<double>("strain_end");
    p.strain_ramp = cli.GetAsType<double>("strain_ramp");
    p.gravity = cli.GetAsType<double>("gravity");
    p.t_place = cli.GetAsType<double>("t_place");
    p.t_shear = cli.GetAsType<double>("t_shear");
    p.approach_speed = cli.GetAsType<double>("approach_speed");
    p.servo_gain = cli.GetAsType<double>("servo_gain");
    p.servo_vmax = cli.GetAsType<double>("servo_vmax");
    p.history = cli.GetAsType<std::string>("history");
    HistoryModel(p.history);

    p.dt = cli.GetAsType<double>("dt");
    p.threads = cli.GetAsType<int>("threads");
    p.seed = cli.GetAsType<int>("seed");
    p.bins = cli.GetAsType<std::vector<double>>("bins");
    if (p.bins.size() != 3 || p.bins[0] < 1 || p.bins[1] < 1 || p.bins[2] < 1) {
        std::cerr << "--bins needs three positive counts" << endl;
        return false;
    }

    p.out_dir = cli.GetAsType<std::string>("out");
    p.out_dt = cli.GetAsType<double>("out_dt");
    p.out_strain = cli.GetAsType<double>("out_strain");
    p.snap_from = cli.GetAsType<double>("snap_from");
    p.snap_every = cli.GetAsType<double>("snap_every");
    p.save_fill = cli.GetAsType<std::string>("save_fill");
    p.load_fill = cli.GetAsType<std::string>("load_fill");
    p.fill_only = cli.GetAsType<bool>("fill_only");
    p.save_state = cli.GetAsType<std::string>("save_state");
    p.load_state = cli.GetAsType<std::string>("load_state");
    p.prepare_only = cli.GetAsType<bool>("prepare_only");
    p.snap_states = cli.GetAsType<bool>("snap_states");
    p.provenance = cli.GetAsType<std::string>("provenance");

    if (p.num_beads <= 0 || p.radius <= 0 || p.radius_spread < 0 || p.radius_spread >= 1) {
        std::cerr << "Need beads with a positive radius and a spread below 1" << endl;
        return false;
    }
    if (!(p.r_in > 0 && p.r_out > p.r_in + 4 * p.radius)) {
        std::cerr << "Need 0 < r_in and r_out > r_in + 2 diameters" << endl;
        return false;
    }
    if (!(p.t_place < p.t_shear) || p.approach_speed <= 0) {
        std::cerr << "Need t_place < t_shear and a positive approach speed" << endl;
        return false;
    }
    if (p.I_nominal <= 0 || p.strain_end <= 0 || p.strain_ramp < 0 || p.strain_ramp >= p.strain_end || p.snap_every <= 0 || p.out_strain <= 0) {
        std::cerr << "Need I_nominal, strain_end, snap_every and out_strain positive, and 0 <= strain_ramp < strain_end" << endl;
        return false;
    }
    if (p.fill_only && p.save_fill.empty()) {
        std::cerr << "--fill_only needs --save_fill" << endl;
        return false;
    }
    if (p.prepare_only && p.save_state.empty()) {
        std::cerr << "--prepare_only needs --save_state" << endl;
        return false;
    }

    return true;
}

// -----------------------------------------------------------------------------
// Random sequential addition of the glued spheres of a plate: centers in the annulus between the walls, at least one
// diameter plus the gap apart, as many as fit in the given number of attempts
// -----------------------------------------------------------------------------

std::vector<ChVector3d> GluedLayer(std::mt19937& rng, double r_in, double r_out, double r, double gap, int attempts) {
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    double rho_min = r_in + r + gap;
    double rho_max = r_out - r - gap;
    double dmin = 2 * r + gap;
    int n = (int)std::ceil(2 * r_out / dmin) + 1;
    std::vector<std::vector<int>> grid(n * n);
    std::vector<ChVector3d> centers;
    for (int a = 0; a < attempts; a++) {
        double rho = std::sqrt(rho_min * rho_min + (rho_max * rho_max - rho_min * rho_min) * uni(rng));
        double th = 2 * CH_PI * uni(rng);
        double x = rho * std::cos(th);
        double z = rho * std::sin(th);
        int ci = (int)std::floor((x + r_out) / dmin);
        int ck = (int)std::floor((z + r_out) / dmin);
        bool ok = true;
        for (int i = ci - 1; i <= ci + 1 && ok; i++) {
            for (int k = ck - 1; k <= ck + 1 && ok; k++) {
                if (i < 0 || k < 0 || i >= n || k >= n)
                    continue;
                for (int m : grid[i * n + k]) {
                    double dx = centers[m].x() - x;
                    double dz = centers[m].z() - z;
                    if (dx * dx + dz * dz < dmin * dmin) {
                        ok = false;
                        break;
                    }
                }
            }
        }
        if (ok) {
            grid[ci * n + ck].push_back((int)centers.size());
            centers.push_back(ChVector3d(x, 0, z));
        }
    }
    return centers;
}

// Values saved with a state, after the bodies: phase, start of the shear, relative angular speed, gap at the start of
// the shear, nominal inertial number, total bead volume (to check that a state matches the beads)
static const size_t kStateExtra = 6;

// =============================================================================

int main(int argc, char* argv[]) {
    Params p;
    if (!GetParams(argc, argv, p))
        return 1;

    SetChronoDataPath(CHRONO_DATA_DIR);

    const std::string out_dir = p.out_dir;
    const std::string particles_dir = out_dir + "/particles";
    const std::string contacts_dir = out_dir + "/contacts";
    const std::string states_dir = out_dir + "/states";
    std::vector<std::string> dirs = {out_dir, particles_dir, contacts_dir};
    if (p.snap_states)
        dirs.push_back(states_dir);
    for (const auto& dir : dirs) {
        std::error_code ec;
        std::filesystem::create_directories(dir, ec);
        if (!std::filesystem::is_directory(dir)) {
            std::cerr << "Error creating directory " << dir << endl;
            return 1;
        }
    }

    auto wall_start = std::chrono::steady_clock::now();
    std::string start_date = Now();

    const double r_glued = p.radius;
    const double r_mid = 0.5 * (p.r_in + p.r_out);
    const double area = CH_PI * (p.r_out * p.r_out - p.r_in * p.r_in);
    const double torque_arm = 2 * CH_PI * (std::pow(p.r_out, 3) - std::pow(p.r_in, 3)) / 3;  // torque of a uniform unit shear stress

    // Bead sizes and the glued layers
    std::mt19937 rng((unsigned int)p.seed);
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    std::vector<double> radius_of(p.num_beads);
    for (auto& r : radius_of)
        r = p.radius * (1 + p.radius_spread * (2 * uni(rng) - 1));
    double r_max = *std::max_element(radius_of.begin(), radius_of.end());
    double d_mean = 0;
    double solid_volume = 0;
    for (double r : radius_of) {
        d_mean += 2 * r;
        solid_volume += (4.0 / 3) * CH_PI * r * r * r;
    }
    d_mean /= p.num_beads;

    auto glued_bottom = GluedLayer(rng, p.r_in, p.r_out, r_glued, p.glued_gap, 200000);
    auto glued_top = GluedLayer(rng, p.r_in, p.r_out, r_glued, p.glued_gap, 200000);
    double glued_coverage = glued_bottom.size() * CH_PI * r_glued * r_glued / area;

    // Annular column in which the beads are generated, above the bottom glued layer. The column holds the beads, each
    // with half the generation gap around it, at the solid fraction gen_fraction.
    double gen_gap = 0.0002;
    double rho_lo = p.r_in + r_max + gen_gap;
    double rho_hi = p.r_out - r_max - gen_gap;
    double y_lo = 2 * r_glued + gen_gap;
    double gen_volume = 0;
    for (double r : radius_of)
        gen_volume += (4.0 / 3) * CH_PI * std::pow(r + gen_gap / 2, 3);
    double h_gen = gen_volume / (p.gen_fraction * CH_PI * (rho_hi * rho_hi - rho_lo * rho_lo));
    if (y_lo + h_gen > p.h_wall - 2 * r_max) {
        std::cerr << "The generation column (" << y_lo + h_gen << " m) is taller than the side walls; use a larger --h_wall or --gen_fraction" << endl;
        return 1;
    }

    // Random positions in the column, apart from one another (random sequential addition, with a grid of cells one
    // largest diameter wide)
    std::vector<ChVector3d> positions;
    {
        double cell = 2 * r_max + gen_gap;
        int nxz = (int)std::ceil(2 * p.r_out / cell) + 1;
        int ny = (int)std::ceil(h_gen / cell) + 1;
        std::vector<std::vector<int>> grid(nxz * nxz * ny);
        auto cell_index = [&](int i, int j, int k) { return (j * nxz + i) * nxz + k; };
        for (int n = 0; n < p.num_beads; n++) {
            double r = radius_of[n];
            bool placed = false;
            for (int attempt = 0; attempt < 100000 && !placed; attempt++) {
                double rho = std::sqrt(rho_lo * rho_lo + (rho_hi * rho_hi - rho_lo * rho_lo) * uni(rng));
                double th = 2 * CH_PI * uni(rng);
                ChVector3d q(rho * std::cos(th), y_lo + r + (h_gen - 2 * r) * uni(rng), rho * std::sin(th));
                int ci = (int)std::floor((q.x() + p.r_out) / cell);
                int cj = (int)std::floor((q.y() - y_lo) / cell);
                int ck = (int)std::floor((q.z() + p.r_out) / cell);
                bool ok = true;
                for (int i = ci - 1; i <= ci + 1 && ok; i++) {
                    for (int j = cj - 1; j <= cj + 1 && ok; j++) {
                        for (int k = ck - 1; k <= ck + 1 && ok; k++) {
                            if (i < 0 || j < 0 || k < 0 || i >= nxz || j >= ny || k >= nxz)
                                continue;
                            for (int m : grid[cell_index(i, j, k)]) {
                                double dmin = r + radius_of[m] + gen_gap;
                                if ((positions[m] - q).Length2() < dmin * dmin) {
                                    ok = false;
                                    break;
                                }
                            }
                        }
                    }
                }
                if (ok) {
                    grid[cell_index(ci, cj, ck)].push_back(n);
                    positions.push_back(q);
                    placed = true;
                }
            }
            if (!placed) {
                std::cerr << "Cannot place bead " << n << " in the column; use a smaller --gen_fraction" << endl;
                return 1;
            }
        }
    }

    // Rotation of -90 degrees around x: z axis to y axis
    ChQuaterniond z2y;
    z2y.SetFromAngleAxis(-CH_PI / 2, ChVector3d(1, 0, 0));

    // Create the system

    ChSystemMulticoreSMC* sys = new ChSystemMulticoreSMC();
    sys->SetCollisionSystemType(ChCollisionSystem::Type::MULTICORE);
    sys->SetGravitationalAcceleration(ChVector3d(0, -p.gravity, 0));

    int threads = std::min(p.threads, omp_get_num_procs());
    sys->SetNumThreads(threads);

    sys->GetSettings()->solver.use_full_inertia_tensor = false;
    sys->GetSettings()->solver.tolerance = 0.01;
    sys->GetSettings()->solver.max_iteration_bilateral = 100;
    sys->GetSettings()->solver.clamp_bilaterals = false;
    sys->GetSettings()->solver.bilateral_clamp_speed = 0.1;

    sys->GetSettings()->solver.contact_force_model = ChSystemSMC::ContactForceModel::Hertz;
    sys->GetSettings()->solver.tangential_displ_mode = HistoryModel(p.history);

    sys->GetSettings()->collision.bins_per_axis = vec3((int)p.bins[0], (int)p.bins[1], (int)p.bins[2]);
    sys->GetSettings()->collision.narrowphase_algorithm = ChNarrowphase::Algorithm::HYBRID;

    // Materials: beads, rough plates (glued spheres and disks, glass on glass), frictionless side walls. The friction of
    // a contact is the smaller of its two materials' values (the default composition).

    auto mat_bead = chrono_types::make_shared<ChContactMaterialSMC>();
    mat_bead->SetYoungModulus((float)p.E);
    mat_bead->SetPoissonRatio((float)p.nu);
    mat_bead->SetRestitution((float)p.cor);
    mat_bead->SetFriction((float)p.mu_fill);

    auto mat_plate = chrono_types::make_shared<ChContactMaterialSMC>();
    mat_plate->SetYoungModulus((float)p.E);
    mat_plate->SetPoissonRatio((float)p.nu);
    mat_plate->SetRestitution((float)p.cor);
    mat_plate->SetFriction((float)p.mu);

    auto mat_wall = chrono_types::make_shared<ChContactMaterialSMC>();
    mat_wall->SetYoungModulus((float)p.E);
    mat_wall->SetPoissonRatio((float)p.nu);
    mat_wall->SetRestitution((float)p.cor);
    mat_wall->SetFriction(0.0f);

    // The plates and the walls are created before the beads, so that the history of their contacts is stored on the
    // beads (the body with the larger index). Plates (family 1) and walls (family 2) do not collide with each other.

    // Bottom plate (body 0): face at y = 0, glued spheres on it
    auto bottom = chrono_types::make_shared<ChBody>();
    bottom->SetPos(ChVector3d(0, 0, 0));
    bottom->SetFixed(true);
    bottom->EnableCollision(true);
    utils::AddCylinderGeometry(bottom.get(), mat_plate, p.r_out, p.h_plate, ChVector3d(0, -p.h_plate / 2, 0), z2y, true);
    for (const auto& g : glued_bottom)
        utils::AddSphereGeometry(bottom.get(), mat_plate, r_glued, ChVector3d(g.x(), r_glued, g.z()), QUNIT, true);
    bottom->GetCollisionModel()->SetFamily(1);
    bottom->GetCollisionModel()->DisallowCollisionsWith(2);
    sys->AddBody(bottom);

    // Top plate (body 1): face (its underside) at y_face, glued spheres hanging below it; parked above the walls while
    // the beads settle
    double y_park = p.h_wall + 2 * r_glued + 0.05;
    auto top = chrono_types::make_shared<ChBody>();
    top->SetPos(ChVector3d(0, y_park, 0));
    top->SetFixed(true);
    top->EnableCollision(true);
    utils::AddCylinderGeometry(top.get(), mat_plate, p.r_out, p.h_plate, ChVector3d(0, p.h_plate / 2, 0), z2y, true);
    for (const auto& g : glued_top)
        utils::AddSphereGeometry(top.get(), mat_plate, r_glued, ChVector3d(g.x(), -r_glued, g.z()), QUNIT, true);
    top->GetCollisionModel()->SetFamily(1);
    top->GetCollisionModel()->DisallowCollisionsWith(2);
    sys->AddBody(top);

    // Inner wall (body 2): a solid cylinder
    auto inner = chrono_types::make_shared<ChBody>();
    inner->SetPos(ChVector3d(0, p.h_wall / 2, 0));
    inner->SetFixed(true);
    inner->EnableCollision(true);
    utils::AddCylinderGeometry(inner.get(), mat_wall, p.r_in, p.h_wall, ChVector3d(0, 0, 0), z2y, true);
    inner->GetCollisionModel()->SetFamily(2);
    inner->GetCollisionModel()->DisallowCollisionsWith(1);
    sys->AddBody(inner);

    // Outer wall (body 3): a ring of box plates, open at both ends
    auto outer = utils::CreateCylindricalContainerFromBoxes(sys, mat_wall, p.r_out, p.h_wall, p.wall, p.plates, ChVector3d(0, 0, 0), z2y, true, false, false, true, true, false);
    outer->GetCollisionModel()->SetFamily(2);
    outer->GetCollisionModel()->DisallowCollisionsWith(1);

    // Beads
    std::vector<std::shared_ptr<ChBody>> beads;
    int first_particle = (int)sys->GetBodies().size();
    for (int i = 0; i < p.num_beads; i++) {
        double r = radius_of[i];
        double mass = p.density * (4.0 / 3) * CH_PI * r * r * r;
        auto bead = chrono_types::make_shared<ChBody>();
        bead->SetMass(mass);
        bead->SetInertiaXX(0.4 * mass * r * r * ChVector3d(1, 1, 1));
        bead->SetPos(positions[i]);
        bead->EnableCollision(true);
        auto ct_shape = chrono_types::make_shared<ChCollisionShapeSphere>(mat_bead, r);
        bead->AddCollisionShape(ct_shape);
        bead->GetCollisionModel()->SetFamily(3);
        auto vis_shape = chrono_types::make_shared<ChVisualShapeSphere>(r);
        vis_shape->SetColor(ChColor(0.8f, 0.8f, 0.8f));
        bead->AddVisualShape(vis_shape);
        sys->AddBody(bead);
        beads.push_back(bead);
    }

    // Phase and reference values
    // phase 0: fill, 1: consolidate (top plate placed, gravity off), 2: shear

    int phase = 0;
    double y_face = NAN;  // face of the top plate (its underside)
    double v_top = 0;     // vertical speed of the top plate
    double t0 = NAN;      // start of the shear
    double Omega = 0;     // relative angular speed of the plates
    double H0 = NAN;      // distance between the centers of the glued layers at the start of the shear
    double t_ramp = 0;    // duration of the ramp of the angular speed
    double I_run = p.I_nominal;
    bool approaching = false;  // the top plate moves down at approach_speed until it carries a tenth of the pressure
    double t_contact = NAN;    // end of the approach

    auto gap = [&]() { return y_face - 2 * r_glued; };

    // Relative rotation of the plates (top minus bottom) and its rate at a time of the shear
    auto rotation = [&](double time, double& angle, double& omega) {
        double s = time - t0;
        if (phase < 2 || s <= 0) {
            angle = 0;
            omega = 0;
        } else if (s < t_ramp) {
            omega = Omega * s / t_ramp;
            angle = 0.5 * Omega * s * s / t_ramp;
        } else {
            omega = Omega;
            angle = Omega * (s - 0.5 * t_ramp);
        }
    };
    auto strain_at = [&](double time) {
        double angle, omega;
        rotation(time, angle, omega);
        return phase == 2 ? angle * r_mid / H0 : 0.0;
    };
    // Time of the shear at which the nominal strain reaches g (g at least strain_ramp)
    auto time_of_strain = [&](double g) { return t0 + g * H0 / (Omega * r_mid) + 0.5 * t_ramp; };

    auto start_shear = [&](double time) {
        t0 = time;
        H0 = gap();
        Omega = I_run * std::sqrt(p.pressure / p.density) * H0 / (d_mean * r_mid);
        t_ramp = 2 * p.strain_ramp * H0 / (Omega * r_mid);
        phase = 2;
    };

    // Servo target: the pressure, once the top plate is placed
    auto target = [&]() { return phase == 0 ? 0.0 : p.pressure; };

    // The top plate is placed with its glued spheres just above the highest bead
    auto place_top = [&]() {
        double y_max = 0;
        for (size_t i = 0; i < beads.size(); i++)
            y_max = std::max(y_max, beads[i]->GetPos().y() + radius_of[i]);
        y_face = y_max + 2 * r_glued + 0.0005;
        v_top = 0;
    };

    // Kinematics of the plates at a time
    auto move_plates = [&](double time) {
        double angle, omega;
        rotation(time, angle, omega);
        ChQuaterniond qb, qt;
        qb.SetFromAngleAxis(-0.5 * angle, ChVector3d(0, 1, 0));
        qt.SetFromAngleAxis(0.5 * angle, ChVector3d(0, 1, 0));
        bottom->SetRot(qb);
        bottom->SetAngVelParent(ChVector3d(0, -0.5 * omega, 0));
        if (phase >= 1) {
            top->SetPos(ChVector3d(0, y_face, 0));
            top->SetLinVel(ChVector3d(0, v_top, 0));
        }
        top->SetRot(qt);
        top->SetAngVelParent(ChVector3d(0, 0.5 * omega, 0));
    };

    auto save = [&](const std::string& file) { return SaveState(*sys, file, {(double)phase, t0, Omega, H0, I_run, solid_volume}); };

    // Start from a settled sample, or from a consolidated or sheared state
    if (!p.load_fill.empty() || !p.load_state.empty()) {
        const std::string& file = p.load_state.empty() ? p.load_fill : p.load_state;
        std::vector<double> extra;
        if (!LoadState(*sys, file, extra) || extra.size() != kStateExtra || std::abs(extra[5] - solid_volume) > 1e-9 * solid_volume) {
            std::cerr << "Cannot load state from " << file << " (missing file, or different beads)" << endl;
            return 1;
        }
        phase = (int)extra[0];
        if (phase >= 1) {
            y_face = top->GetPos().y();
            v_top = top->GetPosDt().y();
            sys->SetGravitationalAcceleration(ChVector3d(0, 0, 0));
            mat_bead->SetFriction((float)p.mu);
        }
        if (phase == 1 && !p.load_state.empty()) {
            // A consolidated sample: the shear starts now, with this run's nominal inertial number
            start_shear(sys->GetChTime());
        } else if (phase == 2) {
            // A snapshot state: the shear resumes with the saved angular speed
            t0 = extra[1];
            Omega = extra[2];
            H0 = extra[3];
            I_run = extra[4];
            t_ramp = 2 * p.strain_ramp * H0 / (Omega * r_mid);
        }
        cout << "Loaded state at t = " << sys->GetChTime() << " s (phase " << phase << ") from " << file << endl;
    }

    auto end_time = [&]() { return phase == 2 ? time_of_strain(p.strain_end) : 1e30; };

    // Run metadata (rewritten at the start of the shear and at the end)

    auto write_meta = [&](bool finished, double wall_time, unsigned long long overflows, const std::string& status) {
        JsonObject m;
        m.Add("program", "test_VDEM_annularShear");
        m.Add("description",
              "Steady shear in an annular cell between counter-rotating rough plates, pressure-controlled, Chrono::Multicore SMC (Hertz), Y up along the axis, SI units");
        m.Add("chrono_version", CHRONO_VERSION);
        m.Add("start", start_date);
        m.Add("status", status);
        m.Add("radius", p.radius);
        m.Add("radius_spread", p.radius_spread);
        m.Add("num_beads", p.num_beads);
        m.Add("mean_diameter", d_mean);
        m.Add("solid_volume", solid_volume);
        m.Add("density", p.density);
        m.Add("E", p.E);
        m.Add("nu", p.nu);
        m.Add("cor", p.cor);
        m.Add("mu_bead", p.mu);
        m.Add("mu_plate", p.mu);
        m.Add("mu_wall", 0.0);
        m.Add("mu_fill", p.mu_fill);
        m.Add("friction_composition", "min");
        m.Add("r_in", p.r_in);
        m.Add("r_out", p.r_out);
        m.Add("r_mid", r_mid);
        m.Add("area", area);
        m.Add("torque_arm", torque_arm);
        m.Add("h_wall", p.h_wall);
        m.Add("outer_wall_plates", p.plates);
        m.Add("h_plate", p.h_plate);
        m.Add("glued_radius", r_glued);
        m.Add("glued_gap", p.glued_gap);
        m.Add("glued_bottom", (int)glued_bottom.size());
        m.Add("glued_top", (int)glued_top.size());
        m.Add("glued_coverage", glued_coverage);
        m.Add("gen_fraction", p.gen_fraction);
        m.Add("gen_height", h_gen);
        m.Add("pressure", p.pressure);
        m.Add("I_nominal", I_run);
        m.Add("strain_end", p.strain_end);
        m.Add("strain_ramp", p.strain_ramp);
        m.Add("gravity_fill", p.gravity);
        m.Add("gravity_test", 0.0);
        m.Add("t_place", p.t_place);
        m.Add("approach_speed", p.approach_speed);
        m.Add("t_contact", t_contact);
        m.Add("t_shear", p.t_shear);
        m.Add("t_shear_start", t0);
        m.Add("t_ramp", t_ramp);
        m.Add("t_end", phase == 2 ? end_time() : NAN);
        m.Add("gap_shear_start", H0);
        m.Add("omega", Omega);
        m.Add("servo_gain", p.servo_gain);
        m.Add("servo_vmax", p.servo_vmax);
        m.Add("history", p.history);
        m.Add("contact_force_model", "Hertz");
        m.Add("narrowphase", "HYBRID");
        m.Add("dt", p.dt);
        m.Add("threads", threads);
        m.Add("seed", p.seed);
        m.Add("bins", p.bins);
        m.Add("out_dt", p.out_dt);
        m.Add("out_strain", p.out_strain);
        m.Add("snap_from", p.snap_from);
        m.Add("snap_every", p.snap_every);
        m.Add("save_fill", p.save_fill);
        m.Add("load_fill", p.load_fill);
        m.Add("save_state", p.save_state);
        m.Add("load_state", p.load_state);
        m.Add("snap_states", p.snap_states);
        m.Add("first_particle_index", first_particle);
        m.AddRaw("walls", "{\"bottom plate\": {\"body\": 0}, \"top plate\": {\"body\": 1}, \"inner wall\": {\"body\": 2}, \"outer wall\": {\"body\": 3}}");
        if (!p.provenance.empty())
            m.AddRaw("provenance", p.provenance);
        m.Add("finished", finished);
        m.Add("wall_clock_s", wall_time);
        m.Add("history_overflows_total", overflows);
        m.Write(out_dir + "/meta.json");
    };
    write_meta(false, 0, 0, "running");

    // Time series

    std::ofstream ts(out_dir + "/timeseries.csv");
    ts << "# Annular shear cell. SI units. Y up along the axis; torques about +y. phase: 0 fill, 1 consolidate, 2 shear.\n";
    ts << "# strain: nominal (relative rotation times r_mid / gap at the start of the shear); omega: relative angular speed; "
          "gap: distance between the centers of the glued layers.\n";
    ts << "# p_top = Fy(top)/area, p_bottom = -Fy(bottom)/area; tau = |torque| / torque_arm (uniform shear stress); mu = tau / "
          "p_top; I_nom: nominal inertial number at the mid radius, from the gap and the target pressure.\n";
    ts << "# phi_center: solid fraction of the beads whose centers lie in the middle half of the gap and 2 diameters or more "
          "from the side walls; f_inner, f_outer: radial contact force on the side walls.\n";
    ts << "time,phase,strain,omega,angle,gap,y_face,v_top,p_target,p_top,p_bottom,torque_top,torque_bottom,tau_top,tau_bottom,"
          "mu_top,mu_bottom,I_nom,phi_center,f_inner,f_outer,kinetic_energy,"
          "contacts,contacts_pp,contacts_pw,coord,coord_mech,rattlers,frac_sliding,mean_fn,"
          "overflows_step,overflows_total,wall_clock,step_time,"
          "t_collision,t_broad,t_narrow,t_update,t_advance,t_contacts\n";
    ts << std::setprecision(9);

    // Snapshots: at the nominal strains snap_from + k * snap_every, up to strain_end
    std::vector<double> snaps;
    for (double g = p.snap_from; g <= p.strain_end + 1e-9; g += p.snap_every)
        snaps.push_back(g);

    // Simulation loop

    double next_out_time = (std::floor(sys->GetChTime() / p.out_dt + 1e-9) + 1) * p.out_dt;
    double next_out_strain = phase == 2 ? (std::floor(strain_at(sys->GetChTime()) / p.out_strain + 1e-9) + 1) * p.out_strain : 0;
    size_t next_snap = 0;
    if (phase == 2) {
        double g0 = strain_at(sys->GetChTime());
        while (next_snap < snaps.size() && snaps[next_snap] < g0 - 1e-9)
            next_snap++;
    }
    double last_wall = 0;
    long long last_step = 0;
    long long step = 0;
    unsigned long long overflows_total = 0;

    auto write_snapshot = [&](const std::string& name, double strain) {
        auto contacts = GetContacts(*sys);
        std::vector<int> counts;
        ComputeContactStats(contacts, first_particle, p.num_beads, &counts);
        double angle, omega;
        rotation(sys->GetChTime(), angle, omega);
        std::ostringstream comment;
        comment << "time = " << std::setprecision(10) << sys->GetChTime() << " s, strain = " << strain << ", omega = " << omega << " rad/s, angle = " << angle
                << " rad, y_face = " << y_face << " m, pressure = " << p.pressure << " Pa";
        WriteParticles(particles_dir + "/" + name + ".csv", beads, radius_of, counts, comment.str());
        WriteContacts(contacts_dir + "/" + name + ".csv", contacts, comment.str());
    };

    while (sys->GetChTime() < end_time() - 1e-9) {
        double time = sys->GetChTime();

        // Consolidate: save the settled sample, switch gravity off and place the top plate
        if (phase == 0 && time >= p.t_place - 1e-9) {
            if (!p.save_fill.empty() && p.load_fill.empty()) {
                if (save(p.save_fill))
                    cout << "Saved settled sample at t = " << time << " s to " << p.save_fill << endl;
                else
                    std::cerr << "Cannot save state to " << p.save_fill << endl;
            }
            if (p.fill_only) {
                double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall_start).count();
                write_meta(true, wall, overflows_total, "filled");
                cout << "Stopped after filling" << endl;
                return 0;
            }
            phase = 1;
            sys->SetGravitationalAcceleration(ChVector3d(0, 0, 0));
            mat_bead->SetFriction((float)p.mu);
            place_top();
            approaching = true;
        }

        // Shear: from the consolidated sample
        if (phase == 1 && time >= p.t_shear - 1e-9) {
            if (!p.save_state.empty()) {
                if (save(p.save_state))
                    cout << "Saved consolidated sample at t = " << time << " s to " << p.save_state << endl;
                else
                    std::cerr << "Cannot save state to " << p.save_state << endl;
            }
            if (p.prepare_only) {
                double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall_start).count();
                write_meta(true, wall, overflows_total, "consolidated");
                cout << "Stopped after consolidating" << endl;
                return 0;
            }
            start_shear(time);
            next_out_strain = p.out_strain;
            write_meta(false, 0, overflows_total, "running");
        }

        // Top plate: the approach, then the stress servo on its height
        if (phase >= 1) {
            const auto& map = sys->data_manager->host_data.ct_body_map;
            double p_now = (map.size() > top->GetIndex() ? sys->GetBodyContactForce(top).y : 0.0) / area;
            if (approaching && p_now >= 0.1 * p.pressure) {
                approaching = false;
                t_contact = time;
                cout << "Top plate in contact at t = " << time << " s" << endl;
            }
            v_top = approaching ? -p.approach_speed : std::clamp(-p.servo_gain * (target() - p_now), -p.servo_vmax, p.servo_vmax);
            y_face += v_top * p.dt;
        }
        move_plates(time);

        sys->DoStepDynamics(p.dt);
        step++;
        overflows_total = sys->GetTotalContactHistoryOverflows();

        time = sys->GetChTime();
        double strain = strain_at(time);

        // Snapshots
        if (phase == 2) {
            while (next_snap < snaps.size() && strain >= snaps[next_snap] - 1e-9) {
                char name[64];
                snprintf(name, sizeof(name), "snap_g%05.2f", snaps[next_snap]);
                write_snapshot(name, strain);
                if (p.snap_states && !save(states_dir + "/" + name + ".bin"))
                    std::cerr << "Cannot save state to " << states_dir + "/" + name + ".bin" << endl;
                next_snap++;
            }
        }

        // Time series: every out_dt before the shear, every out_strain during it
        bool output = phase == 2 ? strain >= next_out_strain - 1e-9 : time >= next_out_time - 1e-9;
        if (output) {
            if (phase == 2) {
                while (next_out_strain <= strain + 1e-9)
                    next_out_strain += p.out_strain;
            } else {
                while (next_out_time <= time + 1e-9)
                    next_out_time += p.out_dt;
            }

            auto contacts = GetContacts(*sys);
            auto st = ComputeContactStats(contacts, first_particle, p.num_beads);

            // Contact forces and torques on the plates, radial forces on the walls
            ChVector3d f_bottom(0), f_top(0);
            double t_bottom = 0, t_top = 0, f_inner = 0, f_outer = 0;
            for (const auto& c : contacts) {
                // Force on body 2 is F; on body 1 it is -F
                ChVector3d F = c.fn * c.n + c.ft;
                ChVector3d pc = 0.5 * (c.p1 + c.p2);
                auto add = [&](int b, const ChVector3d& f) {
                    double torque_y = pc.z() * f.x() - pc.x() * f.z();
                    double rho = std::sqrt(pc.x() * pc.x() + pc.z() * pc.z());
                    double f_radial = rho > 0 ? (pc.x() * f.x() + pc.z() * f.z()) / rho : 0.0;
                    if (b == 0) {
                        f_bottom += f;
                        t_bottom += torque_y;
                    } else if (b == 1) {
                        f_top += f;
                        t_top += torque_y;
                    } else if (b == 2) {
                        f_inner += f_radial;
                    } else if (b == 3) {
                        f_outer += f_radial;
                    }
                };
                add(c.b1, -F);
                add(c.b2, F);
            }

            double p_top = phase >= 1 ? f_top.y() / area : NAN;
            double p_bottom = -f_bottom.y() / area;
            double tau_top = std::abs(t_top) / torque_arm;
            double tau_bottom = std::abs(t_bottom) / torque_arm;
            double angle, omega;
            rotation(time, angle, omega);
            double H = phase >= 1 ? gap() : NAN;
            double I_nom = phase == 2 ? omega * r_mid / H * d_mean / std::sqrt(p.pressure / p.density) : 0.0;

            // Solid fraction in the middle half of the gap, 2 diameters or more from the side walls (beads counted by
            // their centers)
            double phi_center = NAN;
            if (phase >= 1) {
                double y_a = 0.25 * y_face;
                double y_b = 0.75 * y_face;
                double rho_a = p.r_in + 2 * d_mean;
                double rho_b = p.r_out - 2 * d_mean;
                double v_solid = 0;
                for (size_t i = 0; i < beads.size(); i++) {
                    const auto& c = beads[i]->GetPos();
                    double rho2 = c.x() * c.x() + c.z() * c.z();
                    if (c.y() > y_a && c.y() < y_b && rho2 > rho_a * rho_a && rho2 < rho_b * rho_b)
                        v_solid += (4.0 / 3) * CH_PI * radius_of[i] * radius_of[i] * radius_of[i];
                }
                phi_center = v_solid / (CH_PI * (rho_b * rho_b - rho_a * rho_a) * (y_b - y_a));
            }

            double ke = 0;
            for (size_t i = 0; i < beads.size(); i++) {
                const auto& b = beads[i];
                ke += 0.5 * b->GetMass() * b->GetPosDt().Length2();
                ke += 0.5 * 0.4 * b->GetMass() * radius_of[i] * radius_of[i] * b->GetAngVelParent().Length2();
            }

            double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall_start).count();
            double step_time = (step > last_step) ? (wall - last_wall) / (step - last_step) : 0;
            last_wall = wall;
            last_step = step;

            ts << time << "," << phase << "," << strain << "," << omega << "," << angle << "," << H << "," << y_face << "," << v_top << "," << target() << "," << p_top << ","
               << p_bottom << "," << t_top << "," << t_bottom << "," << tau_top << "," << tau_bottom << "," << tau_top / p_top << "," << tau_bottom / p_top << "," << I_nom << ","
               << phi_center << "," << f_inner << "," << f_outer << "," << ke << "," << st.num_contacts << "," << st.num_pp << "," << st.num_pw << "," << st.coord << ","
               << st.coord_mech << "," << st.num_rattlers << "," << st.frac_sliding << "," << st.mean_fn << "," << sys->GetNumContactHistoryOverflows() << "," << overflows_total
               << "," << wall << "," << step_time << "," << sys->GetTimerCollision() << "," << sys->GetTimerCollisionBroad() << "," << sys->GetTimerCollisionNarrow() << ","
               << sys->GetTimerUpdate() << "," << sys->GetTimerAdvance() << "," << sys->GetTimerProcessContact() << std::endl;

            static int printed = 0;
            if (printed++ % 10 == 0) {
                cout << std::fixed << std::setprecision(4) << "t = " << time << "  phase " << phase << "  strain = " << std::setprecision(3) << strain
                     << "  p_top = " << std::setprecision(1) << p_top << " Pa  mu_top = " << std::setprecision(3) << tau_top / p_top << "  phi = " << phi_center
                     << "  gap = " << std::setprecision(5) << H << "  contacts = " << st.num_contacts << "  step = " << std::setprecision(2) << 1e3 * step_time
                     << " ms  overflows = " << overflows_total << endl;
                cout.unsetf(std::ios::fixed);
            }
        }
    }

    double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall_start).count();
    write_meta(true, wall, overflows_total, "finished");
    cout << std::setprecision(8) << "Finished at t = " << sys->GetChTime() << " s, wall clock " << wall << " s, history overflows " << overflows_total << endl;

    return 0;
}

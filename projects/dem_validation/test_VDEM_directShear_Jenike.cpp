// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2026 projectchrono.org
// All right reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Author: Radu Serban
// =============================================================================
//
// Soft-sphere (SMC) direct shear test in a Jenike shear cell, with Chrono::Multicore: a replication of the Hartl and
// Ooi (2008) tests on 6 mm glass beads.
//
// The cell, after Hartl and Ooi: 143 mm inner diameter; a bottom ring 19 mm high, closed by a base, and a top ring
// 24 mm high, open at both ends, resting on the bottom ring (no gap). The round walls are rings of box plates. A lid
// sits on the sample inside the top ring. Lid and base are smooth, or grooved (--grooves): straight grooves across the
// shear direction with the sawtooth profile of Hartl and Ooi's Fig. 1 (90-degree teeth, as deep as half the pitch).
//
// The run has four phases:
// - fill  : beads rain into the cell in layers, from random positions above the top ring
// - level : the beads settle; at t_level those whose centers are above the rim are retired (fixed, without collision,
//           and parked above the cell)
// - load  : the lid is placed on the sample and released. Its weight plus the weight of the beads above the shear
//           plane (the top of the bottom ring) gives the normal stress on the shear plane. The lid moves with the top
//           ring horizontally; it is free vertically and, with --lid_tilt, free to tilt.
// - shear : from t_shear, the top ring and the lid move at constant speed in +x. The top ring starts offset by 3 mm in
//           -x (as in Hartl and Ooi's DEM), so the shear displacement runs from -3 mm of offset onwards.
//
// The shear stress is the x contact force of the beads on the top ring and the lid, over the cell area.
//
// The bead friction is mu_fill until the lid is placed, then mu. Rolling friction (--mu_roll, off by default) acts
// from the lid placement on, at bead-bead and bead-wall contacts: Chrono::Multicore's torque mu_roll R F_N against the
// relative rolling of the two surfaces, applied while the rolling speed exceeds --min_roll_vel.
//
// The global reference frame has Y up; the shear is along +x. All units SI. Run-time visualization with Chrono::VSG
// is enabled with --vsg.
//
// Outputs, in the output directory:
//   meta.json           all parameters, run information and summary
//   timeseries.csv      one row per output interval (written as the run goes)
//   particles/*.csv     particle snapshots of the beads in the cell (at the shear displacements in --snaps, and frames)
//   contacts/*.csv      contact snapshots, at the same times
//   state_fill.bin      the leveled sample before the lid is placed (with --save_fill)
//   state_shear.bin     the loaded sample at the start of shearing (with --save_state)
//   states/*.bin        the state at each snapshot or frame (with --snap_states or --frame_states)
// A run can start from a leveled sample (--load_fill: one filling for several normal stresses) or resume shearing
// from a loaded one (--load_state).
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
#include "chrono/physics/ChLinkLock.h"
#include "chrono/physics/ChLinkMate.h"
#include "chrono/utils/ChUtilsCreators.h"

#include "chrono_multicore/physics/ChSystemMulticore.h"

#include "chrono_thirdparty/cxxopts/ChCLI.h"

#ifdef CHRONO_VSG
    #include "chrono_vsg/ChVisualSystemVSG.h"
    #include "chrono_multicore/visualization/ChMulticoreVisualizationVSG.h"
#endif

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
    double radius;   // m
    double density;  // kg/m^3
    double E;        // Young's modulus (beads and walls), Pa
    double nu;       // Poisson ratio
    double cor;      // coefficient of restitution
    double mu;       // bead friction (glass on glass)
    double mu_wall;  // wall friction (glass on steel; the composite friction is the smaller of the two)
    double mu_fill;  // bead friction while filling and leveling
    double mu_roll;  // rolling friction of beads and walls, from the lid placement on

    // Cell
    double diameter;  // inner diameter, m
    double h_bottom;  // height of the bottom ring, m
    double h_top;     // height of the top ring, m
    double offset;    // initial offset of the top ring in -x, m
    double wall;      // wall thickness, m
    int plates;       // plates per ring
    double h_lid;     // lid thickness, m
    double lid_gap;   // radial clearance between lid and top ring, m
    bool lid_tilt;    // let the lid tilt
    double grooves;   // pitch of the sawtooth grooves on lid and base, m (0: smooth)

    // Filling
    int num_beads;       // beads rained in (0: 15% more than the cell holds at porosity 0.4)
    double rain_height;  // height of the release plane above the rim, m
    double rain_speed;   // initial downward speed, m/s
    double rain_dt;      // time between layers, s
    double rain_gap;     // smallest gap between beads of a layer, m
    int rain_layer;      // most beads in a layer (0: as many as fit)

    // Loading
    double sigma;               // normal stress on the shear plane, Pa
    double gravity;             // m/s^2
    double speed;               // shear speed, m/s
    double travel;              // shear travel, m
    double t_level;             // retire the beads above the rim and place the lid, s
    double t_shear;             // start of shearing, s
    std::string history;        // tangential displacement model: multistep, onestep, none
    std::string history_shear;  // model while shearing ("same": as --history)

    // Numerics
    double dt;
    double min_roll_vel;  // rolling speed below which rolling friction is off, m/s
    int threads;
    int seed;

    // Output
    std::string out_dir;
    double out_dt;              // time series interval, s
    std::vector<double> snaps;  // shear displacements of the snapshots, m
    double frame_disp;          // shear displacement between frames, m (0: none)
    std::string save_fill;      // file for the leveled sample
    std::string load_fill;      // start from a leveled sample
    bool fill_only;             // stop after saving the leveled sample
    std::string save_state;     // file for the state at the start of shearing
    std::string load_state;     // resume shearing from a saved state
    bool snap_states;           // also save the state at each snapshot
    bool frame_states;          // also save the state at each frame
    bool vsg;
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
    ChCLI cli(argv[0], "Direct shear test in a Jenike cell (Chrono::Multicore, SMC)");

    cli.AddOption<double>("Beads", "radius", "Bead radius (m)", "0.003");
    cli.AddOption<double>("Beads", "density", "Bead density (kg/m3)", "2550");
    cli.AddOption<double>("Beads", "E", "Young's modulus of beads and walls (Pa)", "4e7");
    cli.AddOption<double>("Beads", "nu", "Poisson ratio", "0.22");
    cli.AddOption<double>("Beads", "cor", "Coefficient of restitution", "0.87");
    cli.AddOption<double>("Beads", "mu", "Bead friction coefficient (glass on glass)", "0.18");
    cli.AddOption<double>("Beads", "mu_wall", "Wall friction coefficient (glass on steel)", "0.13");
    cli.AddOption<double>("Beads", "mu_fill", "Bead friction while filling and leveling (default: mu)");
    cli.AddOption<double>("Beads", "mu_roll", "Rolling friction coefficient of beads and walls, from the lid placement on", "0");

    cli.AddOption<double>("Cell", "diameter", "Inner diameter (m)", "0.143");
    cli.AddOption<double>("Cell", "h_bottom", "Height of the bottom ring (m)", "0.019");
    cli.AddOption<double>("Cell", "h_top", "Height of the top ring (m)", "0.024");
    cli.AddOption<double>("Cell", "offset", "Initial offset of the top ring in -x (m)", "0.003");
    cli.AddOption<double>("Cell", "wall", "Wall thickness (m)", "0.005");
    cli.AddOption<int>("Cell", "plates", "Plates per ring", "180");
    cli.AddOption<double>("Cell", "h_lid", "Lid thickness (m)", "0.01");
    cli.AddOption<double>("Cell", "lid_gap", "Radial clearance between lid and top ring (m)", "0.0005");
    cli.AddOption<bool>("Cell", "lid_tilt", "Let the lid tilt (default: the lid only moves vertically in the top ring)");
    cli.AddOption<double>("Cell", "grooves", "Pitch of the sawtooth grooves on lid and base (mm; Hartl and Ooi: 1; 0: smooth)", "0");

    cli.AddOption<int>("Filling", "num_beads", "Beads rained in (0: 15% more than the cell holds at porosity 0.4)", "0");
    cli.AddOption<double>("Filling", "rain_height", "Height of the release plane above the rim (m)", "0.02");
    cli.AddOption<double>("Filling", "rain_speed", "Initial downward speed of the beads (m/s)", "0.5");
    cli.AddOption<double>("Filling", "rain_dt", "Time between layers (s)", "0.01");
    cli.AddOption<double>("Filling", "rain_gap", "Smallest gap between the beads of a layer (m)", "0.001");
    cli.AddOption<int>("Filling", "rain_layer", "Most beads in a layer (0: as many as fit)", "40");

    cli.AddOption<double>("Loading", "sigma", "Normal stress on the shear plane (Pa)", "3100");
    cli.AddOption<double>("Loading", "gravity", "Gravitational acceleration (m/s2)", "9.81");
    cli.AddOption<double>("Loading", "speed", "Shear speed (m/s)", "0.001");
    cli.AddOption<double>("Loading", "travel", "Shear travel (m)", "0.012");
    cli.AddOption<double>("Loading", "t_level", "Retire the beads above the rim and place the lid (s)", "1.5");
    cli.AddOption<double>("Loading", "t_shear", "Start of shearing (s)", "2.5");
    cli.AddOption<std::string>("Loading", "history", "Tangential displacement model: multistep, onestep or none", "multistep");
    cli.AddOption<std::string>("Loading", "history_shear", "Tangential displacement model while shearing (same: as --history)", "same");

    cli.AddOption<double>("Numerics", "dt", "Time step (s)", "1e-5");
    cli.AddOption<double>("Numerics", "min_roll_vel", "Rolling speed below which rolling friction is off (m/s)", "1e-6");
    cli.AddOption<int>("Numerics", "threads", "Number of OpenMP threads", "20");
    cli.AddOption<int>("Numerics", "seed", "Seed of the random bead positions", "1");

    cli.AddOption<std::string>("Output", "out", "Output directory", "directShear_Jenike");
    cli.AddOption<double>("Output", "out_dt", "Time series interval (s)", "0.01");
    cli.AddOption<std::vector<double>>("Output", "snaps", "Shear displacements of the snapshots (mm)", "0,0.5,1,2,4,8,12");
    cli.AddOption<double>("Output", "frame_disp", "Shear displacement between particle and contact frames (mm, 0: none)", "0");
    cli.AddOption<std::string>("Output", "save_fill", "Save the leveled sample (before the lid is placed) to this file", "");
    cli.AddOption<std::string>("Output", "load_fill", "Start from a leveled sample saved with --save_fill", "");
    cli.AddOption<bool>("Output", "fill_only", "Stop after the leveled sample is saved");
    cli.AddOption<std::string>("Output", "save_state", "Save the sample at the start of shearing to this file", "");
    cli.AddOption<std::string>("Output", "load_state", "Resume shearing from a state saved with --save_state or --snap_states", "");
    cli.AddOption<bool>("Output", "snap_states", "Also save the state at each snapshot (states/*.bin)");
    cli.AddOption<bool>("Output", "frame_states", "Also save the state at each frame (states/*.bin)");
    cli.AddOption<bool>("Output", "vsg", "Run-time visualization with Chrono::VSG");
    cli.AddOption<std::string>("Output", "provenance", "JSON text recorded in meta.json", "");

    if (!cli.Parse(argc, argv, true))
        return false;

    p.radius = cli.GetAsType<double>("radius");
    p.density = cli.GetAsType<double>("density");
    p.E = cli.GetAsType<double>("E");
    p.nu = cli.GetAsType<double>("nu");
    p.cor = cli.GetAsType<double>("cor");
    p.mu = cli.GetAsType<double>("mu");
    p.mu_wall = cli.GetAsType<double>("mu_wall");
    p.mu_fill = cli.CheckOption("mu_fill") ? cli.GetAsType<double>("mu_fill") : p.mu;
    p.mu_roll = cli.GetAsType<double>("mu_roll");

    p.diameter = cli.GetAsType<double>("diameter");
    p.h_bottom = cli.GetAsType<double>("h_bottom");
    p.h_top = cli.GetAsType<double>("h_top");
    p.offset = cli.GetAsType<double>("offset");
    p.wall = cli.GetAsType<double>("wall");
    p.plates = cli.GetAsType<int>("plates");
    p.h_lid = cli.GetAsType<double>("h_lid");
    p.lid_gap = cli.GetAsType<double>("lid_gap");
    p.lid_tilt = cli.GetAsType<bool>("lid_tilt");
    p.grooves = 1e-3 * cli.GetAsType<double>("grooves");

    p.num_beads = cli.GetAsType<int>("num_beads");
    p.rain_height = cli.GetAsType<double>("rain_height");
    p.rain_speed = cli.GetAsType<double>("rain_speed");
    p.rain_dt = cli.GetAsType<double>("rain_dt");
    p.rain_gap = cli.GetAsType<double>("rain_gap");
    p.rain_layer = cli.GetAsType<int>("rain_layer");

    p.sigma = cli.GetAsType<double>("sigma");
    p.gravity = cli.GetAsType<double>("gravity");
    p.speed = cli.GetAsType<double>("speed");
    p.travel = cli.GetAsType<double>("travel");
    p.t_level = cli.GetAsType<double>("t_level");
    p.t_shear = cli.GetAsType<double>("t_shear");
    p.history = cli.GetAsType<std::string>("history");
    p.history_shear = cli.GetAsType<std::string>("history_shear");
    if (p.history_shear == "same")
        p.history_shear = p.history;
    HistoryModel(p.history);
    HistoryModel(p.history_shear);

    p.dt = cli.GetAsType<double>("dt");
    p.min_roll_vel = cli.GetAsType<double>("min_roll_vel");
    p.threads = cli.GetAsType<int>("threads");
    p.seed = cli.GetAsType<int>("seed");

    p.out_dir = cli.GetAsType<std::string>("out");
    p.out_dt = cli.GetAsType<double>("out_dt");
    p.snaps = cli.GetAsType<std::vector<double>>("snaps");
    for (auto& s : p.snaps)
        s *= 1e-3;
    p.frame_disp = 1e-3 * cli.GetAsType<double>("frame_disp");
    p.save_fill = cli.GetAsType<std::string>("save_fill");
    p.load_fill = cli.GetAsType<std::string>("load_fill");
    p.fill_only = cli.GetAsType<bool>("fill_only");
    p.save_state = cli.GetAsType<std::string>("save_state");
    p.load_state = cli.GetAsType<std::string>("load_state");
    p.snap_states = cli.GetAsType<bool>("snap_states");
    p.frame_states = cli.GetAsType<bool>("frame_states");
    p.vsg = cli.GetAsType<bool>("vsg");
    p.provenance = cli.GetAsType<std::string>("provenance");

    if (p.fill_only && p.save_fill.empty()) {
        std::cerr << "--fill_only needs --save_fill" << endl;
        return false;
    }

    return true;
}

// -----------------------------------------------------------------------------
// Sawtooth grooves across the x direction on a disk of the given radius, centered under (or over) the body, on the
// plane y = y_base (global frame). Each tooth is a square box rotated 45 degrees about the groove axis and centered on
// the plane: its half on the sample side is a 90-degree tooth as deep as half the pitch; the other half is inside the
// base or the lid.
// -----------------------------------------------------------------------------

void AddGrooves(std::shared_ptr<ChBody> body, std::shared_ptr<ChContactMaterial> mat, double pitch, double radius, double y_base) {
    const auto& X = body->GetFrameRefToAbs();
    ChQuaterniond tooth_rot = QuatFromAngleZ(CH_PI / 4);
    double side = pitch / std::sqrt(2.0);
    int n = (int)std::floor(2 * radius / pitch);
    for (int k = 0; k < n; k++) {
        double x = -radius + pitch * (k + 0.5);
        double x_out = std::abs(x) + 0.5 * pitch;
        if (x_out >= radius)
            continue;
        double half_chord = std::sqrt(radius * radius - x_out * x_out);
        ChVector3d pos(X.GetPos().x() + x, y_base, X.GetPos().z());
        utils::AddBoxGeometry(body.get(), mat, ChVector3d(side, side, 2 * half_chord), X.TransformPointParentToLocal(pos), X.GetRot().GetConjugate() * tooth_rot, true);
    }
}

// -----------------------------------------------------------------------------

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
    if (p.snap_states || p.frame_states)
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

    double r = p.radius;
    double R = p.diameter / 2;
    double area = CH_PI * R * R;
    double bead_volume = (4.0 / 3) * CH_PI * r * r * r;
    double bead_mass = p.density * bead_volume;
    double y_joint = p.h_bottom;          // shear plane: top of the bottom ring (base top at y = 0)
    double y_rim = p.h_bottom + p.h_top;  // top of the top ring
    double end_time = p.t_shear + p.travel / p.speed;

    int num_beads = p.num_beads;
    if (num_beads <= 0)
        num_beads = (int)std::ceil(1.15 * 0.6 * area * y_rim / bead_volume);

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
    sys->GetSettings()->solver.min_roll_vel = p.min_roll_vel;

    sys->GetSettings()->collision.bins_per_axis = vec3(10, 10, 10);
    sys->GetSettings()->collision.narrowphase_algorithm = ChNarrowphase::Algorithm::HYBRID;

    // Materials: beads (glass), and the cell walls, base and lid (steel)

    auto material = chrono_types::make_shared<ChContactMaterialSMC>();
    material->SetYoungModulus((float)p.E);
    material->SetPoissonRatio((float)p.nu);
    material->SetRestitution((float)p.cor);
    material->SetFriction((float)p.mu_fill);

    auto mat_ext = chrono_types::make_shared<ChContactMaterialSMC>();
    mat_ext->SetYoungModulus((float)p.E);
    mat_ext->SetPoissonRatio((float)p.nu);
    mat_ext->SetRestitution((float)p.cor);
    mat_ext->SetFriction((float)p.mu_wall);

    // The cell walls and the lid are created before the beads, so that the history of the wall contacts is stored on
    // the beads (the body with the larger index)

    // Bottom ring with the base (body 0); base top at y = 0. The plates in one quadrant (x > 0, z < 0) are not drawn.
    auto bottom = utils::CreateCylindricalContainerFromBoxes(sys, mat_ext, R, p.h_bottom, p.wall, p.plates, ChVector3d(0, 0, 0), z2y, true, false, false, true, true, true);
    if (p.grooves > 0)
        AddGrooves(bottom, mat_ext, p.grooves, R, 0);
    bottom->GetCollisionModel()->SetFamily(1);
    bottom->GetCollisionModel()->DisallowCollisionsWith(2);
    bottom->GetCollisionModel()->DisallowCollisionsWith(3);

    // Top ring, open at both ends (body 1), offset in -x
    double x_top0 = -p.offset;
    auto top = utils::CreateCylindricalContainerFromBoxes(sys, mat_ext, R, p.h_top, p.wall, p.plates, ChVector3d(x_top0, y_joint, 0), z2y, true, false, false, true, true, false);
    top->GetCollisionModel()->SetFamily(2);
    top->GetCollisionModel()->DisallowCollisionsWith(1);
    top->GetCollisionModel()->DisallowCollisionsWith(3);

    // Lid (body 2), parked above the rain until it is placed; mass set when it is placed
    double y_park_lid = y_rim + p.rain_height + 0.1;
    auto lid = chrono_types::make_shared<ChBody>();
    lid->SetMass(1);
    lid->SetPos(ChVector3d(x_top0, y_park_lid, 0));
    lid->SetFixed(true);
    lid->EnableCollision(true);
    utils::AddCylinderGeometry(lid.get(), mat_ext, R - p.lid_gap, p.h_lid, ChVector3d(0, 0, 0), z2y, true);
    if (p.grooves > 0)
        AddGrooves(lid, mat_ext, p.grooves, R - p.lid_gap, y_park_lid - p.h_lid / 2);
    lid->GetCollisionModel()->SetFamily(3);
    lid->GetCollisionModel()->DisallowCollisionsWith(1);
    lid->GetCollisionModel()->DisallowCollisionsWith(2);
    sys->AddBody(lid);

    // The lid moves with the top ring horizontally and is free vertically (with --lid_tilt, also free to tilt)
    if (p.lid_tilt) {
        auto joint = chrono_types::make_shared<ChLinkMateGeneric>(true, false, true, false, true, false);
        joint->Initialize(lid, top, ChFrame<>(lid->GetPos(), QUNIT));
        sys->AddLink(joint);
    } else {
        auto joint = chrono_types::make_shared<ChLinkLockPrismatic>();
        joint->Initialize(lid, top, ChFrame<>(lid->GetPos(), z2y));
        sys->AddLink(joint);
    }

    // Beads, all created now, fixed and without collision, and released into the rain layer by layer (unless a saved
    // state gives their positions). A bead cannot leave the Multicore collision system once it is in it, so a retired
    // bead keeps its collision shape: it is fixed and parked on a grid above the cell, apart from the other beads.

    double y_park = y_rim + p.rain_height + 0.05;
    double park_pitch = 2 * r + 0.002;
    int park_row = (int)std::floor(2 * R / park_pitch);
    int num_parked = 0;
    std::vector<std::shared_ptr<ChBody>> beads;
    int first_particle = (int)sys->GetBodies().size();

    auto retire = [&](std::shared_ptr<ChBody> b) {
        int k = num_parked++;
        int i = k % park_row;
        int j = (k / park_row) % park_row;
        int l = k / (park_row * park_row);
        b->SetPos(ChVector3d(-R + park_pitch * (i + 0.5), y_park + park_pitch * l, -R + park_pitch * (j + 0.5)));
        b->SetPosDt(VNULL);
        b->SetAngVelParent(VNULL);
        b->SetFixed(true);
    };
    auto release = [&](std::shared_ptr<ChBody> b, const ChVector3d& pos) {
        b->SetPos(pos);
        b->SetPosDt(ChVector3d(0, -p.rain_speed, 0));
        b->SetFixed(false);
        b->EnableCollision(true);
    };
    auto is_retired = [&](const std::shared_ptr<ChBody>& b) { return b->GetPos().y() > y_park - 0.5 * park_pitch; };

    for (int i = 0; i < num_beads; i++) {
        auto bead = chrono_types::make_shared<ChBody>();
        bead->SetMass(bead_mass);
        bead->SetInertiaXX(0.4 * bead_mass * r * r * ChVector3d(1, 1, 1));
        auto ct_shape = chrono_types::make_shared<ChCollisionShapeSphere>(material, r);
        bead->AddCollisionShape(ct_shape);
        bead->GetCollisionModel()->SetFamily(4);
        auto sphere = chrono_types::make_shared<ChVisualShapeSphere>(r);
        sphere->SetColor(ChColor(0.8f, 0.8f, 0.8f));
        bead->AddVisualShape(sphere);
        bead->SetPos(ChVector3d(0, y_park + 1, 0));
        bead->SetFixed(true);
        bead->EnableCollision(false);
        sys->AddBody(bead);
        beads.push_back(bead);
    }

    // Rain: layers of beads at random positions in a disk above the top ring, apart from one another and from beads
    // that bounced up to the release plane
    std::mt19937 rng((unsigned int)p.seed);
    std::uniform_real_distribution<double> uni(-1.0, 1.0);
    int next_bead = 0;
    int next_layer = 0;
    double rain_radius = R - r - p.rain_gap;
    double y_rain = y_rim + p.rain_height;
    auto rain_layer = [&]() {
        double dmin = 2 * r + p.rain_gap;
        std::vector<ChVector3d> taken;
        for (int i = 0; i < next_bead; i++) {
            const auto& q = beads[i]->GetPos();
            if (!beads[i]->IsFixed() && std::abs(q.y() - y_rain) < dmin)
                taken.push_back(q);
        }
        size_t num_taken = taken.size();
        int layer_max = num_beads - next_bead;
        if (p.rain_layer > 0)
            layer_max = std::min(layer_max, p.rain_layer);
        for (int attempt = 0; attempt < 20000 && (int)(taken.size() - num_taken) < layer_max; attempt++) {
            double x = rain_radius * uni(rng);
            double z = rain_radius * uni(rng);
            if (x * x + z * z > rain_radius * rain_radius)
                continue;
            ChVector3d pos(x_top0 + x, y_rain, z);
            bool ok = true;
            for (const auto& q : taken) {
                if ((q - pos).Length2() < dmin * dmin) {
                    ok = false;
                    break;
                }
            }
            if (ok)
                taken.push_back(pos);
        }
        for (size_t k = num_taken; k < taken.size(); k++)
            release(beads[next_bead++], taken[k]);
        next_layer++;
    };

    // Phase flags and reference values
    // phase 0: fill, 1: level and settle, 2: load (lid released), 3: shear

    int phase = 0;
    double lid_mass = 0;      // set when the lid is placed
    double top_mass = 0;      // mass of the beads above the shear plane when the lid is placed
    int num_active = 0;       // beads in the cell after leveling
    double solid_volume = 0;  // their volume
    double shear_Height = 0;  // lid height at the start of shearing
    double fill_time = NAN;   // time at which the last bead was released

    auto place_lid = [&]() {
        // Weight of the beads above the shear plane, and the lid mass that completes the normal stress
        top_mass = 0;
        double y_top = -1e10;
        for (const auto& b : beads) {
            if (is_retired(b))
                continue;
            if (b->GetPos().y() > y_joint)
                top_mass += b->GetMass();
            y_top = std::max(y_top, b->GetPos().y() + r);
        }
        lid_mass = p.sigma * area / p.gravity - top_mass;
        if (lid_mass <= 0) {
            std::cerr << "The beads above the shear plane alone exceed the normal stress; using a lid of 1 g" << endl;
            lid_mass = 1e-3;
        }
        double y_lid = std::max(y_top, y_rim) + p.h_lid / 2 + 0.0005;
        if (p.grooves > 0)
            y_lid += 0.5 * p.grooves;
        lid->SetMass(lid_mass);
        lid->SetInertiaXX(lid_mass * ChVector3d(0.25 * R * R, 0.5 * R * R, 0.25 * R * R));
        lid->SetPos(ChVector3d(top->GetPos().x(), y_lid, 0));
        lid->SetFixed(false);
        material->SetFriction((float)p.mu);
        material->SetRollingFriction((float)p.mu_roll);
        mat_ext->SetRollingFriction((float)p.mu_roll);
    };

    // Start from a leveled sample, or resume shearing from a loaded one
    if (!p.load_fill.empty() || !p.load_state.empty()) {
        const std::string& file = p.load_state.empty() ? p.load_fill : p.load_state;
        std::vector<double> extra;
        if (!LoadState(*sys, file, extra) || extra.size() != 4) {
            std::cerr << "Cannot load state from " << file << " (missing file or different cell or bead count)" << endl;
            return 1;
        }
        next_bead = num_beads;
        fill_time = extra[0];
        num_active = 0;
        for (auto& b : beads) {
            if (!is_retired(b)) {
                b->SetFixed(false);
                b->EnableCollision(true);
                num_active++;
            }
        }
        num_parked = num_beads - num_active;
        solid_volume = num_active * bead_volume;
        if (p.load_state.empty()) {
            phase = 1;
            cout << "Loaded leveled sample at t = " << sys->GetChTime() << " s from " << file << " (" << num_active << " beads)" << endl;
        } else {
            phase = 3;
            lid_mass = extra[1];
            top_mass = extra[2];
            shear_Height = extra[3];
            lid->SetMass(lid_mass);
            lid->SetInertiaXX(lid_mass * ChVector3d(0.25 * R * R, 0.5 * R * R, 0.25 * R * R));
            lid->SetFixed(false);
            material->SetFriction((float)p.mu);
            material->SetRollingFriction((float)p.mu_roll);
            mat_ext->SetRollingFriction((float)p.mu_roll);
            sys->GetSettings()->solver.tangential_displ_mode = HistoryModel(p.history_shear);
            cout << "Loaded state at t = " << sys->GetChTime() << " s from " << file << " (" << num_active << " beads)" << endl;
        }
    }

    // Run metadata (rewritten at leveling and at the end)

    auto write_meta = [&](bool finished, double wall_time, unsigned long long overflows, const std::string& status) {
        JsonObject m;
        m.Add("program", "test_VDEM_directShear_Jenike");
        m.Add("description", "Direct shear in a Jenike cell, Chrono::Multicore SMC (Hertz), Y up, shear along +x, SI units");
        m.Add("chrono_version", CHRONO_VERSION);
        m.Add("start", start_date);
        m.Add("status", status);
        m.Add("radius", p.radius);
        m.Add("density", p.density);
        m.Add("bead_mass", bead_mass);
        m.Add("E", p.E);
        m.Add("nu", p.nu);
        m.Add("cor", p.cor);
        m.Add("mu_bead", p.mu);
        m.Add("mu_wall", p.mu_wall);
        m.Add("mu_bead_wall_composite", std::min(p.mu, p.mu_wall));
        m.Add("mu_fill", p.mu_fill);
        m.Add("mu_roll", p.mu_roll);
        m.Add("min_roll_vel", p.min_roll_vel);
        m.Add("cell_diameter", p.diameter);
        m.Add("h_bottom", p.h_bottom);
        m.Add("h_top", p.h_top);
        m.Add("ring_gap", 0.0);
        m.Add("top_ring_offset", p.offset);
        m.Add("wall_thickness", p.wall);
        m.Add("plates_per_ring", p.plates);
        m.Add("h_lid", p.h_lid);
        m.Add("lid_gap", p.lid_gap);
        m.Add("lid_tilt", p.lid_tilt);
        m.Add("grooves_pitch", p.grooves);
        m.Add("nominal_area", area);
        m.Add("num_beads_rained", num_beads);
        m.Add("num_beads_in_cell", num_active);
        m.Add("rain_height", p.rain_height);
        m.Add("rain_speed", p.rain_speed);
        m.Add("rain_dt", p.rain_dt);
        m.Add("rain_gap", p.rain_gap);
        m.Add("rain_layer", p.rain_layer);
        m.Add("rain_layers", next_layer);
        m.Add("fill_time", fill_time);
        m.Add("sigma", p.sigma);
        m.Add("lid_mass", lid_mass);
        m.Add("top_ring_bead_mass", top_mass);
        m.Add("gravity", p.gravity);
        m.Add("shear_speed", p.speed);
        m.Add("shear_travel", p.travel);
        m.Add("t_level", p.t_level);
        m.Add("t_shear", p.t_shear);
        m.Add("t_end", end_time);
        m.Add("history_prepare", p.history);
        m.Add("history_shear", p.history_shear);
        m.Add("contact_force_model", "Hertz");
        m.Add("narrowphase", "HYBRID");
        m.Add("dt", p.dt);
        m.Add("threads", threads);
        m.Add("seed", p.seed);
        m.Add("out_dt", p.out_dt);
        m.Add("snapshots_shear_disp", p.snaps);
        m.Add("frame_disp", p.frame_disp);
        m.Add("save_fill", p.save_fill);
        m.Add("load_fill", p.load_fill);
        m.Add("save_state", p.save_state);
        m.Add("load_state", p.load_state);
        m.Add("snap_states", p.snap_states);
        m.Add("frame_states", p.frame_states);
        m.Add("first_particle_index", first_particle);
        m.AddRaw("walls", "{\"bottom ring and base\": {\"body\": 0}, \"top ring\": {\"body\": 1}, \"lid\": {\"body\": 2}}");
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
    ts << "# Direct shear, Jenike cell. SI units. Y up; the top ring and the lid move in +x.\n";
    ts << "# phase: 0 fill, 1 level, 2 load, 3 shear. Forces are contact forces of the beads on the named body.\n";
    ts << "# tau = -(Fx(top ring) + Fx(lid))/A; ratio = tau/sigma; tau_bottom = Fx(bottom ring and base)/A.\n";
    ts << "# e_lid: void ratio from the lid height over the cell area; n_active: beads in the cell.\n";
    ts << "time,phase,shear_disp,lid_y,dilation,lid_tilt,"
          "top_fx,top_fy,top_fz,lid_fx,lid_fy,lid_fz,bottom_fx,bottom_fy,bottom_fz,"
          "tau,ratio,tau_bottom,e_lid,n_active,kinetic_energy,"
          "contacts,contacts_pp,contacts_pw,coord,coord_mech,rattlers,frac_sliding,mean_fn,"
          "overflows_step,overflows_total,wall_clock,step_time,"
          "t_collision,t_broad,t_narrow,t_update,t_advance,t_contacts,t_ls_setup,t_ls_solve\n";
    ts << std::setprecision(9);

    std::ofstream frames_index(out_dir + "/frames.csv");
    frames_index << "frame,time,shear_disp,phase\n";

#ifdef CHRONO_VSG
    std::shared_ptr<vsg3d::ChVisualSystemVSG> vis;
    if (p.vsg) {
        auto mcore_vis = chrono_types::make_shared<ChMulticoreVisualizationVSG>(sys);
        mcore_vis->SetNetworkColorRamp(ChColor(0.957f, 0.694f, 0.514f), ChColor(0.518f, 0.235f, 0.047f));
        mcore_vis->SetBodySelection(beads);  // the beads are parked and inactive at initialization

        vis = chrono_types::make_shared<vsg3d::ChVisualSystemVSG>();
        vis->AttachSystem(sys);
        vis->AttachPlugin(mcore_vis);
        vis->SetWindowTitle("Direct shear, Jenike cell");
        vis->SetWindowSize(1280, 720);
        vis->SetCameraVertical(CameraVerticalDir::Y);
        vis->AddCamera(ChVector3d(0, 2.5 * y_rim, -3 * p.diameter), ChVector3d(0, 0.5 * y_rim, 0));
        vis->SetBackgroundColor(ChColor(1.0f, 1.0f, 1.0f));
        vis->SetLightIntensity(1.0f);
        vis->Initialize();
    }
#else
    if (p.vsg)
        cout << "Chrono::VSG not available; running without visualization" << endl;
#endif

    // Simulation loop

    int data_out_frame = (int)std::floor(sys->GetChTime() / p.out_dt + 1e-9);
    int frame = 0;
    size_t next_snap = 0;
    double last_wall = 0;
    long long last_step = 0;
    long long step = 0;
    unsigned long long overflows_total = 0;

    auto shear_disp = [&]() { return phase == 3 ? top->GetPos().x() - x_top0 : 0.0; };

    // A run resumed from a saved state skips the snapshots and frames before the loaded shear displacement
    if (phase == 3) {
        double d0 = shear_disp() - 0.5 * p.speed * p.dt;
        while (next_snap < p.snaps.size() && p.snaps[next_snap] < d0)
            next_snap++;
        if (p.frame_disp > 0)
            frame = (int)std::ceil(d0 / p.frame_disp);
    }

    auto active_beads = [&]() {
        std::vector<std::shared_ptr<ChBody>> list;
        for (const auto& b : beads)
            if (!is_retired(b))
                list.push_back(b);
        return list;
    };

    auto write_snapshot = [&](const std::string& name) {
        auto contacts = GetContacts(*sys);
        std::vector<int> counts;
        ComputeContactStats(contacts, first_particle, num_beads, &counts);
        std::vector<std::shared_ptr<ChBody>> list;
        std::vector<int> list_counts;
        for (int i = 0; i < num_beads; i++) {
            if (!is_retired(beads[i])) {
                list.push_back(beads[i]);
                list_counts.push_back(counts[i]);
            }
        }
        std::vector<double> list_radii(list.size(), r);
        WriteParticles(particles_dir + "/" + name + ".csv", list, list_radii, list_counts, sys->GetChTime(), shear_disp());
        WriteContacts(contacts_dir + "/" + name + ".csv", contacts, sys->GetChTime(), shear_disp());
    };

    auto save = [&](const std::string& file) { return SaveState(*sys, file, {fill_time, lid_mass, top_mass, shear_Height}); };

    while (sys->GetChTime() < end_time - 1e-9) {
        double time = sys->GetChTime();

        // Fill: one layer every rain_dt until all beads are released
        if (phase == 0) {
            if (next_bead < num_beads && time >= next_layer * p.rain_dt - 1e-9) {
                rain_layer();
                if (next_bead == num_beads) {
                    fill_time = time;
                    phase = 1;
                    if (fill_time > 0.7 * p.t_level)
                        std::cerr << "Warning: the last layer was released at t = " << fill_time << " s, close to t_level" << endl;
                }
            }
        }

        // Beads that bounced out of the cell while filling are retired once they fall below the base
        if (phase <= 1) {
            for (int i = 0; i < next_bead; i++) {
                if (!beads[i]->IsFixed() && beads[i]->GetPos().y() < -p.wall - 2 * r)
                    retire(beads[i]);
            }
        }

        // Level: retire the beads above the rim, save the leveled sample, then place and release the lid
        if (phase <= 1 && time >= p.t_level - 1e-9) {
            if (phase == 0)
                std::cerr << "Warning: leveling at t = " << time << " s with " << num_beads - next_bead << " beads not released" << endl;
            if (p.load_fill.empty()) {
                for (auto& b : beads) {
                    if (b->IsFixed())
                        continue;
                    const auto& q = b->GetPos();
                    double r_bottom = std::sqrt(q.x() * q.x() + q.z() * q.z());
                    double r_top = std::sqrt((q.x() - x_top0) * (q.x() - x_top0) + q.z() * q.z());
                    bool outside = r_bottom > R + 0.5 * r && r_top > R + 0.5 * r;
                    if (q.y() > y_rim || outside)
                        retire(b);
                }
            }
            num_active = (int)active_beads().size();
            solid_volume = num_active * bead_volume;
            if (!p.save_fill.empty() && p.load_fill.empty()) {
                if (save(p.save_fill))
                    cout << "Saved leveled sample at t = " << time << " s to " << p.save_fill << " (" << num_active << " beads)" << endl;
                else
                    std::cerr << "Cannot save state to " << p.save_fill << endl;
            }
            if (p.fill_only) {
                double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall_start).count();
                write_meta(true, wall, overflows_total, "filled");
                cout << "Stopped after filling (" << num_active << " beads in the cell)" << endl;
                return 0;
            }
            place_lid();
            phase = 2;
            write_meta(false, 0, overflows_total, "running");
        }

        // Shear: the top ring (and the lid with it) moves at constant speed
        if (phase == 2 && time > p.t_shear) {
            material->SetFriction((float)p.mu);
            shear_Height = lid->GetPos().y();
            phase = 3;
            if (!p.save_state.empty() && p.load_state.empty()) {
                if (save(p.save_state))
                    cout << "Saved state at t = " << time << " s to " << p.save_state << endl;
                else
                    std::cerr << "Cannot save state to " << p.save_state << endl;
            }
            sys->GetSettings()->solver.tangential_displ_mode = HistoryModel(p.history_shear);
        }

        if (phase == 3) {
            top->SetPos(ChVector3d(x_top0 + p.speed * (time - p.t_shear), y_joint, 0));
            top->SetLinVel(ChVector3d(p.speed, 0, 0));
        } else {
            top->SetPos(ChVector3d(x_top0, y_joint, 0));
            top->SetLinVel(VNULL);
        }
        top->SetRot(z2y);

        sys->DoStepDynamics(p.dt);
        step++;
        overflows_total = sys->GetTotalContactHistoryOverflows();

#ifdef CHRONO_VSG
        if (vis && step % 200 == 0) {
            if (!vis->Run())
                break;
            vis->Render();
        }
#endif

        time = sys->GetChTime();

        // Snapshots at the requested shear displacements, and frames
        if (phase == 3) {
            double d = shear_disp();
            while (next_snap < p.snaps.size() && d >= p.snaps[next_snap] - 0.5 * p.speed * p.dt) {
                char name[64];
                snprintf(name, sizeof(name), "snap_d%05.2fmm", 1e3 * p.snaps[next_snap]);
                write_snapshot(name);
                if (p.snap_states && !save(states_dir + "/" + name + ".bin"))
                    std::cerr << "Cannot save state to " << states_dir + "/" + name + ".bin" << endl;
                next_snap++;
            }
            if (p.frame_disp > 0 && d >= frame * p.frame_disp - 0.5 * p.speed * p.dt) {
                char name[64];
                snprintf(name, sizeof(name), "frame_%05d", frame);
                write_snapshot(name);
                if (p.frame_states && !save(states_dir + "/" + name + ".bin"))
                    std::cerr << "Cannot save state to " << states_dir + "/" + name + ".bin" << endl;
                frames_index << frame << "," << time << "," << d << ",3" << std::endl;
                frame++;
            }
        }

        // Time series
        if (time >= data_out_frame * p.out_dt) {
            auto contacts = GetContacts(*sys);
            auto st = ComputeContactStats(contacts, first_particle, num_beads);

            // Retired and unreleased beads (fixed) have no contacts: count them out of the coordination numbers
            int n_in = 0;
            for (const auto& b : beads)
                n_in += b->IsFixed() ? 0 : 1;
            st.num_rattlers -= num_beads - n_in;
            st.coord = n_in > 0 ? st.coord * num_beads / n_in : 0;

            // Contact forces on the cell parts
            ChVector3d f_bottom(0), f_top(0), f_lid(0);
            for (const auto& c : contacts) {
                // Force on body 2 is F; on body 1 it is -F
                ChVector3d F = c.fn * c.n + c.ft;
                auto add = [&](int b, const ChVector3d& f) {
                    if (b == 0)
                        f_bottom += f;
                    else if (b == 1)
                        f_top += f;
                    else if (b == 2)
                        f_lid += f;
                };
                add(c.b1, -F);
                add(c.b2, F);
            }

            double lid_y = lid->GetPos().y();
            double dilation = phase == 3 ? lid_y - shear_Height : 0.0;
            double tilt = std::acos(std::min(1.0, lid->GetRotMat().GetAxisY().y()));
            double tau = -(f_top.x() + f_lid.x()) / area;
            double tau_bottom = f_bottom.x() / area;

            // Void ratio from the lid height (lid bottom to base top) over the cell area
            // (with grooves, the teeth on the lid and on the base each fill a layer a quarter of the pitch thick)
            double lid_bottom = lid_y - p.h_lid / 2 - 0.5 * p.grooves;
            double e_lid = (phase >= 2 && solid_volume > 0) ? (area * lid_bottom - solid_volume) / solid_volume : NAN;

            double ke = 0;
            for (const auto& b : beads) {
                if (b->IsFixed())
                    continue;
                ke += 0.5 * b->GetMass() * b->GetPosDt().Length2();
                ke += 0.5 * 0.4 * b->GetMass() * r * r * b->GetAngVelParent().Length2();
            }

            double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall_start).count();
            double step_time = (step > last_step) ? (wall - last_wall) / (step - last_step) : 0;
            last_wall = wall;
            last_step = step;

            ts << time << "," << phase << "," << shear_disp() << "," << lid_y << "," << dilation << "," << tilt << "," << f_top.x() << "," << f_top.y() << "," << f_top.z() << ","
               << f_lid.x() << "," << f_lid.y() << "," << f_lid.z() << "," << f_bottom.x() << "," << f_bottom.y() << "," << f_bottom.z() << "," << tau << "," << tau / p.sigma
               << "," << tau_bottom << "," << e_lid << "," << n_in << "," << ke << "," << st.num_contacts << "," << st.num_pp << "," << st.num_pw << "," << st.coord << ","
               << st.coord_mech << "," << st.num_rattlers << "," << st.frac_sliding << "," << st.mean_fn << "," << sys->GetNumContactHistoryOverflows() << "," << overflows_total
               << "," << wall << "," << step_time << "," << sys->GetTimerCollision() << "," << sys->GetTimerCollisionBroad() << "," << sys->GetTimerCollisionNarrow() << ","
               << sys->GetTimerUpdate() << "," << sys->GetTimerAdvance() << "," << sys->GetTimerProcessContact() << "," << sys->GetTimerLSsetup() << "," << sys->GetTimerLSsolve()
               << std::endl;

            if (data_out_frame % 10 == 0) {
                cout << std::fixed << std::setprecision(3) << "t = " << time << "  phase " << phase << "  disp = " << std::setprecision(4) << 1e3 * shear_disp()
                     << " mm  ratio = " << tau / p.sigma << "  e_lid = " << e_lid << "  beads = " << n_in << "  contacts = " << st.num_contacts
                     << "  step = " << std::setprecision(2) << 1e3 * step_time << " ms  overflows = " << overflows_total << endl;
                cout.unsetf(std::ios::fixed);
            }
            data_out_frame++;
        }
    }

    double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall_start).count();
    write_meta(true, wall, overflows_total, "finished");
    cout << std::setprecision(8) << "Finished at t = " << sys->GetChTime() << " s, wall clock " << wall << " s, history overflows " << overflows_total << endl;

    return 0;
}

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
// Soft-sphere (SMC) drained triaxial test on a cylindrical specimen of steel spheres, with Chrono::Multicore: a
// replication of the Cui, O'Sullivan and O'Neill (2007) tests (5 mm spheres, or a mixture of 4, 5 and 6 mm spheres,
// under a cell pressure of 80 kPa).
//
// The specimen is a cylinder 101 mm across between two rigid platens, wider than the specimen so that no sphere can
// leave it past a platen edge as it bulges (in the laboratory, the membrane is sealed to the platens). Its lateral
// boundary is a stress-controlled membrane, a set of forces towards the axis on the spheres at the surface (force
// accumulators), not a body. With --membrane voronoi (Cui et al.'s rule), the membrane spheres are those near the
// surface that touch no sphere farther out, and each carries the cell pressure on its Voronoi cell of the unrolled
// membrane. With --membrane bins, the spheres are binned by angle and height, and the one that reaches farthest out in
// each bin carries the pressure on the bin's area. The forces are updated every --membrane_every steps.
//
// The run has three phases:
// - fill        : the spheres, generated at random positions in a tall column inside a rigid mold (a ring of box
//                 plates standing on the bottom platen), settle under gravity
// - consolidate : at t_place the top platen is placed on the sample, driven by a stress servo, and the membrane is
//                 switched on; the membrane pressure and the servo's target rise from zero to the consolidation
//                 pressure by t_mold, when the mold, unloaded by the membrane, is lifted clear of the specimen
// - test        : from t_test, either
//                 triaxial:  the top platen moves down at constant speed (the axial strain rate times the height at
//                            the start of the test), the bottom platen is fixed, the membrane holds sigma3
//                 isotropic: the platen servo and the membrane follow a list of pressures, from each to the next at a
//                            constant rate, with a hold at each (compression with unload-reload loops)
//
// The friction of a contact is the larger of its two materials' values, so that the platens' friction (larger than the
// spheres') applies at the sphere-platen contacts; the mold has no friction of its own. The sphere friction is mu_fill
// while the spheres settle in the mold (the deposition), then mu from the placement of the top platen on.
//
// Stresses and strains are positive in compression. The axial stress is the contact force of the spheres on a platen
// over the current cross-section: the volume enclosed by the membrane over the distance between the platens. The
// global reference frame has Y up, along the specimen axis; the bottom platen's face is at y = 0. All units SI.
// Run-time visualization with Chrono::VSG is enabled with --vsg.
//
// Outputs, in the output directory:
//   meta.json           all parameters, run information and summary
//   timeseries.csv      one row per output interval (written as the run goes)
//   particles/*.csv     particle snapshots (triaxial: at the axial strains in --snaps, and frames; isotropic: at the end
//                       of each hold)
//   contacts/*.csv      contact snapshots, at the same times
//   states/*.bin        the state at each snapshot or frame (with --snap_states or --frame_states)
// A run can start from a settled sample (--save_fill, --load_fill: one filling for several tests) or resume the test
// from a saved state (--save_state at the start of the test, or a snapshot state; --load_state).
//
// With --render, the program renders saved states (one file, or all state files in a directory) to images with
// Chrono::VSG and exits, without simulating. The specimen options must match those of the run that saved the states.
// The camera is fixed on the -z side, framing the specimen of the first state. Each state gives one image per layer:
//   grains/   the spheres behind the vertical mid-plane z = 0 (the cut face), light gray or in gray stripes by height
//   network/  the contact force network in a slab behind the mid-plane; cylinder radius and color (light to dark
//             orange) proportional to the normal force, up to 3 times a reference force
//   ratio/    the same network colored by the Coulomb ratio |f_t| / (mu f_n), from 0 to 1 (with --render_ratio)
//   walls/    the platens alone (to separate the network from the platens when compositing)
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
#include <unordered_map>
#include <vector>

#include "chrono/ChConfig.h"
#include "chrono/ChVersion.h"
#include "chrono/assets/ChVisualShapeSphere.h"
#include "chrono/core/ChDataPath.h"
#include "chrono/physics/ChContactMaterial.h"
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
    // Spheres
    std::vector<double> radii;      // sphere radii, m
    std::vector<double> fractions;  // number fractions of the radii (normalized)
    int num_spheres;                // spheres in the specimen
    double density;                 // kg/m^3
    double E;                       // Young's modulus (spheres, platens and mold), Pa
    double nu;                      // Poisson ratio
    double cor;                     // coefficient of restitution
    double mu;                      // sphere friction
    double mu_platen;               // platen friction
    double mu_fill;                 // sphere friction while the spheres settle in the mold

    // Specimen
    double diameter;         // specimen diameter (inner diameter of the mold), m
    double h_platen;         // platen thickness, m
    double platen_diameter;  // platen diameter, wider than the specimen, m
    double wall;             // mold wall thickness, m
    int plates;              // plates in the mold ring
    double gen_fraction;     // solid fraction of the column in which the spheres are generated
    double bin_size;         // size of the membrane bins, in mean sphere diameters
    std::string membrane;    // membrane model: voronoi or bins
    int membrane_every;      // steps between updates of the membrane forces

    // Loading
    std::string test;          // triaxial or isotropic
    double sigma3;             // cell pressure of the triaxial test, Pa
    double strain_rate;        // axial strain rate, 1/s
    double strain_end;         // axial strain at the end of the triaxial test
    std::vector<double> path;  // pressures of the isotropic test, Pa (the first is the consolidation pressure)
    double ramp;               // pressure rate of the isotropic test, Pa/s
    double hold;               // hold at each pressure of the isotropic test, s
    double gravity;            // m/s^2
    double t_place;            // place the top platen and switch the membrane on, s
    double t_mold;             // lift the mold clear of the specimen, s
    double t_test;             // start of the test, s
    double servo_gain;         // platen speed of the stress servo per unit stress error, m/s/Pa
    double servo_vmax;         // largest platen speed of the stress servo, m/s
    std::string history;       // tangential displacement model: multistep, onestep, none
    std::string history_test;  // model during the test ("same": as --history)

    // Numerics
    double dt;
    int threads;
    int seed;

    // Output
    std::string out_dir;
    double out_dt;              // time series interval, s
    std::vector<double> snaps;  // axial strains of the snapshots of the triaxial test
    double frame_strain;        // axial strain between frames of the triaxial test (0: none)
    std::string save_fill;      // file for the settled sample
    std::string load_fill;      // start from a settled sample
    bool fill_only;             // stop after saving the settled sample
    std::string save_state;     // file for the state at the start of the test
    std::string load_state;     // resume the test from a saved state
    bool snap_states;           // also save the state at each snapshot
    bool frame_states;          // also save the state at each frame
    bool vsg;
    std::string provenance;  // JSON text added to meta.json (git state, run command)

    // Rendering of saved states
    std::string render;      // state file, or directory of state files, to render
    std::string render_out;  // image directory
    int render_width;        // image size, pixels
    int render_height;
    double render_slab;     // thickness of the force-network slab behind the mid-plane, m
    double render_fref;     // reference normal force of the network, N (0: mean of the first state)
    double render_fmin;     // contacts below this fraction of the reference force are not drawn
    double render_stripes;  // height of the gray stripes of the spheres, m (0: no stripes)
    bool render_ratio;      // also render the network colored by the Coulomb ratio
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
    ChCLI cli(argv[0], "Drained triaxial test with a stress-controlled membrane (Chrono::Multicore, SMC)");

    cli.AddOption<std::vector<double>>("Spheres", "radii", "Sphere radii (mm)", "2.5");
    cli.AddOption<std::vector<double>>("Spheres", "fractions", "Number fractions of the radii (default: equal)", "1");
    cli.AddOption<int>("Spheres", "num_spheres", "Spheres in the specimen", "15400");
    cli.AddOption<double>("Spheres", "density", "Sphere density (kg/m3)", "7800");
    cli.AddOption<double>("Spheres", "E", "Young's modulus of spheres, platens and mold (Pa)", "2e8");
    cli.AddOption<double>("Spheres", "nu", "Poisson ratio", "0.28");
    cli.AddOption<double>("Spheres", "cor", "Coefficient of restitution", "0.597");
    cli.AddOption<double>("Spheres", "mu", "Sphere friction coefficient", "0.096");
    cli.AddOption<double>("Spheres", "mu_platen", "Platen friction coefficient", "0.228");
    cli.AddOption<double>("Spheres", "mu_fill", "Sphere friction while the spheres settle in the mold (default: mu)");

    cli.AddOption<double>("Specimen", "diameter", "Specimen diameter: inner diameter of the mold (m)", "0.101");
    cli.AddOption<double>("Specimen", "h_platen", "Platen thickness (m)", "0.01");
    cli.AddOption<double>("Specimen", "platen_diameter", "Platen diameter (m, 0: 1.5 times the specimen diameter)", "0");
    cli.AddOption<double>("Specimen", "wall", "Mold wall thickness (m)", "0.005");
    cli.AddOption<int>("Specimen", "plates", "Plates in the mold ring", "120");
    cli.AddOption<double>("Specimen", "gen_fraction", "Solid fraction of the column in which the spheres are generated", "0.3");
    cli.AddOption<double>("Specimen", "bin_size", "Size of the membrane bins; for voronoi, height of the layers that find the membrane spheres (mean sphere diameters)", "1");
    cli.AddOption<std::string>("Specimen", "membrane", "Membrane model: voronoi (Cui et al. 2007) or bins", "voronoi");
    cli.AddOption<int>("Specimen", "membrane_every", "Steps between updates of the membrane forces", "10");

    cli.AddOption<std::string>("Loading", "test", "Test: triaxial or isotropic", "triaxial");
    cli.AddOption<double>("Loading", "sigma3", "Cell pressure of the triaxial test (Pa)", "80000");
    cli.AddOption<double>("Loading", "strain_rate", "Axial strain rate of the triaxial test (1/s)", "0.05");
    cli.AddOption<double>("Loading", "strain_end", "Axial strain at the end of the triaxial test", "0.25");
    cli.AddOption<std::vector<double>>("Loading", "path", "Pressures of the isotropic test, the first being the consolidation pressure (Pa)", "1e4,1e5,2.5e4,1e5,4e5,1e5,4e5,1e4");
    cli.AddOption<double>("Loading", "ramp", "Pressure rate of the isotropic test (Pa/s)", "2e5");
    cli.AddOption<double>("Loading", "hold", "Hold at each pressure of the isotropic test (s)", "0.1");
    cli.AddOption<double>("Loading", "gravity", "Gravitational acceleration (m/s2)", "9.81");
    cli.AddOption<double>("Loading", "t_place", "Place the top platen and switch the membrane on (s)", "1.0");
    cli.AddOption<double>("Loading", "t_mold", "Lift the mold clear of the specimen (s)", "1.2");
    cli.AddOption<double>("Loading", "t_test", "Start of the test (s)", "2.0");
    cli.AddOption<double>("Loading", "servo_gain", "Platen speed of the stress servo per unit stress error (m/s/Pa)", "1e-6");
    cli.AddOption<double>("Loading", "servo_vmax", "Largest platen speed of the stress servo (m/s)", "0.05");
    cli.AddOption<std::string>("Loading", "history", "Tangential displacement model: multistep, onestep or none", "multistep");
    cli.AddOption<std::string>("Loading", "history_test", "Tangential displacement model during the test (same: as --history)", "same");

    cli.AddOption<double>("Numerics", "dt", "Time step (s)", "1e-5");
    cli.AddOption<int>("Numerics", "threads", "Number of OpenMP threads", "20");
    cli.AddOption<int>("Numerics", "seed", "Seed of the random sphere positions and sizes", "1");

    cli.AddOption<std::string>("Output", "out", "Output directory", "triaxial_membrane");
    cli.AddOption<double>("Output", "out_dt", "Time series interval (s)", "0.01");
    cli.AddOption<std::vector<double>>("Output", "snaps", "Axial strains of the snapshots of the triaxial test (%)", "0,1,2,3,5,7.5,10,15,20,25");
    cli.AddOption<double>("Output", "frame_strain", "Axial strain between particle and contact frames of the triaxial test (%, 0: none)", "0");
    cli.AddOption<std::string>("Output", "save_fill", "Save the settled sample (before the top platen is placed) to this file", "");
    cli.AddOption<std::string>("Output", "load_fill", "Start from a settled sample saved with --save_fill", "");
    cli.AddOption<bool>("Output", "fill_only", "Stop after the settled sample is saved");
    cli.AddOption<std::string>("Output", "save_state", "Save the sample at the start of the test to this file", "");
    cli.AddOption<std::string>("Output", "load_state", "Resume from a state saved with --save_state or --snap_states", "");
    cli.AddOption<bool>("Output", "snap_states", "Also save the state at each snapshot (states/*.bin)");
    cli.AddOption<bool>("Output", "frame_states", "Also save the state at each frame (states/*.bin)");
    cli.AddOption<bool>("Output", "vsg", "Run-time visualization with Chrono::VSG");
    cli.AddOption<std::string>("Output", "provenance", "JSON text recorded in meta.json", "");

    cli.AddOption<std::string>("Render", "render", "Render a saved state, or all states in a directory, to images and exit", "");
    cli.AddOption<std::string>("Render", "render_out", "Image directory (default: render, in the output directory)", "");
    cli.AddOption<int>("Render", "render_width", "Image width (pixels)", "1080");
    cli.AddOption<int>("Render", "render_height", "Image height (pixels)", "1080");
    cli.AddOption<double>("Render", "render_slab", "Thickness of the force-network slab behind the mid-plane (mean diameters)", "1.5");
    cli.AddOption<double>("Render", "render_fref", "Reference normal force of the network (N, 0: mean of the first state)", "0");
    cli.AddOption<double>("Render", "render_fmin", "Contacts below this fraction of the reference force are not drawn", "0.5");
    cli.AddOption<double>("Render", "render_stripes", "Height of the gray stripes of the spheres, at the first state (mm, 0: none)", "0");
    cli.AddOption<bool>("Render", "render_ratio", "Also render the network colored by the Coulomb ratio");

    if (!cli.Parse(argc, argv, true))
        return false;

    p.radii = cli.GetAsType<std::vector<double>>("radii");
    for (auto& r : p.radii)
        r *= 1e-3;
    p.fractions = cli.GetAsType<std::vector<double>>("fractions");
    if (p.fractions.size() != p.radii.size())
        p.fractions.assign(p.radii.size(), 1.0);
    double sum = 0;
    for (double f : p.fractions)
        sum += f;
    for (auto& f : p.fractions)
        f /= sum;
    p.num_spheres = cli.GetAsType<int>("num_spheres");
    p.density = cli.GetAsType<double>("density");
    p.E = cli.GetAsType<double>("E");
    p.nu = cli.GetAsType<double>("nu");
    p.cor = cli.GetAsType<double>("cor");
    p.mu = cli.GetAsType<double>("mu");
    p.mu_platen = cli.GetAsType<double>("mu_platen");
    p.mu_fill = cli.CheckOption("mu_fill") ? cli.GetAsType<double>("mu_fill") : p.mu;

    p.diameter = cli.GetAsType<double>("diameter");
    p.h_platen = cli.GetAsType<double>("h_platen");
    p.platen_diameter = cli.GetAsType<double>("platen_diameter");
    if (p.platen_diameter <= 0)
        p.platen_diameter = 1.5 * p.diameter;
    p.wall = cli.GetAsType<double>("wall");
    p.plates = cli.GetAsType<int>("plates");
    p.gen_fraction = cli.GetAsType<double>("gen_fraction");
    p.bin_size = cli.GetAsType<double>("bin_size");
    p.membrane = cli.GetAsType<std::string>("membrane");
    p.membrane_every = cli.GetAsType<int>("membrane_every");

    p.test = cli.GetAsType<std::string>("test");
    p.sigma3 = cli.GetAsType<double>("sigma3");
    p.strain_rate = cli.GetAsType<double>("strain_rate");
    p.strain_end = cli.GetAsType<double>("strain_end");
    p.path = cli.GetAsType<std::vector<double>>("path");
    p.ramp = cli.GetAsType<double>("ramp");
    p.hold = cli.GetAsType<double>("hold");
    p.gravity = cli.GetAsType<double>("gravity");
    p.t_place = cli.GetAsType<double>("t_place");
    p.t_mold = cli.GetAsType<double>("t_mold");
    p.t_test = cli.GetAsType<double>("t_test");
    p.servo_gain = cli.GetAsType<double>("servo_gain");
    p.servo_vmax = cli.GetAsType<double>("servo_vmax");
    p.history = cli.GetAsType<std::string>("history");
    p.history_test = cli.GetAsType<std::string>("history_test");
    if (p.history_test == "same")
        p.history_test = p.history;
    HistoryModel(p.history);
    HistoryModel(p.history_test);

    p.dt = cli.GetAsType<double>("dt");
    p.threads = cli.GetAsType<int>("threads");
    p.seed = cli.GetAsType<int>("seed");

    p.out_dir = cli.GetAsType<std::string>("out");
    p.out_dt = cli.GetAsType<double>("out_dt");
    p.snaps = cli.GetAsType<std::vector<double>>("snaps");
    for (auto& s : p.snaps)
        s *= 1e-2;
    p.frame_strain = 1e-2 * cli.GetAsType<double>("frame_strain");
    p.save_fill = cli.GetAsType<std::string>("save_fill");
    p.load_fill = cli.GetAsType<std::string>("load_fill");
    p.fill_only = cli.GetAsType<bool>("fill_only");
    p.save_state = cli.GetAsType<std::string>("save_state");
    p.load_state = cli.GetAsType<std::string>("load_state");
    p.snap_states = cli.GetAsType<bool>("snap_states");
    p.frame_states = cli.GetAsType<bool>("frame_states");
    p.vsg = cli.GetAsType<bool>("vsg");
    p.provenance = cli.GetAsType<std::string>("provenance");

    p.render = cli.GetAsType<std::string>("render");
    p.render_out = cli.GetAsType<std::string>("render_out");
    if (p.render_out.empty())
        p.render_out = p.out_dir + "/render";
    p.render_width = cli.GetAsType<int>("render_width");
    p.render_height = cli.GetAsType<int>("render_height");
    p.render_slab = cli.GetAsType<double>("render_slab");  // in mean diameters; converted once the sizes are known
    p.render_fref = cli.GetAsType<double>("render_fref");
    p.render_fmin = cli.GetAsType<double>("render_fmin");
    p.render_stripes = 1e-3 * cli.GetAsType<double>("render_stripes");
    p.render_ratio = cli.GetAsType<bool>("render_ratio");

    if (p.radii.empty() || p.num_spheres <= 0) {
        std::cerr << "No spheres" << endl;
        return false;
    }
    if (p.test != "triaxial" && p.test != "isotropic") {
        std::cerr << "Unknown test: " << p.test << endl;
        return false;
    }
    if (p.membrane != "voronoi" && p.membrane != "bins") {
        std::cerr << "Unknown membrane model: " << p.membrane << endl;
        return false;
    }
    if (p.membrane_every < 1) {
        std::cerr << "--membrane_every must be at least 1" << endl;
        return false;
    }
    if (p.test == "isotropic" && p.path.empty()) {
        std::cerr << "The isotropic test needs at least one pressure in --path" << endl;
        return false;
    }
    if (!(p.t_place < p.t_mold && p.t_mold <= p.t_test)) {
        std::cerr << "Need t_place < t_mold <= t_test" << endl;
        return false;
    }
    if (p.fill_only && p.save_fill.empty()) {
        std::cerr << "--fill_only needs --save_fill" << endl;
        return false;
    }

    return true;
}

// -----------------------------------------------------------------------------
// Friction of a contact: the larger of its two materials' values
// -----------------------------------------------------------------------------

class FrictionMax : public ChContactMaterialCompositionStrategy {
  public:
    virtual float CombineFriction(float a1, float a2) const override { return std::max(a1, a2); }
};

// -----------------------------------------------------------------------------
// Pressure path of the isotropic test: from each pressure to the next at a constant rate, then a hold
// -----------------------------------------------------------------------------

struct PressurePath {
    std::vector<double> pressures;
    double rate;
    double hold;

    /// Time from the start of the test to the end of the hold after pressure k + 1.
    double SegmentEnd(size_t k) const {
        double t = 0;
        for (size_t j = 0; j <= k && j + 1 < pressures.size(); j++)
            t += std::abs(pressures[j + 1] - pressures[j]) / rate + hold;
        return t;
    }

    double Duration() const { return pressures.size() > 1 ? SegmentEnd(pressures.size() - 2) : 0.0; }

    /// Pressure at time t after the start of the test.
    double Pressure(double t) const {
        for (size_t k = 0; k + 1 < pressures.size(); k++) {
            double ramp_time = std::abs(pressures[k + 1] - pressures[k]) / rate;
            if (t < ramp_time)
                return pressures[k] + (pressures[k + 1] - pressures[k]) * t / ramp_time;
            t -= ramp_time;
            if (t < hold)
                return pressures[k + 1];
            t -= hold;
        }
        return pressures.back();
    }
};

// -----------------------------------------------------------------------------
// Membrane: forces on the spheres at the specimen's surface, standing in for the cell pressure on a latex membrane,
// through the spheres' first force accumulator. Two models:
// - VORONOI (Cui et al. 2007): the membrane spheres are those within one mean diameter of the outer surface of their
//   layer that touch no sphere farther out (no contact normal points radially outwards, that is, within 60 degrees of
//   the outward radial direction; neighbors on the surface touch at about 90 degrees). Each carries the pressure on
//   its Voronoi cell of the unrolled membrane (arc length around the y axis, periodic, and height between the platen
//   faces), found on a raster of cells an eighth of a mean diameter wide, at its own outer radius.
// - BINS: the spheres are binned by angle and by height; in each bin the sphere that reaches farthest out carries the
//   pressure on the bin's area at its outer radius. The bins of a layer that hold no sphere share their angle out among
//   the layer's other bins, so that each layer carries its full load.
// The forces point to the axis, through the sphere centers, and stay in the accumulators until the next update. Spheres
// beyond the platen faces carry none.
// -----------------------------------------------------------------------------

class Membrane {
  public:
    enum class Model { VORONOI, BINS };

    Membrane(ChSystemMulticoreSMC* sys,
             const std::vector<std::shared_ptr<ChBody>>& spheres,
             const std::vector<double>& radii,
             int first_particle,
             double radius,
             double d_mean,
             double bin_size,
             Model model)
        : m_sys(sys), m_spheres(spheres), m_radii(radii), m_first(first_particle), m_radius(radius), m_d_mean(d_mean), m_bin_size(bin_size), m_model(model) {}

    /// Apply the pressure on the spheres between the platen faces at y0 and y1 (pressure 0: clear the forces), and
    /// compute the enclosed volume and the membrane statistics.
    void Apply(double pressure, double y0, double y1) {
        size_t n = m_spheres.size();
        m_rho.resize(n);
        m_reach.resize(n);
        m_inside.resize(n);
        num_outside = 0;
        for (size_t i = 0; i < n; i++) {
            m_spheres[i]->EmptyAccumulator(0);
            const auto& q = m_spheres[i]->GetPos();
            m_inside[i] = q.y() >= y0 && q.y() <= y1;
            if (!m_inside[i])
                num_outside++;
            m_rho[i] = std::sqrt(q.x() * q.x() + q.z() * q.z());
            m_reach[i] = m_rho[i] + m_radii[i];
        }
        volume = 0;
        force = 0;
        num_carriers = 0;
        num_bins = 0;
        if (m_model == Model::VORONOI)
            ApplyVoronoi(pressure, y0, y1);
        else
            ApplyBins(pressure, y0, y1);
    }

    double volume = 0;     // volume enclosed by the membrane, m^3
    double force = 0;      // sum of the membrane forces, N
    int num_carriers = 0;  // spheres carrying the membrane load
    int num_outside = 0;   // spheres above the top platen face or below the bottom one
    int num_bins = 0;      // bins (BINS) or raster cells (VORONOI)

  private:
    // Force f on sphere i, towards the axis
    void Push(size_t i, double f) {
        if (m_rho[i] < 1e-9 || f <= 0)
            return;
        const auto& q = m_spheres[i]->GetPos();
        m_spheres[i]->AccumulateForce(0, ChVector3d(-q.x() / m_rho[i], 0, -q.z() / m_rho[i]) * f, q, false);
        force += f;
        num_carriers++;
    }

    void ApplyBins(double pressure, double y0, double y1) {
        double h = y1 - y0;
        int nt = std::max(8, (int)std::lround(2 * CH_PI * m_radius / m_bin_size));
        int nz = std::max(1, (int)std::lround(h / m_bin_size));
        double dth = 2 * CH_PI / nt;
        double dz = h / nz;

        std::vector<int> carrier(nt * nz, -1);
        for (size_t i = 0; i < m_spheres.size(); i++) {
            if (!m_inside[i])
                continue;
            const auto& q = m_spheres[i]->GetPos();
            int it = std::clamp((int)std::floor((std::atan2(q.z(), q.x()) + CH_PI) / dth), 0, nt - 1);
            int iz = std::clamp((int)std::floor((q.y() - y0) / dz), 0, nz - 1);
            int b = iz * nt + it;
            if (carrier[b] < 0 || m_reach[i] > m_reach[carrier[b]])
                carrier[b] = (int)i;
        }

        num_bins = nt * nz;
        for (int iz = 0; iz < nz; iz++) {
            int count = 0;
            double sum = 0;
            for (int it = 0; it < nt; it++) {
                int c = carrier[iz * nt + it];
                if (c >= 0) {
                    count++;
                    sum += m_reach[c];
                }
            }
            if (count == 0)
                continue;
            double mean = sum / count;
            double dth_carrier = 2 * CH_PI / count;
            for (int it = 0; it < nt; it++) {
                int c = carrier[iz * nt + it];
                double rb = c >= 0 ? m_reach[c] : mean;
                volume += 0.5 * rb * rb * dth * dz;
                if (c >= 0 && pressure > 0)
                    Push(c, pressure * m_reach[c] * dth_carrier * dz);
            }
        }
    }

    void ApplyVoronoi(double pressure, double y0, double y1) {
        double h = y1 - y0;
        size_t n = m_spheres.size();

        // Candidates: the spheres within one mean diameter of the outer surface of their layer (layers one bin high)
        int nl = std::max(1, (int)std::lround(h / m_bin_size));
        double dl = h / nl;
        std::vector<int> layer(n, 0);
        std::vector<double> layer_reach(nl, 0.0);
        for (size_t i = 0; i < n; i++) {
            if (!m_inside[i])
                continue;
            layer[i] = std::clamp((int)std::floor((m_spheres[i]->GetPos().y() - y0) / dl), 0, nl - 1);
            layer_reach[layer[i]] = std::max(layer_reach[layer[i]], m_reach[i]);
        }
        std::vector<char> member(n, 0);
        for (size_t i = 0; i < n; i++)
            member[i] = m_inside[i] && m_reach[i] >= layer_reach[layer[i]] - m_d_mean;

        // A candidate that touches a sphere farther out, the contact normal within 60 degrees of its own outward radial
        // direction, is not a membrane sphere (contacts of the last step)
        auto cd = m_sys->data_manager->cd_data;
        if (cd) {
            for (uint k = 0; k < cd->num_rigid_contacts; k++) {
                if (cd->dpth_rigid_rigid[k] >= 0)
                    continue;
                int a = cd->bids_rigid_rigid[k].x - m_first;
                int b = cd->bids_rigid_rigid[k].y - m_first;
                if (a < 0 || b < 0 || a >= (int)n || b >= (int)n)
                    continue;
                for (auto [i, j] : {std::pair<int, int>(a, b), std::pair<int, int>(b, a)}) {
                    if (!member[i] || m_rho[i] < 1e-9)
                        continue;
                    const auto& qi = m_spheres[i]->GetPos();
                    const auto& qj = m_spheres[j]->GetPos();
                    ChVector3d nij = qj - qi;
                    if ((nij.x() * qi.x() + nij.z() * qi.z()) / m_rho[i] > 0.5 * nij.Length())
                        member[i] = 0;
                }
            }
        }

        // Membrane spheres on the unrolled membrane: arc length at their mean outer radius, and height
        std::vector<int> ids;
        double r_mean = 0;
        for (size_t i = 0; i < n; i++) {
            if (member[i]) {
                ids.push_back((int)i);
                r_mean += m_reach[i];
            }
        }
        if (ids.empty())
            return;
        r_mean /= ids.size();
        double L = 2 * CH_PI * r_mean;
        std::vector<double> s(ids.size());
        std::vector<double> z(ids.size());
        for (size_t m = 0; m < ids.size(); m++) {
            const auto& q = m_spheres[ids[m]]->GetPos();
            s[m] = r_mean * (std::atan2(q.z(), q.x()) + CH_PI);
            z[m] = q.y() - y0;
        }

        // Buckets of membrane spheres, about one mean diameter wide, periodic around the axis
        int bs = std::max(1, (int)std::floor(L / m_d_mean));
        int bz = std::max(1, (int)std::floor(h / m_d_mean));
        double ws = L / bs;
        double wz = h / bz;
        double w_min = std::min(ws, wz);
        std::vector<std::vector<int>> buckets(bs * bz);
        for (size_t m = 0; m < ids.size(); m++) {
            int ib = std::clamp((int)std::floor(s[m] / ws), 0, bs - 1);
            int jb = std::clamp((int)std::floor(z[m] / wz), 0, bz - 1);
            buckets[jb * bs + ib].push_back((int)m);
        }

        // Raster: each cell goes to the nearest membrane sphere, so the cells of a sphere make up its Voronoi cell.
        // The search goes out ring by ring of buckets; spheres in rings not yet searched are at least ring * w_min away.
        int ns = std::max(16, (int)std::lround(8 * L / m_d_mean));
        int nz = std::max(1, (int)std::lround(8 * h / m_d_mean));
        double hs = L / ns;
        double hz = h / nz;
        num_bins = ns * nz;
        std::vector<int> cells(ids.size(), 0);
        for (int jz = 0; jz < nz; jz++) {
            double zc = (jz + 0.5) * hz;
            int jb0 = std::clamp((int)std::floor(zc / wz), 0, bz - 1);
            for (int is = 0; is < ns; is++) {
                double sc = (is + 0.5) * hs;
                int ib0 = std::clamp((int)std::floor(sc / ws), 0, bs - 1);
                int best = -1;
                double best_d2 = 0;
                for (int ring = 0; ring <= bs + bz; ring++) {
                    for (int dj = -ring; dj <= ring; dj++) {
                        int jb = jb0 + dj;
                        if (jb < 0 || jb >= bz)
                            continue;
                        for (int di = -ring; di <= ring; di++) {
                            if (std::max(std::abs(di), std::abs(dj)) != ring)
                                continue;
                            int ib = ((ib0 + di) % bs + bs) % bs;
                            for (int m : buckets[jb * bs + ib]) {
                                double ds = s[m] - sc;
                                ds -= L * std::round(ds / L);
                                double d2 = ds * ds + (z[m] - zc) * (z[m] - zc);
                                if (best < 0 || d2 < best_d2) {
                                    best = m;
                                    best_d2 = d2;
                                }
                            }
                        }
                    }
                    if (best >= 0 && best_d2 <= (ring * w_min) * (ring * w_min))
                        break;
                }
                if (best < 0)
                    continue;
                cells[best]++;
                double rb = m_reach[ids[best]];
                volume += 0.5 * rb * rb * (hs / r_mean) * hz;
            }
        }

        // The pressure on each Voronoi cell, scaled from the mean outer radius to the sphere's own
        if (pressure > 0) {
            for (size_t m = 0; m < ids.size(); m++)
                Push(ids[m], pressure * (m_reach[ids[m]] / r_mean) * cells[m] * hs * hz);
        }
    }

    ChSystemMulticoreSMC* m_sys;
    const std::vector<std::shared_ptr<ChBody>>& m_spheres;
    const std::vector<double>& m_radii;
    int m_first;        // index of the first sphere body
    double m_radius;    // nominal specimen radius (sets the number of bins around the axis)
    double m_d_mean;    // mean sphere diameter, m
    double m_bin_size;  // bin size, m (BINS; for VORONOI, the height of the candidate layers)
    Model m_model;
    std::vector<double> m_rho;    // distance of each sphere center from the axis
    std::vector<double> m_reach;  // distance of each sphere's outer side from the axis
    std::vector<char> m_inside;   // whether each sphere's center is between the platen faces
};

// Values saved with a state, after the bodies: phase, start of the test, height, volume and mean radius at the start
// of the test, solid volume (to check that a loaded state has the same spheres)
static const size_t kStateExtra = 6;

#ifdef CHRONO_VSG

// -----------------------------------------------------------------------------
// Rendering of saved states (--render)
// -----------------------------------------------------------------------------

// Sphere colors: light gray, or alternating light and medium gray stripes set from the sphere heights at the first
// rendered state (the stripes then move with the spheres, which shows how the specimen deforms)
class SphereColors : public vsg3d::ChInstancedShapesVSG::ColorCallback {
  public:
    virtual ChColor GetColor(const ChBody& body) const override {
        auto s = stripe.find(&body);
        bool dark = s != stripe.end() && s->second % 2 != 0;
        return dark ? ChColor(0.55f, 0.55f, 0.55f) : ChColor(0.9f, 0.9f, 0.9f);
    }

    std::unordered_map<const ChBody*, int> stripe;
};

// Spheres behind the vertical mid-plane z = 0, seen from the -z side
class SphereCut : public vsg3d::ChInstancedShapesVSG::VisibilityCallback {
  public:
    virtual bool IsVisible(const ChBody& body) const override { return body.GetPos().z() > 0; }
};

int RenderStates(const Params& p,
                 ChSystemMulticoreSMC* sys,
                 std::shared_ptr<ChBody> bottom,
                 std::shared_ptr<ChBody> mold,
                 std::shared_ptr<ChBody> top,
                 const std::vector<std::shared_ptr<ChBody>>& spheres,
                 double radius,
                 double d_mean) {
    namespace fs = std::filesystem;
    using NetworkColoring = ChMulticoreVisualizationVSG::NetworkColoring;

    // States to render: one file, or all state files in a directory, in name order
    std::vector<fs::path> files;
    if (fs::is_directory(p.render)) {
        for (const auto& entry : fs::directory_iterator(p.render))
            if (entry.path().extension() == ".bin")
                files.push_back(entry.path());
        std::sort(files.begin(), files.end());
    } else if (fs::exists(p.render)) {
        files.push_back(p.render);
    }
    if (files.empty()) {
        std::cerr << "No state files in " << p.render << endl;
        return 1;
    }

    std::vector<std::string> layers = {"grains", "network", "walls"};
    if (p.render_ratio)
        layers.push_back("ratio");
    for (const auto& layer : layers) {
        std::error_code ec;
        fs::create_directories(p.render_out + "/" + layer, ec);
        if (!fs::is_directory(p.render_out + "/" + layer)) {
            std::cerr << "Error creating directory " << p.render_out + "/" + layer << endl;
            return 1;
        }
    }

    std::vector<double> extra;
    auto load = [&](const fs::path& file) {
        if (!LoadState(*sys, file.string(), extra) || extra.size() != kStateExtra) {
            std::cerr << "Cannot load state from " << file.string() << " (missing file or different specimen)" << endl;
            return false;
        }
        return true;
    };

    // The first state places the camera
    if (!load(files[0]))
        return 1;
    sys->GetSettings()->solver.tangential_displ_mode = HistoryModel((int)extra[0] == 2 ? p.history_test : p.history);

    // Fixed camera on the -z side, looking at the mid-plane, framing the specimen with room for its bulging (a narrow
    // field of view keeps the perspective small)
    double slab = p.render_slab * d_mean;
    double y_face = top->GetPos().y() - p.h_platen / 2;
    double x0 = -1.3 * radius;
    double x1 = 1.3 * radius;
    double y0 = -p.h_platen;
    double y1 = y_face + p.h_platen;
    double aspect = (double)p.render_width / p.render_height;
    double fov = 6;
    double half_height = 0.54 * std::max(y1 - y0, (x1 - x0) / aspect);
    double dist = half_height / std::tan(0.5 * fov * CH_DEG_TO_RAD);
    ChVector3d target((x0 + x1) / 2, (y0 + y1) / 2, 0);

    // Platens and mold in a light gray; the mold only while it surrounds the specimen (once lifted, it would show at the
    // top of the image)
    bool show_mold = mold->GetPos().y() < 1e-6;
    for (const auto& body : {bottom, mold, top}) {
        for (const auto& instance : body->GetVisualModel()->GetShapeInstances()) {
            instance.shape->SetColor(ChColor(0.85f, 0.85f, 0.85f));
            if (body == mold)
                instance.shape->SetVisible(show_mold);
        }
    }

    auto colors = chrono_types::make_shared<SphereColors>();

    auto mcore_vis = chrono_types::make_shared<ChMulticoreVisualizationVSG>(sys);
    mcore_vis->SetBodySelection(spheres);
    mcore_vis->SetBodyColorCallback(colors);
    mcore_vis->SetBodyVisibilityCallback(chrono_types::make_shared<SphereCut>());
    mcore_vis->SetNetworkMode(vsg3d::ChForceNetworkVSG::Mode::CYLINDERS);
    mcore_vis->SetNetworkColorRamp(ChColor(0.957f, 0.694f, 0.514f), ChColor(0.518f, 0.235f, 0.047f));  // #F4B183 to #843C0C
    mcore_vis->SetRegionSlab(ChVector3d(0, 0, slab / 2), ChVector3d(0, 0, 1), slab / 2);

    auto vis = chrono_types::make_shared<vsg3d::ChVisualSystemVSG>();
    vis->AttachPlugin(mcore_vis);
    vis->SetWindowTitle("Triaxial test");
    vis->SetWindowSize(p.render_width, p.render_height);
    vis->SetCameraVertical(CameraVerticalDir::Y);
    vis->AddCamera(target - ChVector3d(0, 0, dist), target);
    vis->SetCameraAngleDeg(fov);
    vis->SetBackgroundColor(ChColor(1.0f, 1.0f, 1.0f));
    vis->SetLightIntensity(1.0f);
    vis->SetLightDirection(1.1 * CH_PI, CH_PI / 4);
    vis->SetGuiVisibility(false);
    vis->HideLogo();
    vis->Initialize();

    double f_ref = p.render_fref;
    for (size_t k = 0; k < files.size(); k++) {
        if (!vis->Run())
            break;
        if (!load(files[k]))
            return 1;

        // One step computes the contacts and contact forces of the saved state (the membrane forces do not enter them)
        y_face = top->GetPos().y() - p.h_platen / 2;
        sys->DoStepDynamics(p.dt);

        double fn_sum = 0;
        int fn_num = 0;
        for (const auto& c : GetContacts(*sys)) {
            if (c.fn > 0) {
                fn_sum += c.fn;
                fn_num++;
            }
        }
        double fn_mean = fn_num > 0 ? fn_sum / fn_num : 0;
        if (f_ref <= 0)
            f_ref = fn_mean;

        if (k == 0 && p.render_stripes > 0) {
            for (const auto& s : spheres)
                colors->stripe[s.get()] = (int)std::floor(s->GetPos().y() / p.render_stripes);
        }

        mcore_vis->SetNetworkRadius(0, 0.2 * d_mean, 0, 3 * f_ref);
        mcore_vis->SetNetworkColoring(NetworkColoring::NORMAL_FORCE, 0, 3 * f_ref);
        mcore_vis->SetNetworkForceThreshold(p.render_fmin * f_ref);

        // An image request captures the frame of the previous Render() and is written during the next one
        std::string stem = files[k].stem().string();
        auto shot = [&](const std::string& layer, bool show_spheres, bool network) {
            mcore_vis->SetBodyVisibility(show_spheres);
            mcore_vis->SetNetworkVisibility(network);
            vis->Render();
            vis->WriteImageToFile(p.render_out + "/" + layer + "/" + stem + ".bmp");
        };
        shot("grains", true, false);
        shot("network", false, true);
        unsigned int segments = mcore_vis->GetNumNetworkSegments();
        if (p.render_ratio) {
            mcore_vis->SetNetworkColoring(NetworkColoring::COULOMB_RATIO, 0, 1);
            shot("ratio", false, true);
        }
        shot("walls", false, false);
        vis->Render();

        double eps_a = (int)extra[0] == 2 ? (extra[2] - y_face) / extra[2] : 0.0;
        cout << stem << ": t = " << sys->GetChTime() << " s, axial strain " << 100 * eps_a << " %, mean normal force " << fn_mean << " N (reference " << f_ref << " N), "
             << segments << " half-segments drawn" << endl;
    }

    return 0;
}

#endif

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
    if (!p.render.empty())
        dirs.clear();
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

    double R = p.diameter / 2;
    double R_platen = p.platen_diameter / 2;
    bool triaxial = p.test == "triaxial";
    PressurePath path{p.path, p.ramp, p.hold};

    // Sphere sizes: the count of each radius from its fraction, in random order
    std::mt19937 rng((unsigned int)p.seed);
    std::vector<double> radius_of;
    std::vector<int> size_counts;
    {
        int assigned = 0;
        for (size_t k = 0; k < p.radii.size(); k++) {
            int n = (k + 1 < p.radii.size()) ? (int)std::lround(p.fractions[k] * p.num_spheres) : p.num_spheres - assigned;
            n = std::max(0, std::min(n, p.num_spheres - assigned));
            size_counts.push_back(n);
            assigned += n;
            radius_of.insert(radius_of.end(), n, p.radii[k]);
        }
        std::shuffle(radius_of.begin(), radius_of.end(), rng);
    }
    double r_max = *std::max_element(p.radii.begin(), p.radii.end());
    double d_mean = 0;
    double solid_volume = 0;
    for (double r : radius_of) {
        d_mean += 2 * r;
        solid_volume += (4.0 / 3) * CH_PI * r * r * r;
    }
    d_mean /= p.num_spheres;

    // Column in which the spheres are generated, and the mold around it. The column holds the spheres, each with half
    // the generation gap around it, at the solid fraction gen_fraction (below the jamming fraction of random sequential
    // addition, about 0.38).
    double gen_gap = 0.0002;
    double gen_volume = 0;
    for (double r : radius_of)
        gen_volume += (4.0 / 3) * CH_PI * std::pow(r + gen_gap / 2, 3);
    double h_gen = gen_volume / (p.gen_fraction * CH_PI * (R - r_max) * (R - r_max));
    double h_mold = h_gen + 4 * r_max;

    // Random positions in the column, apart from one another (random sequential addition, with a grid of cells one
    // largest diameter wide)
    std::vector<ChVector3d> positions;
    {
        std::uniform_real_distribution<double> uni(0.0, 1.0);
        double cell = 2 * r_max + gen_gap;
        int nxz = (int)std::ceil(2 * R / cell) + 1;
        int ny = (int)std::ceil(h_gen / cell) + 1;
        std::vector<std::vector<int>> grid(nxz * nxz * ny);
        auto cell_index = [&](int i, int j, int k) { return (j * nxz + i) * nxz + k; };
        for (int n = 0; n < p.num_spheres; n++) {
            double r = radius_of[n];
            bool placed = false;
            for (int attempt = 0; attempt < 100000 && !placed; attempt++) {
                double rr = (R - r - gen_gap) * std::sqrt(uni(rng));
                double th = 2 * CH_PI * uni(rng);
                ChVector3d q(rr * std::cos(th), r + gen_gap + (h_gen - 2 * r - gen_gap) * uni(rng), rr * std::sin(th));
                int ci = (int)std::floor((q.x() + R) / cell);
                int cj = (int)std::floor(q.y() / cell);
                int ck = (int)std::floor((q.z() + R) / cell);
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
                std::cerr << "Cannot place sphere " << n << " in the column; use a smaller --gen_fraction" << endl;
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
    sys->SetMaterialCompositionStrategy(std::make_unique<FrictionMax>());

    int threads = std::min(p.threads, omp_get_num_procs());
    sys->SetNumThreads(threads);

    sys->GetSettings()->solver.use_full_inertia_tensor = false;
    sys->GetSettings()->solver.tolerance = 0.01;
    sys->GetSettings()->solver.max_iteration_bilateral = 100;
    sys->GetSettings()->solver.clamp_bilaterals = false;
    sys->GetSettings()->solver.bilateral_clamp_speed = 0.1;

    sys->GetSettings()->solver.contact_force_model = ChSystemSMC::ContactForceModel::Hertz;
    sys->GetSettings()->solver.tangential_displ_mode = HistoryModel(p.history);

    sys->GetSettings()->collision.bins_per_axis = vec3(10, 20, 10);
    sys->GetSettings()->collision.narrowphase_algorithm = ChNarrowphase::Algorithm::HYBRID;

    // Materials: spheres, platens, and the mold (no friction of its own: a sphere-mold contact takes the sphere's)

    auto mat_sphere = chrono_types::make_shared<ChContactMaterialSMC>();
    mat_sphere->SetYoungModulus((float)p.E);
    mat_sphere->SetPoissonRatio((float)p.nu);
    mat_sphere->SetRestitution((float)p.cor);
    mat_sphere->SetFriction((float)p.mu_fill);

    auto mat_platen = chrono_types::make_shared<ChContactMaterialSMC>();
    mat_platen->SetYoungModulus((float)p.E);
    mat_platen->SetPoissonRatio((float)p.nu);
    mat_platen->SetRestitution((float)p.cor);
    mat_platen->SetFriction((float)p.mu_platen);

    auto mat_mold = chrono_types::make_shared<ChContactMaterialSMC>();
    mat_mold->SetYoungModulus((float)p.E);
    mat_mold->SetPoissonRatio((float)p.nu);
    mat_mold->SetRestitution((float)p.cor);
    mat_mold->SetFriction(0.0f);

    // The platens and the mold are created before the spheres, so that the history of their contacts is stored on the
    // spheres (the body with the larger index)

    // Bottom platen (body 0), face at y = 0
    auto bottom = chrono_types::make_shared<ChBody>();
    bottom->SetPos(ChVector3d(0, -p.h_platen / 2, 0));
    bottom->SetFixed(true);
    bottom->EnableCollision(true);
    utils::AddCylinderGeometry(bottom.get(), mat_platen, R_platen, p.h_platen, ChVector3d(0, 0, 0), z2y, true);
    bottom->GetCollisionModel()->SetFamily(1);
    bottom->GetCollisionModel()->DisallowCollisionsWith(2);
    sys->AddBody(bottom);

    // Mold (body 1): a ring of plates standing on the bottom platen, open at both ends
    auto mold = utils::CreateCylindricalContainerFromBoxes(sys, mat_mold, R, h_mold, p.wall, p.plates, ChVector3d(0, 0, 0), z2y, true, false, false, true, true, false);
    mold->GetCollisionModel()->SetFamily(2);
    mold->GetCollisionModel()->DisallowCollisionsWith(1);

    // Top platen (body 2), parked above the mold until it is placed, then moved by the program
    double y_park = h_mold + p.h_platen + 0.05;
    auto top = chrono_types::make_shared<ChBody>();
    top->SetPos(ChVector3d(0, y_park, 0));
    top->SetFixed(true);
    top->EnableCollision(true);
    utils::AddCylinderGeometry(top.get(), mat_platen, R_platen, p.h_platen, ChVector3d(0, 0, 0), z2y, true);
    top->GetCollisionModel()->SetFamily(1);
    top->GetCollisionModel()->DisallowCollisionsWith(2);
    sys->AddBody(top);

    // Spheres, each with a force accumulator for the membrane
    std::vector<std::shared_ptr<ChBody>> spheres;
    int first_particle = (int)sys->GetBodies().size();
    for (int i = 0; i < p.num_spheres; i++) {
        double r = radius_of[i];
        double mass = p.density * (4.0 / 3) * CH_PI * r * r * r;
        auto sphere = chrono_types::make_shared<ChBody>();
        sphere->SetMass(mass);
        sphere->SetInertiaXX(0.4 * mass * r * r * ChVector3d(1, 1, 1));
        sphere->SetPos(positions[i]);
        sphere->EnableCollision(true);
        auto ct_shape = chrono_types::make_shared<ChCollisionShapeSphere>(mat_sphere, r);
        sphere->AddCollisionShape(ct_shape);
        sphere->GetCollisionModel()->SetFamily(3);
        auto vis_shape = chrono_types::make_shared<ChVisualShapeSphere>(r);
        vis_shape->SetColor(ChColor(0.8f, 0.8f, 0.8f));
        sphere->AddVisualShape(vis_shape);
        sphere->AddAccumulator();
        sys->AddBody(sphere);
        spheres.push_back(sphere);
    }

    Membrane membrane(sys, spheres, radius_of, first_particle, R, d_mean, p.bin_size * d_mean, p.membrane == "bins" ? Membrane::Model::BINS : Membrane::Model::VORONOI);

    if (!p.render.empty()) {
#ifdef CHRONO_VSG
        return RenderStates(p, sys, bottom, mold, top, spheres, R, d_mean);
#else
        std::cerr << "Rendering needs Chrono::VSG" << endl;
        return 1;
#endif
    }

    // Phase and reference values
    // phase 0: fill, 1: consolidate (top platen placed, membrane on), 2: test

    int phase = 0;
    bool mold_lifted = false;
    double y_face = NAN;   // face of the top platen (its bottom)
    double v_top = 0;      // vertical speed of the top platen
    double t_start = NAN;  // start of the test
    double H0 = NAN;       // height at the start of the test
    double V0 = NAN;       // volume enclosed by the membrane at the start of the test
    double R0 = NAN;       // mean radius at the start of the test

    auto test_end = [&](double t0) { return t0 + (triaxial ? p.strain_end / p.strain_rate : path.Duration()); };
    double end_time = test_end(p.t_test);

    // Cell pressure: none while filling; rising from zero to the consolidation pressure (sigma3, or the first pressure
    // of the path) between t_place and t_mold, then held; during the test, sigma3 or the path
    double p_consolidate = triaxial ? p.sigma3 : p.path[0];
    auto cell_pressure = [&](double time) {
        if (phase == 0)
            return 0.0;
        if (phase == 1)
            return p_consolidate * std::min(1.0, (time - p.t_place) / (p.t_mold - p.t_place));
        return triaxial ? p.sigma3 : path.Pressure(time - t_start);
    };

    auto place_top = [&]() {
        double y_max = 0;
        for (size_t i = 0; i < spheres.size(); i++)
            y_max = std::max(y_max, spheres[i]->GetPos().y() + radius_of[i]);
        y_face = y_max + 0.0005;
        v_top = 0;
        top->SetPos(ChVector3d(0, y_face + p.h_platen / 2, 0));
    };

    // The mold leaves contact once the membrane holds the sample; it is then lifted clear, its bottom above the top
    // platen (a body cannot leave the Multicore collision system)
    auto lift_mold = [&]() {
        mold->SetPos(ChVector3d(0, y_face + p.h_platen + 0.01, 0));
        mold_lifted = true;
    };

    auto save = [&](const std::string& file) { return SaveState(*sys, file, {(double)phase, t_start, H0, V0, R0, solid_volume}); };

    // Start from a settled sample, or resume from a saved state
    if (!p.load_fill.empty() || !p.load_state.empty()) {
        const std::string& file = p.load_state.empty() ? p.load_fill : p.load_state;
        std::vector<double> extra;
        if (!LoadState(*sys, file, extra) || extra.size() != kStateExtra || std::abs(extra[5] - solid_volume) > 1e-9 * solid_volume) {
            std::cerr << "Cannot load state from " << file << " (missing file, or different spheres)" << endl;
            return 1;
        }
        phase = (int)extra[0];
        if (phase >= 1) {
            y_face = top->GetPos().y() - p.h_platen / 2;
            v_top = top->GetPosDt().y();
            mold_lifted = mold->GetPos().y() > 1e-6;
            mat_sphere->SetFriction((float)p.mu);
        }
        if (phase == 2) {
            t_start = extra[1];
            H0 = extra[2];
            V0 = extra[3];
            R0 = extra[4];
            end_time = test_end(t_start);
            sys->GetSettings()->solver.tangential_displ_mode = HistoryModel(p.history_test);
        }
        cout << "Loaded state at t = " << sys->GetChTime() << " s (phase " << phase << ") from " << file << endl;
    }

    // Run metadata (rewritten at the start of the test and at the end)

    auto write_meta = [&](bool finished, double wall_time, unsigned long long overflows, const std::string& status) {
        JsonObject m;
        m.Add("program", "test_VDEM_triaxial_membrane");
        m.Add("description",
              "Drained triaxial test with a stress-controlled membrane, Chrono::Multicore SMC (Hertz), Y up along the specimen axis, SI units, compression positive");
        m.Add("chrono_version", CHRONO_VERSION);
        m.Add("start", start_date);
        m.Add("status", status);
        m.Add("radii", p.radii);
        m.Add("fractions", p.fractions);
        std::vector<double> counts(size_counts.begin(), size_counts.end());
        m.Add("counts", counts);
        m.Add("num_spheres", p.num_spheres);
        m.Add("mean_diameter", d_mean);
        m.Add("solid_volume", solid_volume);
        m.Add("density", p.density);
        m.Add("E", p.E);
        m.Add("nu", p.nu);
        m.Add("cor", p.cor);
        m.Add("mu_sphere", p.mu);
        m.Add("mu_platen", p.mu_platen);
        m.Add("mu_fill", p.mu_fill);
        m.Add("friction_composition", "max");
        m.Add("diameter", p.diameter);
        m.Add("h_platen", p.h_platen);
        m.Add("platen_diameter", p.platen_diameter);
        m.Add("mold_wall", p.wall);
        m.Add("mold_plates", p.plates);
        m.Add("mold_height", h_mold);
        m.Add("gen_fraction", p.gen_fraction);
        m.Add("gen_height", h_gen);
        m.Add("membrane_model", p.membrane);
        m.Add("membrane_every", p.membrane_every);
        m.Add("bin_size", p.bin_size * d_mean);
        m.Add("test", p.test);
        m.Add("sigma3", p.sigma3);
        m.Add("strain_rate", p.strain_rate);
        m.Add("strain_end", p.strain_end);
        m.Add("path", p.path);
        m.Add("ramp", p.ramp);
        m.Add("hold", p.hold);
        m.Add("gravity", p.gravity);
        m.Add("t_place", p.t_place);
        m.Add("t_mold", p.t_mold);
        m.Add("t_test", p.t_test);
        m.Add("t_start", t_start);
        m.Add("t_end", end_time);
        m.Add("height_test_start", H0);
        m.Add("volume_test_start", V0);
        m.Add("radius_test_start", R0);
        m.Add("void_ratio_test_start", V0 / solid_volume - 1);
        m.Add("servo_gain", p.servo_gain);
        m.Add("servo_vmax", p.servo_vmax);
        m.Add("history_prepare", p.history);
        m.Add("history_test", p.history_test);
        m.Add("contact_force_model", "Hertz");
        m.Add("narrowphase", "HYBRID");
        m.Add("dt", p.dt);
        m.Add("threads", threads);
        m.Add("seed", p.seed);
        m.Add("out_dt", p.out_dt);
        m.Add("snapshots_axial_strain", p.snaps);
        m.Add("frame_strain", p.frame_strain);
        m.Add("save_fill", p.save_fill);
        m.Add("load_fill", p.load_fill);
        m.Add("save_state", p.save_state);
        m.Add("load_state", p.load_state);
        m.Add("snap_states", p.snap_states);
        m.Add("frame_states", p.frame_states);
        m.Add("first_particle_index", first_particle);
        m.AddRaw("walls", "{\"bottom platen\": {\"body\": 0}, \"mold\": {\"body\": 1}, \"top platen\": {\"body\": 2}}");
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
    ts << "# Triaxial test with a stress-controlled membrane. SI units; stresses and strains positive in compression. Y up "
          "along the specimen axis.\n";
    ts << "# phase: 0 fill, 1 consolidate, 2 test. Forces are contact forces of the spheres on the named body.\n";
    ts << "# height: between the platen faces; volume: enclosed by the membrane; area = volume/height; radius = sqrt(area/pi).\n";
    ts << "# sigma1_top = Fy(top platen)/area, sigma1_bottom = -Fy(bottom platen)/area, sigma1 = their mean; sigma3: the "
          "membrane pressure.\n";
    ts << "# p = (sigma1 + 2 sigma3)/3, q = sigma1 - sigma3, eta = q/p, sin_phi = (sigma1 - sigma3)/(sigma1 + sigma3).\n";
    ts << "# membrane_spheres, membrane_force, membrane_bins (bins, or raster cells of the voronoi model), outside (spheres "
          "beyond the platen faces): at the last membrane update.\n";
    ts << "# sigma1_top_a0 = Fy(top platen)/A0, with A0 the area at the start of the test (no area correction, as in Cui et al.'s "
          "laboratory tests); sin_phi_a0 from it.\n";
    ts << "# eps_a, eps_r, eps_v: from the height, radius and volume at the start of the test. e: from the membrane volume; "
          "e_cv: in a central cylinder of half the radius and half the height.\n";
    ts << "time,phase,eps_a,eps_r,eps_v,height,radius,volume,e,e_cv,"
          "sigma1_top,sigma1_bottom,sigma1,sigma3,p,q,eta,sin_phi,sigma1_top_a0,sin_phi_a0,"
          "top_y,top_vy,top_fx,top_fy,top_fz,bottom_fx,bottom_fy,bottom_fz,mold_f,"
          "membrane_spheres,membrane_force,membrane_bins,outside,kinetic_energy,"
          "contacts,contacts_pp,contacts_pw,coord,coord_mech,rattlers,frac_sliding,mean_fn,"
          "overflows_step,overflows_total,wall_clock,step_time,"
          "t_collision,t_broad,t_narrow,t_update,t_advance,t_contacts,t_ls_setup,t_ls_solve\n";
    ts << std::setprecision(9);

    std::ofstream frames_index(out_dir + "/frames.csv");
    frames_index << "frame,time,eps_a,phase\n";

#ifdef CHRONO_VSG
    std::shared_ptr<vsg3d::ChVisualSystemVSG> vis;
    if (p.vsg) {
        auto mcore_vis = chrono_types::make_shared<ChMulticoreVisualizationVSG>(sys);
        mcore_vis->SetNetworkColorRamp(ChColor(0.957f, 0.694f, 0.514f), ChColor(0.518f, 0.235f, 0.047f));
        mcore_vis->SetBodySelection(spheres);

        vis = chrono_types::make_shared<vsg3d::ChVisualSystemVSG>();
        vis->AttachSystem(sys);
        vis->AttachPlugin(mcore_vis);
        vis->SetWindowTitle("Triaxial test");
        vis->SetWindowSize(1280, 720);
        vis->SetCameraVertical(CameraVerticalDir::Y);
        vis->AddCamera(ChVector3d(0, 0.5 * h_gen, -4 * p.diameter), ChVector3d(0, 0.4 * h_gen, 0));
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
    int membrane_phase = -1;  // phase of the last membrane update
    // The Voronoi rule reads the contacts of the last step, which a state loaded from a file does not have: the first
    // membrane update then waits for the collision detection of the first step
    bool contacts_known = p.load_fill.empty() && p.load_state.empty();
    double last_wall = 0;
    long long last_step = 0;
    long long step = 0;
    unsigned long long overflows_total = 0;

    auto axial_strain = [&]() { return phase == 2 ? (H0 - y_face) / H0 : 0.0; };

    // A run resumed from a saved state skips the snapshots and frames it has passed
    if (phase == 2) {
        if (triaxial) {
            double e0 = axial_strain() - 0.5 * p.strain_rate * p.dt;
            while (next_snap < p.snaps.size() && p.snaps[next_snap] < e0)
                next_snap++;
            if (p.frame_strain > 0)
                frame = (int)std::ceil(e0 / p.frame_strain);
        } else {
            while (next_snap + 1 < p.path.size() && sys->GetChTime() - t_start >= path.SegmentEnd(next_snap) - 0.5 * p.dt)
                next_snap++;
        }
    }

    auto write_snapshot = [&](const std::string& name) {
        auto contacts = GetContacts(*sys);
        std::vector<int> counts;
        ComputeContactStats(contacts, first_particle, p.num_spheres, &counts);
        std::ostringstream comment;
        comment << "time = " << std::setprecision(8) << sys->GetChTime() << " s, axial strain = " << axial_strain() << ", cell pressure = " << cell_pressure(sys->GetChTime())
                << " Pa";
        WriteParticles(particles_dir + "/" + name + ".csv", spheres, radius_of, counts, comment.str());
        WriteContacts(contacts_dir + "/" + name + ".csv", contacts, comment.str());
    };

    while (sys->GetChTime() < end_time - 1e-9) {
        double time = sys->GetChTime();

        // Consolidate: save the settled sample, place the top platen and switch the membrane on
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
            place_top();
            mat_sphere->SetFriction((float)p.mu);
            phase = 1;
        }

        if (phase >= 1 && !mold_lifted && time >= p.t_mold - 1e-9)
            lift_mold();

        // Test: from the consolidated sample
        if (phase == 1 && time >= p.t_test - 1e-9) {
            t_start = time;
            H0 = y_face;
            V0 = membrane.volume;
            R0 = std::sqrt(V0 / (CH_PI * H0));
            end_time = test_end(t_start);
            phase = 2;
            if (!p.save_state.empty() && p.load_state.empty()) {
                if (save(p.save_state))
                    cout << "Saved state at t = " << time << " s to " << p.save_state << endl;
                else
                    std::cerr << "Cannot save state to " << p.save_state << endl;
            }
            sys->GetSettings()->solver.tangential_displ_mode = HistoryModel(p.history_test);
            write_meta(false, 0, overflows_total, "running");
        }

        // Top platen: the stress servo while consolidating and in the isotropic test, a constant speed in the triaxial
        // test. Then the membrane, between the platen faces: updated every membrane_every steps and when a phase starts.
        double pressure = cell_pressure(time);
        if (phase >= 1) {
            if (phase == 2 && triaxial) {
                v_top = -p.strain_rate * H0;
            } else {
                double area = membrane.volume > 0 ? membrane.volume / y_face : CH_PI * R * R;
                const auto& map = sys->data_manager->host_data.ct_body_map;
                double f_top = map.size() > top->GetIndex() ? sys->GetBodyContactForce(top).y : 0.0;
                v_top = std::clamp(-p.servo_gain * (pressure - f_top / area), -p.servo_vmax, p.servo_vmax);
            }
            y_face += v_top * p.dt;
            top->SetPos(ChVector3d(0, y_face + p.h_platen / 2, 0));
            top->SetLinVel(ChVector3d(0, v_top, 0));
            if (contacts_known && (step % p.membrane_every == 0 || phase != membrane_phase)) {
                membrane.Apply(pressure, 0, y_face);
                membrane_phase = phase;
            }
        }

        sys->DoStepDynamics(p.dt);
        step++;
        if (!contacts_known) {
            contacts_known = true;
            if (phase >= 1) {
                membrane.Apply(cell_pressure(sys->GetChTime()), 0, y_face);
                membrane_phase = phase;
            }
        }
        overflows_total = sys->GetTotalContactHistoryOverflows();

#ifdef CHRONO_VSG
        if (vis && step % 200 == 0) {
            if (!vis->Run())
                break;
            vis->Render();
        }
#endif

        time = sys->GetChTime();

        // Snapshots: at the requested axial strains and frames (triaxial), at the end of each hold (isotropic)
        if (phase == 2 && triaxial) {
            double eps = axial_strain();
            while (next_snap < p.snaps.size() && eps >= p.snaps[next_snap] - 0.5 * p.strain_rate * p.dt) {
                char name[64];
                snprintf(name, sizeof(name), "snap_e%05.2fpc", 100 * p.snaps[next_snap]);
                write_snapshot(name);
                if (p.snap_states && !save(states_dir + "/" + name + ".bin"))
                    std::cerr << "Cannot save state to " << states_dir + "/" + name + ".bin" << endl;
                next_snap++;
            }
            if (p.frame_strain > 0 && eps >= frame * p.frame_strain - 0.5 * p.strain_rate * p.dt) {
                char name[64];
                snprintf(name, sizeof(name), "frame_%05d", frame);
                write_snapshot(name);
                if (p.frame_states && !save(states_dir + "/" + name + ".bin"))
                    std::cerr << "Cannot save state to " << states_dir + "/" + name + ".bin" << endl;
                frames_index << frame << "," << time << "," << eps << ",2" << std::endl;
                frame++;
            }
        } else if (phase == 2) {
            while (next_snap + 1 < p.path.size() && time - t_start >= path.SegmentEnd(next_snap) - 0.5 * p.dt) {
                char name[64];
                snprintf(name, sizeof(name), "snap_k%02d_p%06.1fkPa", (int)next_snap + 1, 1e-3 * p.path[next_snap + 1]);
                write_snapshot(name);
                if (p.snap_states && !save(states_dir + "/" + name + ".bin"))
                    std::cerr << "Cannot save state to " << states_dir + "/" + name + ".bin" << endl;
                next_snap++;
            }
        }

        // Time series
        if (time >= data_out_frame * p.out_dt) {
            auto contacts = GetContacts(*sys);
            auto st = ComputeContactStats(contacts, first_particle, p.num_spheres);

            // Contact forces on the platens and the mold
            ChVector3d f_bottom(0), f_mold(0), f_top(0);
            for (const auto& c : contacts) {
                // Force on body 2 is F; on body 1 it is -F
                ChVector3d F = c.fn * c.n + c.ft;
                auto add = [&](int b, const ChVector3d& f) {
                    if (b == 0)
                        f_bottom += f;
                    else if (b == 1)
                        f_mold += f;
                    else if (b == 2)
                        f_top += f;
                };
                add(c.b1, -F);
                add(c.b2, F);
            }

            double height = phase >= 1 ? y_face : NAN;
            double volume = phase >= 1 ? membrane.volume : NAN;
            double area = volume / height;
            double radius = std::sqrt(area / CH_PI);
            double s3 = cell_pressure(time);
            double s1_top = f_top.y() / area;
            double s1_bottom = -f_bottom.y() / area;
            double s1 = 0.5 * (s1_top + s1_bottom);
            double pm = (s1 + 2 * s3) / 3;
            double q = s1 - s3;
            double eps_a = axial_strain();
            double eps_r = phase == 2 ? (R0 - radius) / R0 : 0.0;
            double eps_v = phase == 2 ? (V0 - volume) / V0 : 0.0;
            double e = volume / solid_volume - 1;
            double s1_top_a0 = phase == 2 ? f_top.y() * H0 / V0 : NAN;

            // Void ratio in a central cylinder of half the radius and half the height (spheres counted by their centers)
            double e_cv = NAN;
            if (phase >= 1) {
                double r_cv = 0.5 * R;
                double y_lo = 0.25 * height;
                double y_hi = 0.75 * height;
                double v_solid = 0;
                for (size_t i = 0; i < spheres.size(); i++) {
                    const auto& c = spheres[i]->GetPos();
                    if (c.y() > y_lo && c.y() < y_hi && c.x() * c.x() + c.z() * c.z() < r_cv * r_cv)
                        v_solid += (4.0 / 3) * CH_PI * radius_of[i] * radius_of[i] * radius_of[i];
                }
                e_cv = v_solid > 0 ? CH_PI * r_cv * r_cv * (y_hi - y_lo) / v_solid - 1 : NAN;
            }

            double ke = 0;
            for (size_t i = 0; i < spheres.size(); i++) {
                const auto& s = spheres[i];
                ke += 0.5 * s->GetMass() * s->GetPosDt().Length2();
                ke += 0.5 * 0.4 * s->GetMass() * radius_of[i] * radius_of[i] * s->GetAngVelParent().Length2();
            }

            double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall_start).count();
            double step_time = (step > last_step) ? (wall - last_wall) / (step - last_step) : 0;
            last_wall = wall;
            last_step = step;

            ts << time << "," << phase << "," << eps_a << "," << eps_r << "," << eps_v << "," << height << "," << radius << "," << volume << "," << e << "," << e_cv << ","
               << s1_top << "," << s1_bottom << "," << s1 << "," << s3 << "," << pm << "," << q << "," << q / pm << "," << (s1 - s3) / (s1 + s3) << "," << s1_top_a0 << ","
               << (s1_top_a0 - s3) / (s1_top_a0 + s3) << "," << top->GetPos().y() << "," << v_top << "," << f_top.x() << "," << f_top.y() << "," << f_top.z() << "," << f_bottom.x()
               << "," << f_bottom.y() << "," << f_bottom.z() << "," << f_mold.Length() << "," << membrane.num_carriers << "," << membrane.force << "," << membrane.num_bins << ","
               << membrane.num_outside << "," << ke << "," << st.num_contacts << "," << st.num_pp << "," << st.num_pw << "," << st.coord << "," << st.coord_mech << ","
               << st.num_rattlers << "," << st.frac_sliding << "," << st.mean_fn << "," << sys->GetNumContactHistoryOverflows() << "," << overflows_total << "," << wall << ","
               << step_time << "," << sys->GetTimerCollision() << "," << sys->GetTimerCollisionBroad() << "," << sys->GetTimerCollisionNarrow() << "," << sys->GetTimerUpdate()
               << "," << sys->GetTimerAdvance() << "," << sys->GetTimerProcessContact() << "," << sys->GetTimerLSsetup() << "," << sys->GetTimerLSsolve() << std::endl;

            if (data_out_frame % 10 == 0) {
                cout << std::fixed << std::setprecision(3) << "t = " << time << "  phase " << phase << "  eps_a = " << std::setprecision(2) << 100 * eps_a
                     << " %  sigma1 = " << std::setprecision(1) << 1e-3 * s1 << " kPa  sigma3 = " << 1e-3 * s3 << " kPa  q/p = " << std::setprecision(3) << q / pm << "  e = " << e
                     << "  contacts = " << st.num_contacts << "  step = " << std::setprecision(2) << 1e3 * step_time << " ms  overflows = " << overflows_total << endl;
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

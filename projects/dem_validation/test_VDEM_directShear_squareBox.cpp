// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All right reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Author: Jonathan Fleischmann, Radu Serban
// =============================================================================
//
// Soft-sphere (SMC) direct shear test in a square box, with Chrono::Multicore.
//
// Two samples are predefined (--sample):
// - hartl   : the 2015 model of the Hartl and Ooi (2008) direct shear tests on 6 mm glass beads: 5000 beads in a
//             12 cm square box (default).
// - history : the 1800-sphere sample of the 2015 contact-history study (Fleischmann et al., JCND 2016):
//             5 mm spheres in a 6 cm square box, compacted with friction 0.05.
// Every parameter of a sample can be changed on the command line.
//
// The run has three phases:
// - fill    : spheres laid on a jittered grid fall into the box; the lid is parked above the box.
// - compact : the lid (mass = sigma A / g) is placed in the upper box and released; a prismatic joint lets it move
//             only vertically. With --dense, the sphere friction is lowered during this phase.
// - shear   : the lower box moves at constant speed in +z (position and velocity imposed); the upper box is fixed.
//
// The global reference frame has Y up. All units SI. Run-time visualization with Chrono::VSG is enabled with --vsg.
//
// The time series also records the timers of the last step (s): collision (broad and narrow phase), update, and
// advance (contact forces, linear solver setup and solve).
//
// Outputs, in the output directory:
//   meta.json           all parameters, run information and summary
//   timeseries.csv      one row per output interval (written as the run goes)
//   particles/*.csv     particle snapshots (at the shear displacements in --snaps, and frames)
//   contacts/*.csv      contact snapshots, at the same times
//   state_shear.bin     the prepared sample at the start of shearing (with --save_state)
// =============================================================================

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "chrono/ChConfig.h"
#include "chrono/ChVersion.h"
#include "chrono/assets/ChVisualShapeBox.h"
#include "chrono/assets/ChVisualShapeSphere.h"
#include "chrono/core/ChDataPath.h"
#include "chrono/physics/ChLinkLock.h"

#include "chrono_multicore/physics/ChSystemMulticore.h"

#include "chrono_thirdparty/cxxopts/ChCLI.h"

#ifdef CHRONO_VSG
    #include "chrono_vsg/ChVisualSystemVSG.h"
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
    std::string sample;

    // Spheres
    int layers, nx, nz;  // grid of falling spheres (layers in y, nx in x, nz in z)
    double radius;       // m
    double density;      // kg/m^3
    double E;            // Young's modulus (spheres and walls), Pa
    double nu;           // Poisson ratio
    double cor;          // coefficient of restitution
    double mu;           // sphere friction
    double mu_wall;      // wall friction (the composite friction is the smaller of the two)
    bool dense;          // lower the sphere friction to mu_dense while the lid compacts the sample
    double mu_dense;

    // Box
    double width;      // inner width (x) and length (z) of the box, m
    double height;     // height of the box walls, m
    double thickness;  // wall thickness, m

    // Loading
    double sigma;               // normal stress, Pa
    double gravity;             // m/s^2
    double speed;               // shear speed, m/s
    double travel;              // shear travel, m
    double t_settle;            // end of the fill phase, s
    double t_shear;             // start of shearing, s
    std::string history;        // tangential displacement model while preparing: multistep, onestep, none
    std::string history_shear;  // model while shearing ("same": as while preparing)

    // Numerics
    double dt;
    int threads;
    int seed;

    // Output
    std::string out_dir;
    double out_dt;              // time series interval, s
    std::vector<double> snaps;  // shear displacements of the snapshots, m
    double frame_disp;          // shear displacement between frames, m (0: none)
    double frame_dt_prep;       // time between frames before shearing, s (0: none)
    std::string save_state;     // file for the state at the start of shearing
    std::string load_state;     // start from a saved state (at the start of shearing)
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
    ChCLI cli(argv[0], "Direct shear test in a square box (Chrono::Multicore, SMC)");

    cli.AddOption<std::string>("Sample", "sample", "Predefined sample: hartl or history", "hartl");
    cli.AddOption<int>("Sample", "layers", "Number of sphere layers in the falling grid");
    cli.AddOption<int>("Sample", "nx", "Spheres per row in x");
    cli.AddOption<int>("Sample", "nz", "Spheres per row in z");
    cli.AddOption<double>("Sample", "radius", "Sphere radius (m)");
    cli.AddOption<double>("Sample", "density", "Sphere density (kg/m3)");
    cli.AddOption<double>("Sample", "E", "Young's modulus of spheres and walls (Pa)");
    cli.AddOption<double>("Sample", "nu", "Poisson ratio");
    cli.AddOption<double>("Sample", "cor", "Coefficient of restitution");
    cli.AddOption<double>("Sample", "mu", "Sphere friction coefficient");
    cli.AddOption<double>("Sample", "mu_wall", "Wall friction coefficient");
    cli.AddOption<bool>("Sample", "dense", "Lower the sphere friction while the lid compacts the sample");
    cli.AddOption<double>("Sample", "mu_dense", "Sphere friction while compacting (with --dense)");
    cli.AddOption<double>("Sample", "width", "Inner width and length of the box (m)");
    cli.AddOption<double>("Sample", "height", "Height of the box walls (m)");

    cli.AddOption<double>("Loading", "sigma", "Normal stress (Pa)");
    cli.AddOption<double>("Loading", "gravity", "Gravitational acceleration (m/s2)", "9.81");
    cli.AddOption<double>("Loading", "speed", "Shear speed (m/s)", "0.001");
    cli.AddOption<double>("Loading", "travel", "Shear travel (m)", "0.01");
    cli.AddOption<double>("Loading", "t_settle", "End of the fill phase (s)");
    cli.AddOption<double>("Loading", "t_shear", "Start of shearing (s)");
    cli.AddOption<std::string>("Loading", "history", "Tangential displacement model: multistep, onestep or none", "multistep");
    cli.AddOption<std::string>("Loading", "history_shear", "Tangential displacement model while shearing (same: as --history)", "same");

    cli.AddOption<double>("Numerics", "dt", "Time step (s)", "1e-5");
    cli.AddOption<int>("Numerics", "threads", "Number of OpenMP threads", "20");
    cli.AddOption<int>("Numerics", "seed", "Seed of the grid jitter (1 reproduces the 2015 program)", "1");

    cli.AddOption<std::string>("Output", "out", "Output directory", "directShear_squareBox");
    cli.AddOption<double>("Output", "out_dt", "Time series interval (s)", "0.01");
    cli.AddOption<std::vector<double>>("Output", "snaps", "Shear displacements of the snapshots (mm)", "0,0.5,1,2,4,8,12");
    cli.AddOption<double>("Output", "frame_disp", "Shear displacement between particle and contact frames (mm, 0: none)", "0");
    cli.AddOption<double>("Output", "frame_dt_prep", "Time between particle and contact frames before shearing (s, 0: none)", "0");
    cli.AddOption<std::string>("Output", "save_state", "Save the sample at the start of shearing to this file", "");
    cli.AddOption<std::string>("Output", "load_state", "Start shearing from a sample saved with --save_state", "");
    cli.AddOption<bool>("Output", "vsg", "Run-time visualization with Chrono::VSG");
    cli.AddOption<std::string>("Output", "provenance", "JSON text recorded in meta.json", "");

    if (!cli.Parse(argc, argv, true))
        return false;

    p.sample = cli.GetAsType<std::string>("sample");
    if (p.sample == "hartl") {
        // test_VDEM_demo_shear.cpp (2015-04-22): Hartl and Ooi's glass beads
        p.layers = 50;
        p.nx = 10;
        p.nz = 10;
        p.radius = 0.003;
        p.density = 2550;
        p.E = 4.0e7;
        p.nu = 0.22;
        p.cor = 0.87;
        p.mu = 0.18;
        p.mu_wall = 0.13;
        p.dense = false;
        p.mu_dense = 0.01;
        p.width = 0.12;
        p.height = 0.06;
        p.sigma = 3.1e3;
        p.t_settle = 0.23;
        p.t_shear = 2.0;
    } else if (p.sample == "history") {
        // demo_shearDEM.cpp (2015-04-13): the contact-history study
        p.layers = 50;
        p.nx = 6;
        p.nz = 6;
        p.radius = 0.0025;
        p.density = 2500;
        p.E = 8.0e6;
        p.nu = 0.3;
        p.cor = 0.9;
        p.mu = 0.5;
        p.mu_wall = 0.5;
        p.dense = true;
        p.mu_dense = 0.05;
        p.width = 0.06;
        p.height = 0.06;
        p.sigma = 1.0e3;
        p.t_settle = 0.2;
        p.t_shear = 10.0;
    } else {
        std::cerr << "Unknown sample: " << p.sample << endl;
        return false;
    }
    p.thickness = 0.01;

    auto set_int = [&](const char* name, int& v) {
        if (cli.CheckOption(name))
            v = cli.GetAsType<int>(name);
    };
    auto set_dbl = [&](const char* name, double& v) {
        if (cli.CheckOption(name))
            v = cli.GetAsType<double>(name);
    };
    set_int("layers", p.layers);
    set_int("nx", p.nx);
    set_int("nz", p.nz);
    set_dbl("radius", p.radius);
    set_dbl("density", p.density);
    set_dbl("E", p.E);
    set_dbl("nu", p.nu);
    set_dbl("cor", p.cor);
    set_dbl("mu", p.mu);
    set_dbl("mu_wall", p.mu_wall);
    if (cli.CheckOption("dense"))
        p.dense = true;
    set_dbl("mu_dense", p.mu_dense);
    set_dbl("width", p.width);
    set_dbl("height", p.height);
    set_dbl("sigma", p.sigma);
    set_dbl("t_settle", p.t_settle);
    set_dbl("t_shear", p.t_shear);

    p.gravity = cli.GetAsType<double>("gravity");
    p.speed = cli.GetAsType<double>("speed");
    p.travel = cli.GetAsType<double>("travel");
    p.history = cli.GetAsType<std::string>("history");
    p.history_shear = cli.GetAsType<std::string>("history_shear");
    if (p.history_shear == "same")
        p.history_shear = p.history;
    HistoryModel(p.history);
    HistoryModel(p.history_shear);

    p.dt = cli.GetAsType<double>("dt");
    p.threads = cli.GetAsType<int>("threads");
    p.seed = cli.GetAsType<int>("seed");

    p.out_dir = cli.GetAsType<std::string>("out");
    p.out_dt = cli.GetAsType<double>("out_dt");
    p.snaps = cli.GetAsType<std::vector<double>>("snaps");
    for (auto& s : p.snaps)
        s *= 1e-3;
    p.frame_disp = 1e-3 * cli.GetAsType<double>("frame_disp");
    p.frame_dt_prep = cli.GetAsType<double>("frame_dt_prep");
    p.save_state = cli.GetAsType<std::string>("save_state");
    p.load_state = cli.GetAsType<std::string>("load_state");
    p.vsg = cli.GetAsType<bool>("vsg");
    p.provenance = cli.GetAsType<std::string>("provenance");

    return true;
}

// -----------------------------------------------------------------------------
// Utility for adding (visible or invisible) walls, given half-dimensions
// -----------------------------------------------------------------------------

void AddWall(std::shared_ptr<ChBody>& body,
             std::shared_ptr<ChContactMaterial> mat,
             const ChVector3d& dim,
             const ChVector3d& loc,
             bool visible) {
    auto ct_shape = chrono_types::make_shared<ChCollisionShapeBox>(mat, 2.0 * dim);
    body->AddCollisionShape(ct_shape, ChFrame<>(loc, QUNIT));

    if (visible) {
        auto box = chrono_types::make_shared<ChVisualShapeBox>(2.0 * dim);
        box->SetColor(ChColor(0.75f, 0.75f, 0.75f));
        body->AddVisualShape(box, ChFrame<>(loc));
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
    for (const auto& dir : {out_dir, particles_dir, contacts_dir}) {
        std::error_code ec;
        std::filesystem::create_directories(dir, ec);
        if (!std::filesystem::is_directory(dir)) {
            std::cerr << "Error creating directory " << dir << endl;
            return 1;
        }
    }

    auto wall_start = std::chrono::steady_clock::now();
    std::string start_date = Now();

    double gravity = p.gravity;
    double radius = p.radius;
    double mass = p.density * (4.0 / 3) * CH_PI * radius * radius * radius;
    double width = p.width;
    double length = p.width;
    double height = p.height;
    double thickness = p.thickness;
    double area = width * length;
    double end_time = p.t_shear + p.travel / p.speed;

    // Define a rotation of -90 degrees around x (z2y)
    ChQuaternion<> z2y;
    z2y.SetFromAngleAxis(-CH_PI / 2, ChVector3d(1, 0, 0));

    // Create the system

    ChSystemMulticoreSMC* sys = new ChSystemMulticoreSMC();
    sys->SetCollisionSystemType(ChCollisionSystem::Type::MULTICORE);
    sys->SetGravitationalAcceleration(ChVector3d(0, -gravity, 0));

    int threads = std::min(p.threads, omp_get_num_procs());
    sys->SetNumThreads(threads);

    sys->GetSettings()->solver.use_full_inertia_tensor = false;
    sys->GetSettings()->solver.tolerance = 0.01;
    sys->GetSettings()->solver.max_iteration_bilateral = 100;
    sys->GetSettings()->solver.clamp_bilaterals = false;
    sys->GetSettings()->solver.bilateral_clamp_speed = 0.1;

    sys->GetSettings()->solver.contact_force_model = ChSystemSMC::ContactForceModel::Hertz;
    sys->GetSettings()->solver.tangential_displ_mode = HistoryModel(p.history);

    sys->GetSettings()->collision.bins_per_axis = vec3(10, 10, 10);
    sys->GetSettings()->collision.narrowphase_algorithm = ChNarrowphase::Algorithm::HYBRID;

    // Materials: one for the spheres, one for the box walls and the lid

    auto material = chrono_types::make_shared<ChContactMaterialSMC>();
    material->SetYoungModulus((float)p.E);
    material->SetPoissonRatio((float)p.nu);
    material->SetRestitution((float)p.cor);
    material->SetFriction((float)p.mu);

    auto mat_ext = chrono_types::make_shared<ChContactMaterialSMC>();
    mat_ext->SetYoungModulus((float)p.E);
    mat_ext->SetPoissonRatio((float)p.nu);
    mat_ext->SetRestitution((float)p.cor);
    mat_ext->SetFriction((float)p.mu_wall);

    // Lower bin: floor (shape 0) at mid-height of its walls (shapes 1 to 4)

    auto bin = chrono_types::make_shared<ChBody>();
    bin->SetMass(1);
    bin->SetPos(ChVector3d(0, -height / 2, 0));
    bin->SetFixed(true);
    bin->EnableCollision(true);

    AddWall(bin, mat_ext, ChVector3d(width / 2, thickness / 2, length / 2), ChVector3d(0, 0, 0), true);
    AddWall(bin, mat_ext, ChVector3d(thickness / 2, height / 2, length / 2 + thickness),
            ChVector3d(-width / 2 - thickness / 2, 0, 0), true);
    AddWall(bin, mat_ext, ChVector3d(thickness / 2, height / 2, length / 2 + thickness),
            ChVector3d(width / 2 + thickness / 2, 0, 0), false);
    AddWall(bin, mat_ext, ChVector3d(width / 2 + thickness, height / 2, thickness / 2),
            ChVector3d(0, 0, -length / 2 - thickness / 2), true);
    AddWall(bin, mat_ext, ChVector3d(width / 2 + thickness, height / 2, thickness / 2),
            ChVector3d(0, 0, length / 2 + thickness / 2), true);
    bin->GetCollisionModel()->SetFamily(1);
    bin->GetCollisionModel()->DisallowCollisionsWith(2);
    bin->GetCollisionModel()->DisallowCollisionsWith(3);
    bin->GetCollisionModel()->AllowCollisionsWith(4);

    sys->AddBody(bin);

    // Upper shear box (shapes 0 to 3), raised by one sphere radius above the lower bin

    auto box = chrono_types::make_shared<ChBody>();
    box->SetMass(1);
    box->SetPos(ChVector3d(0, height / 2 + radius, 0));
    box->SetFixed(true);
    box->EnableCollision(true);

    AddWall(box, mat_ext, ChVector3d(thickness / 2, height / 2, length / 2 + thickness),
            ChVector3d(-width / 2 - thickness / 2, 0, 0), true);
    AddWall(box, mat_ext, ChVector3d(thickness / 2, height / 2, length / 2 + thickness),
            ChVector3d(width / 2 + thickness / 2, 0, 0), false);
    AddWall(box, mat_ext, ChVector3d(width / 2 + thickness, height / 2, thickness / 2),
            ChVector3d(0, 0, -length / 2 - thickness / 2), true);
    AddWall(box, mat_ext, ChVector3d(width / 2 + thickness, height / 2, thickness / 2),
            ChVector3d(0, 0, length / 2 + thickness / 2), true);
    box->GetCollisionModel()->SetFamily(2);
    box->GetCollisionModel()->DisallowCollisionsWith(1);
    box->GetCollisionModel()->DisallowCollisionsWith(3);
    box->GetCollisionModel()->AllowCollisionsWith(4);

    sys->AddBody(box);

    // Lid (load plate), parked above the falling spheres

    auto plate = chrono_types::make_shared<ChBody>();
    plate->SetMass(p.sigma * area / gravity);
    plate->SetPos(ChVector3d(0, 2.0 * radius * float(p.layers) + thickness, 0));
    plate->SetFixed(true);
    plate->EnableCollision(true);

    AddWall(plate, mat_ext, ChVector3d(width / 2, thickness / 2, length / 2), ChVector3d(0, 0, 0), true);
    plate->GetCollisionModel()->SetFamily(3);
    plate->GetCollisionModel()->DisallowCollisionsWith(1);
    plate->GetCollisionModel()->DisallowCollisionsWith(2);
    plate->GetCollisionModel()->AllowCollisionsWith(4);

    sys->AddBody(plate);

    // Spheres: a grid of (layers x nx x nz), pitch 4r horizontally and 2r vertically, with a horizontal jitter
    // of up to 0.99 r from rand() (seed 1, the C default, reproduces the 2015 program)

    std::srand((unsigned int)p.seed);
    std::vector<std::shared_ptr<ChBody>> balls;
    std::vector<double> radii;
    int first_particle = (int)sys->GetBodies().size();

    for (int i = 0; i < p.layers; i++) {
        for (int j = 0; j < p.nx; j++) {
            for (int k = 0; k < p.nz; k++) {
                double ball_y = 2.0 * radius * float(i);
                double ball_x = 4.0 * radius * (float(j - p.nx / 2) + 0.5) + 0.99 * radius * (float(std::rand() % 100) / 50 - 1.0);
                double ball_z = 4.0 * radius * (float(k - p.nz / 2) + 0.5) + 0.99 * radius * (float(std::rand() % 100) / 50 - 1.0);

                auto ball = chrono_types::make_shared<ChBody>();
                ball->SetMass(mass);
                ball->SetInertiaXX((2.0 / 5.0) * mass * radius * radius * ChVector3d(1, 1, 1));
                ball->SetPos(ChVector3d(ball_x, ball_y, ball_z));
                ball->SetFixed(false);
                ball->EnableCollision(true);

                auto ball_ct_shape = chrono_types::make_shared<ChCollisionShapeSphere>(material, radius);
                ball->AddCollisionShape(ball_ct_shape);
                ball->GetCollisionModel()->SetFamily(4);

                auto sphere = chrono_types::make_shared<ChVisualShapeSphere>(radius);
                sphere->SetColor(ChColor(0.8f, 0.8f, 0.8f));
                ball->AddVisualShape(sphere);

                sys->AddBody(ball);
                balls.push_back(ball);
                radii.push_back(radius);
            }
        }
    }
    int num_particles = (int)balls.size();

    // Prismatic joint between the lid and the upper box (translation along the global Y axis)

    auto prismatic_plate_box = chrono_types::make_shared<ChLinkLockPrismatic>();
    prismatic_plate_box->SetName("prismatic_plate_box");
    prismatic_plate_box->Initialize(plate, box, ChFrame<>(ChVector3d(0, 0, 0), z2y));
    sys->AddLink(prismatic_plate_box);

    // Phase flags and shear reference values

    bool settling = true;
    bool shearing = false;
    double shear_Height = 0;  // lid height at the start of shearing
    double shear_Disp = 0;    // bin position at the start of shearing

    if (!p.load_state.empty()) {
        std::vector<double> extra;
        if (!LoadState(*sys, p.load_state, extra)) {
            std::cerr << "Cannot load state from " << p.load_state << " (missing file or different sample)" << endl;
            return 1;
        }
        plate->SetFixed(false);
        settling = false;
        cout << "Loaded state at t = " << sys->GetChTime() << " s from " << p.load_state << endl;
    }

    // Solid volume, for the void ratio
    double solid_volume = num_particles * (4.0 / 3) * CH_PI * radius * radius * radius;

    // Floor top in the bin frame
    double floor_top = thickness / 2;

    // Run metadata (rewritten at the end with the run summary)

    auto write_meta = [&](bool finished, double wall_time, unsigned long long overflows, const std::string& status) {
        JsonObject m;
        m.Add("program", "test_VDEM_directShear_squareBox");
        m.Add("description", "Direct shear in a square box, Chrono::Multicore SMC (Hertz), Y up, SI units");
        m.Add("chrono_version", CHRONO_VERSION);
        m.Add("start", start_date);
        m.Add("status", status);
        m.Add("sample", p.sample);
        m.Add("num_spheres", num_particles);
        m.Add("grid_layers_nx_nz", std::vector<double>{(double)p.layers, (double)p.nx, (double)p.nz});
        m.Add("radius", p.radius);
        m.Add("density", p.density);
        m.Add("sphere_mass", mass);
        m.Add("E", p.E);
        m.Add("nu", p.nu);
        m.Add("cor", p.cor);
        m.Add("mu_sphere", p.mu);
        m.Add("mu_wall", p.mu_wall);
        m.Add("mu_sphere_wall_composite", std::min(p.mu, p.mu_wall));
        m.Add("dense", p.dense);
        m.Add("mu_dense", p.mu_dense);
        m.Add("box_width", width);
        m.Add("box_length", length);
        m.Add("box_wall_height", height);
        m.Add("wall_thickness", thickness);
        m.Add("gap_between_boxes", radius);
        m.Add("nominal_area", area);
        m.Add("sigma", p.sigma);
        m.Add("lid_mass", plate->GetMass());
        m.Add("gravity", p.gravity);
        m.Add("shear_speed", p.speed);
        m.Add("shear_travel", p.travel);
        m.Add("t_settle", p.t_settle);
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
        m.Add("frame_dt_prep", p.frame_dt_prep);
        m.Add("save_state", p.save_state);
        m.Add("load_state", p.load_state);
        m.Add("first_particle_index", first_particle);
        m.AddRaw("walls",
                 "{\"bin\": {\"body\": 0, \"shapes\": {\"0\": \"floor\", \"1\": \"-x\", \"2\": \"+x\", \"3\": \"-z\", "
                 "\"4\": \"+z\"}}, \"box\": {\"body\": 1, \"shapes\": {\"0\": \"-x\", \"1\": \"+x\", \"2\": \"-z\", "
                 "\"3\": \"+z\"}}, \"lid\": {\"body\": 2, \"shapes\": {\"0\": \"lid\"}}}");
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
    ts << "# Direct shear, square box. SI units. Y up; the lower box moves in +z.\n";
    ts << "# phase: 0 fill, 1 compact, 2 shear. Forces are contact forces on the named body or walls.\n";
    ts << "# tau_bin = -Fz(bin)/A (2015 definition); tau_upper = Fz(upper box + lid)/A; ratio = tau/sigma.\n";
    ts << "# e_lid: void ratio from the lid height over the box area; e_cv: void ratio in a central control volume.\n";
    ts << "time,phase,shear_disp,lid_y,dilation,lid_tilt,"
          "bin_fx,bin_fy,bin_fz,floor_fx,floor_fy,floor_fz,binwalls_fx,binwalls_fy,binwalls_fz,"
          "boxwalls_fx,boxwalls_fy,boxwalls_fz,lid_fx,lid_fy,lid_fz,"
          "tau_bin,ratio_bin,tau_upper,ratio_upper,e_lid,e_cv,n_cv,kinetic_energy,"
          "contacts,contacts_pp,contacts_pw,coord,coord_mech,rattlers,frac_sliding,mean_fn,"
          "overflows_step,overflows_total,wall_clock,step_time,"
          "t_collision,t_broad,t_narrow,t_update,t_advance,t_contacts,t_ls_setup,t_ls_solve\n";
    ts << std::setprecision(9);

    std::ofstream frames_index(out_dir + "/frames.csv");
    frames_index << "frame,time,shear_disp,phase\n";

#ifdef CHRONO_VSG
    std::shared_ptr<vsg3d::ChVisualSystemVSG> vis;
    if (p.vsg) {
        vis = chrono_types::make_shared<vsg3d::ChVisualSystemVSG>();
        vis->AttachSystem(sys);
        vis->SetWindowTitle("Direct shear, square box");
        vis->SetWindowSize(1280, 720);
        vis->SetCameraVertical(CameraVerticalDir::Y);
        vis->AddCamera(ChVector3d(3 * width, 0.5 * height, 0.5 * width), ChVector3d(0, 0, 0));
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
    int prep_frame = 0;
    size_t next_snap = 0;
    double last_wall = 0;
    long long last_step = 0;
    long long step = 0;
    unsigned long long overflows_total = 0;

    auto shear_disp = [&]() { return shearing ? bin->GetPos().z() - shear_Disp : 0.0; };
    auto phase = [&]() { return settling ? 0 : (shearing ? 2 : 1); };

    auto write_snapshot = [&](const std::string& name) {
        auto contacts = GetContacts(*sys);
        std::vector<int> counts;
        ComputeContactStats(contacts, first_particle, num_particles, &counts);
        WriteParticles(particles_dir + "/" + name + ".csv", balls, radii, counts, sys->GetChTime(), shear_disp());
        WriteContacts(contacts_dir + "/" + name + ".csv", contacts, sys->GetChTime(), shear_disp());
    };

    while (sys->GetChTime() < end_time - 1e-9) {
        double time = sys->GetChTime();

        if (time > p.t_settle && settling) {
            if (p.dense)
                material->SetFriction((float)p.mu_dense);
            plate->SetPos(ChVector3d(0, height, 0));
            plate->SetFixed(false);
            settling = false;
        }

        if (time > p.t_shear && !shearing) {
            if (p.dense)
                material->SetFriction((float)p.mu);
            shear_Height = plate->GetPos().y();
            shear_Disp = bin->GetPos().z();
            shearing = true;
            if (!p.save_state.empty() && p.load_state.empty()) {
                if (SaveState(*sys, p.save_state, {shear_Height, shear_Disp}))
                    cout << "Saved state at t = " << time << " s to " << p.save_state << endl;
                else
                    std::cerr << "Cannot save state to " << p.save_state << endl;
            }
            sys->GetSettings()->solver.tangential_displ_mode = HistoryModel(p.history_shear);
        }

        if (shearing) {
            bin->SetPos(ChVector3d(0, -height / 2, -p.speed * p.t_shear + p.speed * time));
            bin->SetLinVel(ChVector3d(0, 0, p.speed));  // for tangential contact displacement history model
            bin->SetRot(QUNIT);
        } else {
            bin->SetPos(ChVector3d(0, -height / 2, 0));
            bin->SetLinVel(ChVector3d(0, 0, 0));
            bin->SetRot(QUNIT);
        }

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
        if (shearing) {
            double d = shear_disp();
            while (next_snap < p.snaps.size() && d >= p.snaps[next_snap] - 0.5 * p.speed * p.dt) {
                char name[64];
                snprintf(name, sizeof(name), "snap_d%05.2fmm", 1e3 * p.snaps[next_snap]);
                write_snapshot(name);
                next_snap++;
            }
            if (p.frame_disp > 0 && d >= frame * p.frame_disp - 0.5 * p.speed * p.dt) {
                char name[64];
                snprintf(name, sizeof(name), "frame_%05d", frame);
                write_snapshot(name);
                frames_index << frame << "," << time << "," << d << ",2" << std::endl;
                frame++;
            }
        } else if (p.frame_dt_prep > 0 && time >= prep_frame * p.frame_dt_prep - 0.5 * p.dt) {
            char name[64];
            snprintf(name, sizeof(name), "prep_%05d", prep_frame);
            write_snapshot(name);
            frames_index << "p" << prep_frame << "," << time << ",0," << phase() << std::endl;
            prep_frame++;
        }

        // Time series
        if (time >= data_out_frame * p.out_dt) {
            auto contacts = GetContacts(*sys);
            auto st = ComputeContactStats(contacts, first_particle, num_particles);

            // Contact forces on walls, by shape
            ChVector3d f_floor(0), f_binwalls(0), f_boxwalls(0), f_lid(0);
            for (const auto& c : contacts) {
                // Force on body 2 is F; on body 1 it is -F
                ChVector3d F = c.fn * c.n + c.ft;
                auto add = [&](int b, int l, const ChVector3d& f) {
                    if (b == 0)
                        (l == 0 ? f_floor : f_binwalls) += f;
                    else if (b == 1)
                        f_boxwalls += f;
                    else if (b == 2)
                        f_lid += f;
                };
                add(c.b1, c.l1, -F);
                add(c.b2, c.l2, F);
            }
            ChVector3d f_bin = f_floor + f_binwalls;

            double lid_y = plate->GetPos().y();
            double dilation = shearing ? lid_y - shear_Height : 0.0;
            double tilt = std::acos(std::min(1.0, plate->GetRotMat().GetAxisY().y()));
            double tau_bin = -f_bin.z() / area;
            double tau_upper = (f_boxwalls.z() + f_lid.z()) / area;

            // Void ratio from the lid height (lid bottom to floor top) over the nominal area
            double floor_y = bin->GetPos().y() + floor_top;
            double lid_bottom = lid_y - thickness / 2;
            double e_lid = settling ? NAN : (area * (lid_bottom - floor_y) - solid_volume) / solid_volume;

            // Void ratio in a central control volume: sphere centers in a box 2 diameters from the side walls
            // and from the floor and lid (an unbiased estimate of the solid fraction)
            double e_cv = NAN;
            int n_cv = 0;
            if (!settling) {
                double hx = width / 2 - 4 * radius;
                double y0 = floor_y + 4 * radius;
                double y1 = lid_bottom - 4 * radius;
                double z0 = bin->GetPos().z();  // the volume covers the overlap of the two halves in z
                if (y1 > y0) {
                    for (const auto& b : balls) {
                        const auto& q = b->GetPos();
                        if (std::abs(q.x()) < hx && std::abs(q.z() - 0.5 * z0) < hx - 0.5 * std::abs(z0) && q.y() > y0 && q.y() < y1)
                            n_cv++;
                    }
                    double v_cv = (2 * hx) * (2 * hx - std::abs(z0)) * (y1 - y0);
                    double v_s = n_cv * (4.0 / 3) * CH_PI * radius * radius * radius;
                    e_cv = v_s > 0 ? (v_cv - v_s) / v_s : NAN;
                }
            }

            double ke = 0;
            for (const auto& b : balls) {
                ke += 0.5 * b->GetMass() * b->GetPosDt().Length2();
                ke += 0.5 * 0.4 * b->GetMass() * radius * radius * b->GetAngVelParent().Length2();
            }

            double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall_start).count();
            double step_time = (step > last_step) ? (wall - last_wall) / (step - last_step) : 0;
            last_wall = wall;
            last_step = step;

            ts << time << "," << phase() << "," << shear_disp() << "," << lid_y << "," << dilation << "," << tilt << "," << f_bin.x() << "," << f_bin.y() << "," << f_bin.z() << ","
               << f_floor.x() << "," << f_floor.y() << "," << f_floor.z() << "," << f_binwalls.x() << "," << f_binwalls.y() << "," << f_binwalls.z() << "," << f_boxwalls.x() << ","
               << f_boxwalls.y() << "," << f_boxwalls.z() << "," << f_lid.x() << "," << f_lid.y() << "," << f_lid.z() << "," << tau_bin << "," << tau_bin / p.sigma << ","
               << tau_upper << "," << tau_upper / p.sigma << "," << e_lid << "," << e_cv << "," << n_cv << "," << ke << "," << st.num_contacts << "," << st.num_pp << ","
               << st.num_pw << "," << st.coord << "," << st.coord_mech << "," << st.num_rattlers << "," << st.frac_sliding << "," << st.mean_fn << ","
               << sys->GetNumContactHistoryOverflows() << "," << overflows_total << "," << wall << "," << step_time << "," << sys->GetTimerCollision() << ","
               << sys->GetTimerCollisionBroad() << "," << sys->GetTimerCollisionNarrow() << "," << sys->GetTimerUpdate() << "," << sys->GetTimerAdvance() << ","
               << sys->GetTimerProcessContact() << "," << sys->GetTimerLSsetup() << "," << sys->GetTimerLSsolve() << std::endl;

            if (data_out_frame % 10 == 0) {
                cout << std::fixed << std::setprecision(3) << "t = " << time << "  phase " << phase() << "  disp = " << std::setprecision(4) << 1e3 * shear_disp()
                     << " mm  ratio = " << tau_bin / p.sigma << "  e_lid = " << e_lid << "  contacts = " << st.num_contacts << "  step = " << std::setprecision(2)
                     << 1e3 * step_time << " ms  overflows = " << overflows_total << endl;
                cout.unsetf(std::ios::fixed);
            }
            data_out_frame++;
        }
    }

    double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall_start).count();
    write_meta(true, wall, overflows_total, "finished");
    cout << "Finished at t = " << sys->GetChTime() << " s, wall clock " << wall << " s, history overflows " << overflows_total << endl;

    return 0;
}

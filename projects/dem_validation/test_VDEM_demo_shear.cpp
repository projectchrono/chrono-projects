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
// Soft-sphere (SMC) or hard-sphere (NSC) direct shear box validation code.
// Problem parameters correspond to the Hartl and Ooi (2008) direct shear tests
// on glass beads.
//
// The global reference frame has Y up.
// All units SI.
// =============================================================================

#include <iostream>
#include <vector>
#include <valarray>
#include <string>
#include <sstream>
#include <cmath>

#include "chrono/ChConfig.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsGenerators.h"
#include "chrono/utils/ChUtilsInputOutput.h"
#include "chrono/assets/ChVisualShapeSphere.h"
#include "chrono/assets/ChVisualShapeBox.h"

#include "chrono_multicore/physics/ChSystemMulticore.h"
#include "chrono_multicore/solver/ChSystemDescriptorMulticore.h"

#include "chrono_thirdparty/filesystem/path.h"

#ifdef CHRONO_OPENGL
#include "chrono_opengl/ChVisualSystemOpenGL.h"
#endif

#include "../utils.h"

using namespace chrono;

using std::cout;
using std::flush;
using std::endl;

// -----------------------------------------------------------------------------
// Problem definitions
// -----------------------------------------------------------------------------

// Comment the following line to use NSC contact
#define USE_SMC

// Desired number of OpenMP threads (will be clamped to maximum available)
int threads = 20;

// Solver settings
#ifdef USE_SMC
double time_step = 1e-5;
double tolerance = 0.01;
int max_iteration_bilateral = 100;
#else
double time_step = 1e-4;
double tolerance = 0.1;
int max_iteration_normal = 0;
int max_iteration_sliding = 10000;
int max_iteration_spinning = 0;
int max_iteration_bilateral = 100;
double contact_recovery_speed = 10e30;
#endif

bool clamp_bilaterals = false;
double bilateral_clamp_speed = 0.1;

// Simulation parameters
#ifdef USE_SMC
double settling_time = 0.23;
double begin_shear_time = 2.0;
double end_simulation_time = 12.0;
double shear_speed = 0.001;  // m/s
#else
double settling_time = 0.23;
double begin_shear_time = 0.5;
double end_simulation_time = 2.5;
double shear_speed = 0.005;  // m/s
#endif

// Normal pressure (Pa)
// double normal_pressure = 24.2e3;
// double normal_pressure = 12.5e3;
// double normal_pressure = 6.4e3;
double normal_pressure = 3.1e3;

// Output
#ifdef USE_SMC
const std::string out_dir = "../SHEAR_SMC";
#else
const std::string out_dir = "../SHEAR_NSC";
#endif

const std::string pov_dir = out_dir + "/POVRAY";
const std::string shear_file = out_dir + "/shear_ratio.dat";
const std::string force_file = out_dir + "/shear_force.dat";
const std::string stats_file = out_dir + "/stats.dat";

bool write_povray_data = true;

double data_out_step = 1e-2;    // time interval between data outputs
double visual_out_step = 1e-1;  // time interval between PovRay outputs

// -----------------------------------------------------------------------------
// Utility for adding (visible or invisible) walls
// -----------------------------------------------------------------------------
void AddWall(std::shared_ptr<ChBody>& body,
             std::shared_ptr<ChContactMaterial> mat,
             const ChVector3d& dim,
             const ChVector3d& loc,
             bool visible) {
    auto ct_shape = chrono_types::make_shared<ChCollisionShapeBox>(mat, dim);
    body->AddCollisionShape(ct_shape, ChFrame<>(loc, QUNIT));

    if (visible == true) {
        auto box = chrono_types::make_shared<ChVisualShapeBox>(2.0 * dim);
        box->SetColor(ChColor(1, 0, 0));
        body->AddVisualShape(box, ChFrame<>(loc));
    }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    // Create output directories

    if (!filesystem::create_directory(filesystem::path(out_dir))) {
        cout << "Error creating directory " << out_dir << endl;
        return 1;
    }
    if (!filesystem::create_directory(filesystem::path(pov_dir))) {
        cout << "Error creating directory " << pov_dir << endl;
        return 1;
    }

    // Parameters for the system
    double gravity = 9.81;  // m/s^2

    // Parameters for the balls
    const int a = 50;
    const int b = 10;
    const int c = 10;
    int numballs = a * b * c;  // number of falling balls = (a X b X c)

    bool dense = false;

    double radius = 0.003;  // m
    double density = 2550;  // kg/m^3
    double mass = density * (4.0 / 3) * CH_PI * radius * radius * radius;
    float Y = 4.0e7;  // Pa
    float nu = 0.22f;
    float COR = 0.87f;
    float mu = 0.18f;

    // Parameters for containing bin, shear box, and load plate
    float mu_ext = 0.13f;

    int groundId = 0;
    double width = 0.12;
    double length = 0.12;
    double height = 0.06;
    double thickness = 0.01;

    double shear_Area;
    double shear_Height;
    double shear_Disp;

    ChVector3d pos(0, 0, 0);
    ChQuaternion<> rot(1, 0, 0, 0);
    ChVector3d vel(0, 0, 0);
    real3 force(0, 0, 0);

    // Define two quaternions representing:
    // - a rotation of -90 degrees around x (z2y)
    // - a rotation of +90 degrees around y (z2x)

    ChQuaternion<> z2y;
    ChQuaternion<> z2x;
    z2y.SetFromAngleAxis(-CH_PI / 2, ChVector3d(1, 0, 0));
    z2x.SetFromAngleAxis(CH_PI / 2, ChVector3d(0, 1, 0));

// Create the system

#ifdef USE_SMC
    cout << "Create SMC system" << endl;
    const std::string title = "soft-sphere (SMC) direct shear box test";
    ChSystemMulticoreSMC* sys = new ChSystemMulticoreSMC();
#else
    cout << "Create NSC system" << endl;
    const std::string title = "hard-sphere (NSC) direct shear box test";
    ChSystemMulticoreNSC* sys = new ChSystemMulticoreNSC();
#endif

    sys->SetCollisionSystemType(ChCollisionSystem::Type::MULTICORE);

    sys->SetGravitationalAcceleration(ChVector3d(0, -gravity, 0));

    // Set number of threads

    int max_threads = omp_get_num_procs();
    if (threads > max_threads)
        threads = max_threads;
    sys->SetNumThreads(threads);

    // Edit system settings

    sys->GetSettings()->solver.use_full_inertia_tensor = false;
    sys->GetSettings()->solver.tolerance = tolerance;
    sys->GetSettings()->solver.max_iteration_bilateral = max_iteration_bilateral;
    sys->GetSettings()->solver.clamp_bilaterals = clamp_bilaterals;
    sys->GetSettings()->solver.bilateral_clamp_speed = bilateral_clamp_speed;

#ifdef USE_SMC
    sys->GetSettings()->solver.contact_force_model = ChSystemSMC::ContactForceModel::Hertz;
    sys->GetSettings()->solver.tangential_displ_mode = ChSystemSMC::TangentialDisplacementModel::MultiStep;
#else
    sys->GetSettings()->solver.solver_mode = SolverMode::SLIDING;
    sys->GetSettings()->solver.max_iteration_normal = max_iteration_normal;
    sys->GetSettings()->solver.max_iteration_sliding = max_iteration_sliding;
    sys->GetSettings()->solver.max_iteration_spinning = max_iteration_spinning;
    sys->GetSettings()->solver.alpha = 0;
    sys->GetSettings()->solver.contact_recovery_speed = contact_recovery_speed;
    sys->ChangeSolverType(SolverType::APGD);

    sys->GetSettings()->collision.collision_envelope = 0.05 * radius;
#endif

    sys->GetSettings()->collision.bins_per_axis = vec3(10, 10, 10);
    sys->GetSettings()->collision.narrowphase_algorithm = ChNarrowphase::Algorithm::HYBRID;

// Create a ball material (will be used by balls only)

#ifdef USE_SMC
    auto material = chrono_types::make_shared<ChContactMaterialSMC>();
    material->SetYoungModulus(Y);
    material->SetPoissonRatio(nu);
    material->SetRestitution(COR);
    material->SetFriction(mu);
#else
    auto material = chrono_types::make_shared<ChContactMaterialNSC>();
    material->SetRestitution(COR);
    material->SetFriction(mu);
#endif

// Create a material for all objects other than balls

#ifdef USE_SMC
    auto mat_ext = chrono_types::make_shared<ChContactMaterialSMC>();
    mat_ext->SetYoungModulus(Y);
    mat_ext->SetPoissonRatio(nu);
    mat_ext->SetRestitution(COR);
    mat_ext->SetFriction(mu_ext);
#else
    auto mat_ext = chrono_types::make_shared<ChContactMaterialNSC>();
    mat_ext->SetRestitution(COR);
    mat_ext->SetFriction(mu_ext);
#endif

    // Create lower bin

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

    // Create upper shear box

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

    // Create upper load plate

    auto plate = chrono_types::make_shared<ChBody>();

    shear_Area = width * length;

    plate->SetMass(normal_pressure * shear_Area / gravity);
    plate->SetPos(ChVector3d(0, 2.0 * radius * float(a) + thickness, 0));
    plate->SetFixed(true);
    plate->EnableCollision(true);

    AddWall(plate, mat_ext, ChVector3d(width / 2, thickness / 2, length / 2), ChVector3d(0, 0, 0), true);
    plate->GetCollisionModel()->SetFamily(3);
    plate->GetCollisionModel()->DisallowCollisionsWith(1);
    plate->GetCollisionModel()->DisallowCollisionsWith(2);
    plate->GetCollisionModel()->AllowCollisionsWith(4);

    sys->AddBody(plate);

    // Create (a X b X c) many falling balls

    int i, j, k;
    double ball_x, ball_y, ball_z;

    for (i = 0; i < a; i++) {
        for (j = 0; j < b; j++) {
            for (k = 0; k < c; k++) {
                ball_y = 2.0 * radius * float(i);

                ball_x = 4.0 * radius * (float(j - b / 2) + 0.5) + 0.99 * radius * (float(rand() % 100) / 50 - 1.0);
                ball_z = 4.0 * radius * (float(k - c / 2) + 0.5) + 0.99 * radius * (float(rand() % 100) / 50 - 1.0);

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
                sphere->SetColor(ChColor(1, 0, 1));
                ball->AddVisualShape(sphere);

                sys->AddBody(ball);
            }
        }
    }

    // Create prismatic (translational) joint between load plate and shear box.
    // The translational axis of a prismatic joint is along the Z axis of the
    // specified joint coordinate system.  Here, we apply the 'z2y' rotation to
    // align it with the Y axis of the global reference frame.

    auto prismatic_plate_box = chrono_types::make_shared<ChLinkLockPrismatic>();
    prismatic_plate_box->SetName("prismatic_plate_box");
    prismatic_plate_box->Initialize(plate, box, ChFrame<>(ChVector3d(0, 0, 0), z2y));
    sys->AddLink(prismatic_plate_box);

    // Setup output

    std::ofstream shearStream(shear_file.c_str());
    std::ofstream forceStream(force_file.c_str());
    std::ofstream statsStream(stats_file.c_str());

// Create the OpenGL visualization window

#ifdef CHRONO_OPENGL
    opengl::ChVisualSystemOpenGL vis;
    vis.AttachSystem(sys);
    vis.SetWindowTitle(title);
    vis.SetWindowSize(1280, 720);
    vis.SetRenderMode(opengl::WIREFRAME);
    vis.Initialize();
    vis.AddCamera(ChVector3d(3 * width, 0, 0), ChVector3d(0, 0, 0));
    vis.SetCameraVertical(CameraVerticalDir::Y);
#endif

    // Begin simulation

    bool settling = true;
    bool shearing = false;

    int data_out_frame = 0;
    int visual_out_frame = 0;

    while (sys->GetChTime() < end_simulation_time) {
        if (sys->GetChTime() > settling_time && settling == true) {
            if (dense == true)
                material->SetFriction(0.01f);
            plate->SetPos(ChVector3d(0, height, 0));
            plate->SetFixed(false);
            settling = false;
        }

        if (sys->GetChTime() > begin_shear_time && shearing == false) {
            if (dense == true)
                material->SetFriction(mu);
            shear_Height = plate->GetPos().y();
            shear_Disp = bin->GetPos().z();
            shearing = true;
        }

        if (shearing == true) {
            bin->SetPos(ChVector3d(0, -height / 2, -shear_speed * begin_shear_time + shear_speed * sys->GetChTime()));
            bin->SetLinVel(ChVector3d(0, 0, shear_speed));  // for tangential contact displacement history model
            bin->SetRot(QUNIT);
        } else {
            bin->SetPos(ChVector3d(0, -height / 2, 0));
            bin->SetLinVel(ChVector3d(0, 0, 0));
            bin->SetRot(QUNIT);
        }

//  Do time step

#ifdef CHRONO_OPENGL
        if (vis.Run()) {
            sys->DoStepDynamics(time_step);
            vis.Render();
        } else
            break;
#else
        sys->DoStepDynamics(time_step);
#endif

        TimingOutput(sys, &statsStream);

        //  Output to files

        if (sys->GetChTime() >= data_out_frame * data_out_step) {
#ifndef USE_SMC
            sys->CalculateContactForces();
#endif
            force = sys->GetBodyContactForce(0);

            forceStream << sys->GetChTime() << "\t" << plate->GetPos().y() - bin->GetPos().y() << "\t"
                        << bin->GetPos().x() << "\t" << bin->GetPos().y() << "\t" << bin->GetPos().z() << "\t" << force.x
                        << "\t" << force.y << "\t" << force.z << "\n";

            cout << sys->GetChTime() << "\t" << plate->GetPos().y() - bin->GetPos().y() << "\t" << bin->GetPos().x()
                 << "\t" << bin->GetPos().y() << "\t" << bin->GetPos().z() << "\t" << force.x << "\t" << force.y << "\t"
                 << force.z << "\n";

            //  Output to shear data file

            if (shearing == true) {
                shearStream << (bin->GetPos().z() - shear_Disp) / (2 * radius) << "\t";
                shearStream << -force.z / (shear_Area * normal_pressure) << "\t";
                shearStream << (plate->GetPos().y() - shear_Height) / (2 * radius) << "\n";

                cout << (bin->GetPos().z() - shear_Disp) / (2 * radius) << "\t";
                cout << -force.z / (shear_Area * normal_pressure) << "\t";
                cout << (plate->GetPos().y() - shear_Height) / (2 * radius) << "\n";
            }

            data_out_frame++;
        }

        //  Output to POV-Ray

        if (write_povray_data && sys->GetChTime() >= visual_out_frame * visual_out_step) {
            char filename[100];
            sprintf(filename, "%s/data_%03d.dat", pov_dir.c_str(), visual_out_frame + 1);
            utils::WriteVisualizationAssets(sys, filename, false);

            visual_out_frame++;
        }
    }

    return 0;
}

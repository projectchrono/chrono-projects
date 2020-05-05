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
// Author: Radu Serban
// =============================================================================
//
// ChronoParallel demo program for low velocity cratering studies.
//
// The model simulated here consists of a spherical projectile dropped in a
// bed of granular material, using either penalty or complementarity method for
// frictional contact.
//
// The global reference frame has Z up.
// All units SI.
// =============================================================================

#include <cstdio>
#include <vector>
#include <cmath>

#include "chrono/ChConfig.h"
#include "chrono/core/ChStream.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsGenerators.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_parallel/solver/ChSystemDescriptorParallel.h"

#include "chrono_thirdparty/filesystem/path.h"

// Control use of OpenGL run-time rendering
// Note: CHRONO_OPENGL is defined in ChConfig.h
#undef CHRONO_OPENGL

#ifdef CHRONO_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif

#include "../utils.h"

using namespace chrono;
using namespace chrono::collision;

using std::cout;
using std::flush;
using std::endl;

// -----------------------------------------------------------------------------
// Problem definitions
// -----------------------------------------------------------------------------

// Comment the following line to use NSC contact
#define USE_SMC

enum ProblemType { SETTLING, DROPPING };

ProblemType problem = SETTLING;

// -----------------------------------------------------------------------------
// Global problem definitions
// -----------------------------------------------------------------------------

// Desired number of OpenMP threads (will be clamped to maximum available)
int threads = 20;

// Perform dynamic tuning of number of threads?
bool thread_tuning = false;

// Simulation parameters
double gravity = 9.81;

double time_settling_min = 0.1;
double time_settling_max = 0.8;
double time_dropping = 0.06;

#ifdef USE_SMC
double time_step = 1e-5;
int max_iteration = 20;
#else
double time_step = 1e-4;
int max_iteration_normal = 0;
int max_iteration_sliding = 50;
int max_iteration_spinning = 0;
float contact_recovery_speed = 0.1;
#endif

double tolerance = 1.0;

// Contact force model
#ifdef USE_SMC
ChSystemSMC::ContactForceModel contact_force_model = ChSystemSMC::ContactForceModel::Hooke;
ChSystemSMC::TangentialDisplacementModel tangential_displ_mode = ChSystemSMC::TangentialDisplacementModel::MultiStep;
#endif

// Output
bool povray_output = true;

#ifdef USE_SMC
const std::string out_dir = "../CRATER_SMC";
#else
const std::string out_dir = "../CRATER_NSC";
#endif

const std::string pov_dir = out_dir + "/POVRAY";
const std::string height_file = out_dir + "/height.dat";
const std::string stats_file = out_dir + "/stats.dat";
const std::string checkpoint_file = out_dir + "/settled.dat";

int out_fps_settling = 120;
int out_fps_dropping = 1200;

int timing_frame = -1;  // output detailed step timing at this frame

// Parameters for the granular material
int Id_g = 1;
double r_g = 1e-3 / 2;
double rho_g = 2500;
double vol_g = (4.0 / 3) * CH_C_PI * r_g * r_g * r_g;
double mass_g = rho_g * vol_g;
ChVector<> inertia_g = 0.4 * mass_g * r_g * r_g * ChVector<>(1, 1, 1);

float Y_g = 1e8f;
float mu_g = 0.3f;
float cr_g = 0.1f;

// Parameters for the falling ball
int Id_b = 0;
double R_b = 2.54e-2 / 2;
double rho_b = 700;
double vol_b = (4.0 / 3) * CH_C_PI * R_b * R_b * R_b;
double mass_b = rho_b * vol_b;
ChVector<> inertia_b = 0.4 * mass_b * R_b * R_b * ChVector<>(1, 1, 1);

float Y_b = 1e8f;
float mu_b = 0.3f;
float cr_b = 0.1f;

// Parameters for the containing bin
int binId = -200;
double hDimX = 4e-2;         // length in x direction
double hDimY = 4e-2;         // depth in y direction
double hDimZ = 7.5e-2;       // height in z direction
double hThickness = 0.5e-2;  // wall thickness

float Y_c = 2e6f;
float mu_c = 0.3f;
float cr_c = 0.1f;

// Number of layers and height of one layer for generator domain
int numLayers = 10;
double layerHeight = 1e-2;
////int numLayers = 1;
////double layerHeight = 3e-3;

// Drop height (above surface of settled granular material)
double h = 10e-2;

// -----------------------------------------------------------------------------
// Create the dynamic objects:
// - granular material consisting of identical spheres with specified radius and
//   material properties; the spheres are generated in a number of vertical
//   layers with locations within each layer obtained using Poisson Disk
//   sampling (thus ensuring that no two spheres are closer than twice the
//   radius)
// - a containing bin consisting of five boxes (no top)
// -----------------------------------------------------------------------------
int CreateObjects(ChSystemParallel* system) {
// Create the containing bin
#ifdef USE_SMC
    auto mat_c = chrono_types::make_shared<ChMaterialSurfaceSMC>();
    mat_c->SetYoungModulus(Y_c);
    mat_c->SetFriction(mu_c);
    mat_c->SetRestitution(cr_c);

    utils::CreateBoxContainer(system, binId, mat_c, ChVector<>(hDimX, hDimY, hDimZ), hThickness);
#else
    auto mat_c = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    mat_c->SetFriction(mu_c);

    utils::CreateBoxContainer(system, binId, mat_c, ChVector<>(hDimX, hDimY, hDimZ), hThickness);
#endif

// Create a material for the granular material
#ifdef USE_SMC
    auto mat_g = chrono_types::make_shared<ChMaterialSurfaceSMC>();
    mat_g->SetYoungModulus(Y_g);
    mat_g->SetFriction(mu_g);
    mat_g->SetRestitution(cr_g);
#else
    auto mat_g = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    mat_g->SetFriction(mu_g);
#endif

    // Create a mixture entirely made out of spheres
    utils::Generator gen(system);

    std::shared_ptr<utils::MixtureIngredient> m1 = gen.AddMixtureIngredient(utils::MixtureType::SPHERE, 1.0);
    m1->setDefaultMaterial(mat_g);
    m1->setDefaultDensity(rho_g);
    m1->setDefaultSize(r_g);

    gen.setBodyIdentifier(Id_g);

    double r = 1.01 * r_g;

    for (int i = 0; i < numLayers; i++) {
        double center = r + layerHeight / 2 + i * (2 * r + layerHeight);
        gen.createObjectsBox(utils::SamplingType::POISSON_DISK, 2 * r, ChVector<>(0, 0, center),
                             ChVector<>(hDimX - r, hDimY - r, layerHeight / 2));
        cout << "Layer " << i << "  total bodies: " << gen.getTotalNumBodies() << endl;
    }

    return gen.getTotalNumBodies();
}

// -----------------------------------------------------------------------------
// Create the falling ball such that its bottom point is at the specified height
// and its downward initial velocity has the specified magnitude.
// -----------------------------------------------------------------------------
void CreateFallingBall(ChSystemParallel* system, double z, double vz) {
    // Create a material for the falling ball
#ifdef USE_SMC
    auto mat_b = chrono_types::make_shared<ChMaterialSurfaceSMC>();
    mat_b->SetYoungModulus(1e8f);
    mat_b->SetFriction(0.4f);
    mat_b->SetRestitution(0.1f);
#else
    auto mat_b = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    mat_b->SetFriction(mu_c);
#endif

    // Create the falling ball
    auto ball = chrono_types::make_shared<ChBody>(chrono_types::make_shared<ChCollisionModelParallel>());

    ball->SetIdentifier(Id_b);
    ball->SetMass(mass_b);
    ball->SetInertiaXX(inertia_b);
    ball->SetPos(ChVector<>(0, 0, z + r_g + R_b));
    ball->SetRot(ChQuaternion<>(1, 0, 0, 0));
    ball->SetPos_dt(ChVector<>(0, 0, -vz));
    ball->SetCollide(true);
    ball->SetBodyFixed(false);

    ball->GetCollisionModel()->ClearModel();
    utils::AddSphereGeometry(ball.get(), mat_b, R_b);
    ball->GetCollisionModel()->BuildModel();

    system->AddBody(ball);
}

// -----------------------------------------------------------------------------
// Find the height of the highest and lowest, respectively, sphere in the
// granular mix, respectively.  We only look at bodies whith stricty positive
// identifiers (to exclude the containing bin).
// -----------------------------------------------------------------------------
double FindHighest(ChSystem* sys) {
    double highest = 0;
    for (auto body : sys->Get_bodylist()) {
        if (body->GetIdentifier() > 0 && body->GetPos().z() > highest)
            highest = body->GetPos().z();
    }
    return highest;
}

double FindLowest(ChSystem* sys) {
    double lowest = 1000;
    for (auto body : sys->Get_bodylist()) {
        if (body->GetIdentifier() > 0 && body->GetPos().z() < lowest)
            lowest = body->GetPos().z();
    }
    return lowest;
}

// -----------------------------------------------------------------------------
// Return true if all bodies in the granular mix have a linear velocity whose
// magnitude is below the specified value.
// -----------------------------------------------------------------------------
bool CheckSettled(ChSystem* sys, double threshold) {
    double t2 = threshold * threshold;

    for (auto body : sys->Get_bodylist()) {
        if (body->GetIdentifier() > 0) {
            double vel2 = body->GetPos_dt().Length2();
            if (vel2 > t2)
                return false;
        }
    }

    return true;
}

// -----------------------------------------------------------------------------
int main(int argc, char* argv[]) {
// Create system
#ifdef USE_SMC
    cout << "Create SMC system" << endl;
    ChSystemParallelSMC* msystem = new ChSystemParallelSMC();
#else
    cout << "Create NSC system" << endl;
    ChSystemParallelNSC* msystem = new ChSystemParallelNSC();
#endif

    // Debug log messages.
    ////msystem->SetLoggingLevel(LOG_INFO, true);
    ////msystem->SetLoggingLevel(LOG_TRACE, true);

    // Set number of threads.
    int max_threads = omp_get_num_procs();
    if (threads > max_threads)
        threads = max_threads;
    omp_set_num_threads(threads);
    cout << "Using " << threads << " threads" << endl;

    msystem->GetSettings()->perform_thread_tuning = thread_tuning;

    // Set gravitational acceleration
    msystem->Set_G_acc(ChVector<>(0, 0, -gravity));

    // Edit system settings
    msystem->GetSettings()->solver.use_full_inertia_tensor = false;
    msystem->GetSettings()->solver.tolerance = tolerance;

#ifdef USE_SMC
    msystem->GetSettings()->collision.narrowphase_algorithm = NarrowPhaseType::NARROWPHASE_R;

    msystem->GetSettings()->solver.contact_force_model = contact_force_model;
    msystem->GetSettings()->solver.tangential_displ_mode = tangential_displ_mode;
#else
    msystem->GetSettings()->solver.solver_mode = SolverMode::SLIDING;
    msystem->GetSettings()->solver.max_iteration_normal = max_iteration_normal;
    msystem->GetSettings()->solver.max_iteration_sliding = max_iteration_sliding;
    msystem->GetSettings()->solver.max_iteration_spinning = max_iteration_spinning;
    msystem->GetSettings()->solver.alpha = 0;
    msystem->GetSettings()->solver.contact_recovery_speed = contact_recovery_speed;
    msystem->ChangeSolverType(SolverType::APGDREF);

    msystem->GetSettings()->collision.collision_envelope = 0.05 * r_g;
#endif

    msystem->GetSettings()->collision.bins_per_axis = vec3(20, 20, 20);

    // Depending on problem type:
    // - Select end simulation time
    // - Select output FPS
    // - Create granular material and container
    // - Create falling ball
    double time_end;
    int out_fps;
    std::shared_ptr<ChBody> ball;

    if (problem == SETTLING) {
        time_end = time_settling_max;
        out_fps = out_fps_settling;

        cout << "Create granular material" << endl;
        // Create the fixed falling ball just below the granular material
        CreateFallingBall(msystem, -3 * R_b, 0);
        ball = msystem->Get_bodylist().at(0);
        ball->SetBodyFixed(true);
        CreateObjects(msystem);
    } else {
        time_end = time_dropping;
        out_fps = out_fps_dropping;

        // Create the granular material and the container from the checkpoint file.
        cout << "Read checkpoint data from " << checkpoint_file;
        utils::ReadCheckpoint(msystem, checkpoint_file);
        cout << "  done.  Read " << msystem->Get_bodylist().size() << " bodies." << endl;

        // Move the falling ball just above the granular material with a velocity
        // given by free fall from the specified height and starting at rest.
        double z = FindHighest(msystem);
        double vz = std::sqrt(2 * gravity * h);
        cout << "Move falling ball with center at " << z + R_b + r_g << " and velocity " << vz << endl;
        ball = msystem->Get_bodylist().at(0);
        ball->SetMass(mass_b);
        ball->SetInertiaXX(inertia_b);
        ball->SetPos(ChVector<>(0, 0, z + r_g + R_b));
        ball->SetRot(ChQuaternion<>(1, 0, 0, 0));
        ball->SetPos_dt(ChVector<>(0, 0, -vz));
        ball->SetBodyFixed(false);
    }

    // Number of steps
    int num_steps = (int)std::ceil(time_end / time_step);
    int out_steps = (int)std::ceil((1.0 / time_step) / out_fps);

    // Zero velocity level for settling check
    // (fraction of a grain radius per second)
    double zero_v = 0.1 * r_g;

    // Create output directories.
    if (!filesystem::create_directory(filesystem::path(out_dir))) {
        cout << "Error creating directory " << out_dir << endl;
        return 1;
    }
    if (!filesystem::create_directory(filesystem::path(pov_dir))) {
        cout << "Error creating directory " << pov_dir << endl;
        return 1;
    }

    // Perform the simulation
    double time = 0;
    int sim_frame = 0;
    int out_frame = 0;
    int next_out_frame = 0;
    double exec_time = 0;
    int num_contacts = 0;
    ChStreamOutAsciiFile sfile(stats_file.c_str());
    ChStreamOutAsciiFile hfile(height_file.c_str());

#ifdef CHRONO_OPENGL
    opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
    gl_window.Initialize(1280, 720, "Crater Test", msystem);
    gl_window.SetCamera(ChVector<>(0, -10 * hDimY, hDimZ), ChVector<>(0, 0, hDimZ), ChVector<>(0, 0, 1));
    gl_window.SetRenderMode(opengl::WIREFRAME);
#endif

    while (time < time_end) {
        if (sim_frame == next_out_frame) {
            cout << endl;
            cout << "---- Frame:          " << out_frame << endl;
            cout << "     Sim frame:      " << sim_frame << endl;
            cout << "     Time:           " << time << endl;
            cout << "     Lowest point:   " << FindLowest(msystem) << endl;
            cout << "     Avg. contacts:  " << num_contacts / out_steps << endl;
            cout << "     Execution time: " << exec_time << endl;

            sfile << time << "  " << exec_time << "  " << num_contacts / out_steps << "\n";

            // If enabled, output data for PovRay postprocessing.
            if (povray_output) {
                char filename[100];
                sprintf(filename, "%s/data_%03d.dat", pov_dir.c_str(), out_frame + 1);
                utils::WriteShapesPovray(msystem, filename, false);
            }

            // Create a checkpoint from the current state.
            if (problem == SETTLING) {
                cout << "     Write checkpoint data " << flush;
                utils::WriteCheckpoint(msystem, checkpoint_file);
                cout << msystem->Get_bodylist().size() << " bodies" << endl;
            }

            // Save current projectile height.
            if (problem == DROPPING) {
                hfile << time << "  " << ball->GetPos().z() << "\n";
                cout << "     Ball height:    " << ball->GetPos().z() << endl;
            }

            out_frame++;
            next_out_frame += out_steps;
            num_contacts = 0;
        }

        if (problem == SETTLING && time > time_settling_min && CheckSettled(msystem, zero_v)) {
            cout << "Granular material settled...  time = " << time << endl;
            break;
        }

// Advance simulation by one step
#ifdef CHRONO_OPENGL
        if (gl_window.Active()) {
            gl_window.DoStepDynamics(time_step);
            gl_window.Render();
        } else
            break;
#else
        msystem->DoStepDynamics(time_step);
#endif

        ////progressbar(out_steps + sim_frame - next_out_frame + 1, out_steps);
        TimingOutput(msystem);

        time += time_step;
        sim_frame++;
        exec_time += msystem->GetTimerStep();
        num_contacts += msystem->GetNcontacts();

        // If requested, output detailed timing information for this step
        if (sim_frame == timing_frame)
            msystem->PrintStepStats();
    }

    // Create a checkpoint from the last state
    if (problem == SETTLING) {
        cout << "Write checkpoint data to " << checkpoint_file;
        utils::WriteCheckpoint(msystem, checkpoint_file);
        cout << "  done.  Wrote " << msystem->Get_bodylist().size() << " bodies." << endl;
    }

    // Final stats
    cout << "==================================" << endl;
    cout << "Number of bodies:  " << msystem->Get_bodylist().size() << endl;
    cout << "Lowest position:   " << FindLowest(msystem) << endl;
    cout << "Simulation time:   " << exec_time << endl;
    cout << "Number of threads: " << threads << endl;

    return 0;
}

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
// ChronoParallel demo program for cohesive SMC granular material simulation.
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
//#undef CHRONO_OPENGL

#ifdef CHRONO_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif

#include "../utils.h"

using namespace chrono;
using namespace chrono::collision;

using std::cout;
using std::endl;

// =======================================================================
// Global problem definitions

// Desired number of OpenMP threads (will be clamped to maximum available)
int threads = 8;

// Perform dynamic tuning of number of threads?
bool thread_tuning = false;

// Simulation parameters
double gravity = 9.81;
double time_step = 1e-4;
double time_end = 5;

int max_iteration = 20;

// Output
const std::string out_dir = "../FOAM";
const std::string pov_dir = out_dir + "/POVRAY";
const std::string out_file = out_dir + "/timing.dat";
double out_fps = 50;

// Parameters for the granular material
int Id_g = 1;
double r_g = 0.01;
double rho_g = 2000;

float Y_g = 2e7f;
float mu_g = 0.2f;
float cr_g = 0.1f;
float cohesion_g = 300.0f;

// Parameters for the containing bin
int binId = -200;
double hDimX = 10;        // length in x direction
double hDimY = 10;        // depth in y direction
double hDimZ = 1;         // height in z direction
double hThickness = 0.4;  // wall thickness

float Y_c = 2e6f;
float mu_c = 0.3f;
float cr_c = 0.1f;
float cohesion_c = 5.0f;

// Particle generator
utils::Generator* gen;

double initVel = 5;  // initial particle velocity in negative X direction

int maxNumParticles = 100000;

// =======================================================================

int SpawnParticles() {
    double dist = 2 * 0.99 * r_g;

    ////gen->createObjectsBox(utils::SamplingType::POISSON_DISK,
    ////                     dist,
    ////                     ChVector<>(9, 0, 3),
    ////                     ChVector<>(0, 1, 0.5),
    ////                     ChVector<>(-initVel, 0, 0));
    gen->createObjectsCylinderX(utils::SamplingType::POISSON_DISK, dist, ChVector<>(9, 0, 3), 0.2f, 0, ChVector<>(-initVel, 0, 0));
    cout << "  total bodies: " << gen->getTotalNumBodies() << endl;

    return gen->getTotalNumBodies();
}

// ========================================================================
int main(int argc, char* argv[]) {
    // Create system
    ChSystemParallelSMC* msystem = new ChSystemParallelSMC();

    // Set number of threads.
    int max_threads = omp_get_num_procs();
    if (threads > max_threads)
        threads = max_threads;
    omp_set_num_threads(threads);

    msystem->GetSettings()->perform_thread_tuning = thread_tuning;

    // Set gravitational acceleration
    msystem->Set_G_acc(ChVector<>(0, 0, -gravity));

    // Using constant adhesion model
    msystem->GetSettings()->solver.adhesion_force_model = ChSystemSMC::AdhesionForceModel::Constant;

    // Edit system settings
    msystem->GetSettings()->solver.tolerance = 1e-4;

    msystem->GetSettings()->collision.bins_per_axis = vec3(10, 10, 10);

    msystem->GetSettings()->collision.narrowphase_algorithm = NarrowPhaseType::NARROWPHASE_R;

    // Create a material for the granular material
    auto mat_g = chrono_types::make_shared<ChMaterialSurfaceSMC>();
    mat_g->SetYoungModulus(Y_g);
    mat_g->SetFriction(mu_g);
    mat_g->SetRestitution(cr_g);
    mat_g->SetAdhesion(cohesion_g);

    // Create a material for the container
    auto mat_c = chrono_types::make_shared<ChMaterialSurfaceSMC>();
    mat_c->SetYoungModulus(Y_c);
    mat_c->SetFriction(mu_c);
    mat_c->SetRestitution(cr_c);
    mat_c->SetAdhesion(cohesion_c);

    // Create the containing bin
    utils::CreateBoxContainer(msystem, binId, mat_c, ChVector<>(hDimX, hDimY, hDimZ), hThickness);

    // Create a mixture entirely made out of spheres
    double vol_g = (4.0 / 3) * CH_C_PI * r_g * r_g * r_g;
    double mass_g = rho_g * vol_g;
    ChVector<> inertia_g = 0.4 * mass_g * r_g * r_g * ChVector<>(1, 1, 1);

    gen = new utils::Generator(msystem);

    std::shared_ptr<utils::MixtureIngredient> m1 = gen->AddMixtureIngredient(utils::MixtureType::SPHERE, 1.0);
    m1->setDefaultMaterial(mat_g);
    m1->setDefaultDensity(rho_g);
    m1->setDefaultSize(r_g);

    gen->setBodyIdentifier(Id_g);

    // Number of steps
    int num_steps = (int)std::ceil(time_end / time_step);
    int out_steps = (int)std::ceil((1 / time_step) / out_fps);
    int gen_steps = (int)std::ceil(3 * r_g / initVel / time_step);

    // Create output directories.
    if (!filesystem::create_directory(filesystem::path(out_dir))) {
        cout << "Error creating directory " << out_dir << endl;
        return 1;
    }
    if (!filesystem::create_directory(filesystem::path(pov_dir))) {
        cout << "Error creating directory " << pov_dir << endl;
        return 1;
    }

#ifdef CHRONO_OPENGL
    // Initialize OpenGL
    opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
    gl_window.Initialize(1280, 720, "Foam demo", msystem);
    gl_window.SetCamera(ChVector<>(4, -5, 4), ChVector<>(9, 0, 3), ChVector<>(0, 0, 1));
    gl_window.SetRenderMode(opengl::WIREFRAME);
#endif

    // Perform the simulation
    double time = 0;
    int sim_frame = 0;
    int out_frame = 0;
    double exec_time = 0;
    ChStreamOutAsciiFile ofile(out_file.c_str());

    while (time < time_end) {
        int numParticles = (int)msystem->Get_bodylist().size() - 1;

        if (numParticles < maxNumParticles && sim_frame % gen_steps == 0) {
            SpawnParticles();
        }

        if (sim_frame % out_steps == 0) {
            char filename[100];
            sprintf(filename, "%s/data_%03d.dat", pov_dir.c_str(), out_frame + 1);
            utils::WriteShapesPovray(msystem, filename);

            cout << " --------------------------------- Output frame:   " << out_frame + 1 << endl;
            cout << "                                   Sim frame:      " << sim_frame << endl;
            cout << "                                   Time:           " << time << endl;
            cout << "                                   Execution time: " << exec_time << endl;
            cout << "                                   Num. bodies:    " << numParticles << endl;

            ofile << sim_frame << "  " << time << "  " << exec_time << "  " << numParticles << "  "
                  << msystem->GetNcontacts() << "\n";

            out_frame++;
        }

// Advance dynamics.
#ifdef CHRONO_OPENGL
        if (gl_window.Active()) {
            gl_window.DoStepDynamics(time_step);
            gl_window.Render();
        } else
            break;
#else
        msystem->DoStepDynamics(time_step);
#endif

        time += time_step;
        sim_frame++;
        exec_time += msystem->GetTimerStep();
    }

    // Final stats
    cout << "==================================" << endl;
    cout << "Number of bodies: " << msystem->Get_bodylist().size() << endl;
    cout << "Simulation time: " << exec_time << endl;
    cout << "Number of threads: " << threads << endl;

    return 0;
}

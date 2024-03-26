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
// Chrono::Multicore demo program for cohesive SMC granular material simulation.
//
// The global reference frame has Z up.
// All units SI.
// =============================================================================

#include <cstdio>
#include <vector>
#include <cmath>

#include "chrono/ChConfig.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsGenerators.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#include "chrono_multicore/physics/ChSystemMulticore.h"
#include "chrono_multicore/solver/ChSystemDescriptorMulticore.h"

#include "chrono_thirdparty/filesystem/path.h"

// Control use of OpenGL run-time rendering
// Note: CHRONO_OPENGL is defined in ChConfig.h
//#undef CHRONO_OPENGL

#ifdef CHRONO_OPENGL
#include "chrono_opengl/ChVisualSystemOpenGL.h"
#endif

#include "../utils.h"

using namespace chrono;

using std::cout;
using std::endl;

// =======================================================================
// Global problem definitions

// Desired number of OpenMP threads (will be clamped to maximum available)
int threads = 8;

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
int tag_particles = 0;
double r_g = 0.01;
double rho_g = 2000;

float Y_g = 2e7f;
float mu_g = 0.2f;
float cr_g = 0.1f;
float cohesion_g = 300.0f;

// Parameters for the containing bin
double hDimX = 10;        // length in x direction
double hDimY = 10;        // depth in y direction
double hDimZ = 1;         // height in z direction
double hThickness = 0.4;  // wall thickness

float Y_c = 2e6f;
float mu_c = 0.3f;
float cr_c = 0.1f;
float cohesion_c = 5.0f;

// Particle generator
utils::ChGenerator* gen;

double initVel = 5;  // initial particle velocity in negative X direction

int maxNumParticles = 100000;

// =======================================================================

int SpawnParticles() {
    double dist = 2 * 0.99 * r_g;
    utils::ChPDSampler<double> sampler(dist);

    ////gen->CreateObjectsBox(sampler,
    ////                      ChVector3d(9, 0, 3),
    ////                      ChVector3d(0, 1, 0.5),
    ////                      ChVector3d(-initVel, 0, 0));
    gen->CreateObjectsCylinderX(sampler, ChVector3d(9, 0, 3), 0.2f, 0, ChVector3d(-initVel, 0, 0));
    cout << "  total bodies: " << gen->GetTotalNumBodies() << endl;

    return gen->GetTotalNumBodies();
}

// ========================================================================
int main(int argc, char* argv[]) {
    // Create system
    ChSystemMulticoreSMC* sys = new ChSystemMulticoreSMC();

    // Set number of threads.
    int max_threads = omp_get_num_procs();
    if (threads > max_threads)
        threads = max_threads;
    sys->SetNumThreads(threads);

    // Set gravitational acceleration
    sys->SetGravitationalAcceleration(ChVector3d(0, 0, -gravity));

    // Using constant adhesion model
    sys->GetSettings()->solver.adhesion_force_model = ChSystemSMC::AdhesionForceModel::Constant;

    // Edit system settings
    sys->GetSettings()->solver.tolerance = 1e-4;

    sys->GetSettings()->collision.bins_per_axis = vec3(10, 10, 10);

    sys->GetSettings()->collision.narrowphase_algorithm = ChNarrowphase::Algorithm::PRIMS;

    // Create a material for the granular material
    auto mat_g = chrono_types::make_shared<ChContactMaterialSMC>();
    mat_g->SetYoungModulus(Y_g);
    mat_g->SetFriction(mu_g);
    mat_g->SetRestitution(cr_g);
    mat_g->SetAdhesion(cohesion_g);

    // Create a material for the container
    auto mat_c = chrono_types::make_shared<ChContactMaterialSMC>();
    mat_c->SetYoungModulus(Y_c);
    mat_c->SetFriction(mu_c);
    mat_c->SetRestitution(cr_c);
    mat_c->SetAdhesion(cohesion_c);

    // Create the containing bin
    utils::CreateBoxContainer(sys, mat_c, ChVector3d(hDimX, hDimY, hDimZ), hThickness);

    // Create a mixture entirely made out of spheres
    double vol_g = (4.0 / 3) * CH_PI * r_g * r_g * r_g;
    double mass_g = rho_g * vol_g;
    ChVector3d inertia_g = 0.4 * mass_g * r_g * r_g * ChVector3d(1, 1, 1);

    gen = new utils::ChGenerator(sys);

    std::shared_ptr<utils::ChMixtureIngredient> m1 = gen->AddMixtureIngredient(utils::MixtureType::SPHERE, 1.0);
    m1->SetDefaultMaterial(mat_g);
    m1->SetDefaultDensity(rho_g);
    m1->SetDefaultSize(r_g);

    gen->SetStartTag(tag_particles);

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
    opengl::ChVisualSystemOpenGL vis;
    vis.AttachSystem(sys);
    vis.SetWindowTitle("Foam demo");
    vis.SetWindowSize(1280, 720);
    vis.SetRenderMode(opengl::WIREFRAME);
    vis.Initialize();
    vis.AddCamera(ChVector3d(4, -5, 4), ChVector3d(9, 0, 3));
    vis.SetCameraVertical(CameraVerticalDir::Z);
#endif

    // Perform the simulation
    double time = 0;
    int sim_frame = 0;
    int out_frame = 0;
    double exec_time = 0;
    std::ofstream ofile(out_file.c_str());

    while (time < time_end) {
        int numParticles = (int)sys->GetBodies().size() - 1;

        if (numParticles < maxNumParticles && sim_frame % gen_steps == 0) {
            SpawnParticles();
        }

        if (sim_frame % out_steps == 0) {
            char filename[100];
            sprintf(filename, "%s/data_%03d.dat", pov_dir.c_str(), out_frame + 1);
            utils::WriteVisualizationAssets(sys, filename);

            cout << " --------------------------------- Output frame:   " << out_frame + 1 << endl;
            cout << "                                   Sim frame:      " << sim_frame << endl;
            cout << "                                   Time:           " << time << endl;
            cout << "                                   Execution time: " << exec_time << endl;
            cout << "                                   Num. bodies:    " << numParticles << endl;

            ofile << sim_frame << "  " << time << "  " << exec_time << "  " << numParticles << "  "
                  << sys->GetNumContacts() << "\n";

            out_frame++;
        }

// Advance dynamics.
#ifdef CHRONO_OPENGL
        if (vis.Run()) {
            sys->DoStepDynamics(time_step);
            vis.Render();
        } else
            break;
#else
        sys->DoStepDynamics(time_step);
#endif

        time += time_step;
        sim_frame++;
        exec_time += sys->GetTimerStep();
    }

    // Final stats
    cout << "==================================" << endl;
    cout << "Number of bodies: " << sys->GetBodies().size() << endl;
    cout << "Simulation time: " << exec_time << endl;
    cout << "Number of threads: " << threads << endl;

    return 0;
}

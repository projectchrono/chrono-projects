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
// Author: Radu Serban, Arman Pazouki
// =============================================================================
//
// Chrono::Multicore demo program for low velocity cratering studies.
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
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsGenerators.h"
#include "chrono/utils/ChUtilsInputOutput.h"

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

enum ProblemType { SETTLING, DROPPING };
ProblemType problem = DROPPING;

enum PenetratorGeom { P_SPHERE, P_CONE1, P_CONE2 };
PenetratorGeom penetGeom = P_SPHERE;

// -----------------------------------------------------------------------------
// Global problem definitions
// -----------------------------------------------------------------------------

// Desired number of OpenMP threads (will be clamped to maximum available)
int threads = 2;

// Simulation parameters
double gravity = 9.81;

double time_settling_min = 0.1;
double time_settling_max = 0.8;
double time_dropping = 0.2;

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
ChSystemSMC::ContactForceModel contact_force_model = ChSystemSMC::ContactForceModel::Hertz;
ChSystemSMC::TangentialDisplacementModel tangential_displ_mode = ChSystemSMC::TangentialDisplacementModel::MultiStep;
#endif

// Output
bool povray_output = true;

#ifdef USE_SMC
const std::string out_dir = "../PENETRATOR_SMC";
#else
const std::string out_dir = "../PENETRATOR_NSC";
#endif

const std::string pov_dir = out_dir + "/POVRAY";
const std::string height_file = out_dir + "/height.dat";
const std::string stats_file = out_dir + "/stats.dat";
const std::string checkpoint_file = out_dir + "/settled.dat";

int out_fps_settling = 120;
int out_fps_dropping = 1200;

int timing_frame = -1;  // output detailed step timing at this frame

// Parameters for the granular material
int tag_particles = 0;
double r_g = 3e-3 / 2;
double rho_g = 2500;
double vol_g = (4.0 / 3) * CH_PI * r_g * r_g * r_g;
double mass_g = rho_g * vol_g;
ChVector3d inertia_g = 0.4 * mass_g * r_g * r_g * ChVector3d(1, 1, 1);

float Y_g = 1e8f;
float mu_g = 0.3f;
float cr_g = 0.1f;

// Parameters for the penetrator
double R_b = 2.54e-2 / 2;
double H_bc1 = 34.39e-3; // cone1 height
double R_bc1 = 9.215e-3; // cone1 radius
double H_bc2 = 22.10e-3; // cone1 height
double R_bc2 = 12.76e-3; // cone1 radius
double rho_b = 700;

float Y_b = 1.0e8f;
float mu_b = 0.4f;
float cr_b = 0.1f;

// Parameters for the containing bin
double hDimX = 4e-2;         // length in x direction
double hDimY = 4e-2;         // depth in y direction
double hDimZ = 7.5e-2;       // height in z direction
double hThickness = 0.5e-2;  // wall thickness

float Y_c = 2.0e6f;
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
int CreateObjects(ChSystemMulticore* sys) {
// Create the containing bin
#ifdef USE_SMC
    auto mat_c = chrono_types::make_shared<ChContactMaterialSMC>();
    mat_c->SetYoungModulus(Y_c);
    mat_c->SetFriction(mu_c);
    mat_c->SetRestitution(cr_c);

    utils::CreateBoxContainer(sys, mat_c, ChVector3d(hDimX, hDimY, hDimZ), hThickness);
#else
    auto mat_c = chrono_types::make_shared<ChContactMaterialNSC>();
    mat_c->SetFriction(mu_c);

    utils::CreateBoxContainer(sys, binId, mat_c, ChVector3d(hDimX, hDimY, hDimZ), hThickness);
#endif

// Create a material for the granular material
#ifdef USE_SMC
    auto mat_g = chrono_types::make_shared<ChContactMaterialSMC>();
    mat_g->SetYoungModulus(Y_g);
    mat_g->SetFriction(mu_g);
    mat_g->SetRestitution(cr_g);
#else
    auto mat_g = chrono_types::make_shared<ChContactMaterialNSC>();
    mat_g->SetFriction(mu_g);
#endif

    // Create a mixture entirely made out of spheres
    double r = 1.01 * r_g;
    utils::ChPDSampler<double> sampler(2 * r);
    utils::ChGenerator gen(sys);

    std::shared_ptr<utils::ChMixtureIngredient> m1 = gen.AddMixtureIngredient(utils::MixtureType::SPHERE, 1.0);
    m1->SetDefaultMaterial(mat_g);
    m1->SetDefaultDensity(rho_g);
    m1->SetDefaultSize(r_g);

    gen.SetStartTag(tag_particles);

    for (int i = 0; i < numLayers; i++) {
        double center = r + layerHeight / 2 + i * (2 * r + layerHeight);
        gen.CreateObjectsBox(sampler, ChVector3d(0, 0, center),
                             ChVector3d(hDimX - r, hDimY - r, layerHeight / 2));
        cout << "Layer " << i << "  total bodies: " << gen.GetTotalNumBodies() << endl;
    }

    return gen.GetTotalNumBodies();
}

// -----------------------------------------------------------------------------
// Calculate intertia properties of the falling object
// -----------------------------------------------------------------------------
void CalculatePenetratorInertia(double & mass, ChVector3d & inertia) {
    ChVector3d gyr_b;  // components gyration
    double vol_b;      // components volume
    switch (penetGeom) {
        case P_SPHERE:
            vol_b = ChSphere::GetVolume(R_b);
            gyr_b = ChSphere::GetGyration(R_b).diagonal();
            mass = rho_b * vol_b;
            inertia = mass * gyr_b;
            break;
        case P_CONE1:
            // apex angle = 30 de
            vol_b = ChCone::GetVolume(R_bc1, H_bc1);
            gyr_b = ChCone::GetGyration(R_bc1, H_bc1).diagonal();
            mass = rho_b * vol_b;
            inertia = mass * gyr_b;
            break;
        case P_CONE2:
            // apex angle = 60 deg
            vol_b = ChCone::GetVolume(R_bc2, H_bc2);
            gyr_b = ChCone::GetGyration(R_bc2, H_bc2).diagonal();
            mass = rho_b * vol_b;
            inertia = mass * gyr_b;
            break;
	}
}

// -----------------------------------------------------------------------------
// Create collision geometry of the falling object
// -----------------------------------------------------------------------------
void CreatePenetratorGeometry(std::shared_ptr<ChBody> obj, std::shared_ptr<ChContactMaterial> mat) {
    switch (penetGeom) {
        case P_SPHERE:
            utils::AddSphereGeometry(obj.get(), mat, R_b);
            break;
        case P_CONE1:
            // apex angle = 30 de
            utils::AddConeGeometry(obj.get(), mat, R_bc1, H_bc1);
            break;
        case P_CONE2:
            // apex angle = 60 deg
            utils::AddConeGeometry(obj.get(), mat, R_bc2, H_bc2);
            break;
    }
}

// -----------------------------------------------------------------------------
// Calculate falling object height
// -----------------------------------------------------------------------------
double RecalcPenetratorLocation(double z) {
    double locZ = 0;
    switch (penetGeom) {
        case P_SPHERE:
            locZ = z + R_b + r_g;
            break;
        case P_CONE1:
            // apex angle = 30 de
            locZ = z + 0.75 * H_bc1 + r_g;
            break;
        case P_CONE2:
            // apex angle = 60 deg
            locZ = z + 0.75 * H_bc2 + r_g;
            break;
    }
    return locZ;
}

// -----------------------------------------------------------------------------
// Find the height of the highest and lowest, respectively, sphere in the
// granular mix, respectively.  We only look at bodies whith stricty positive
// identifiers (to exclude the containing bin).
// -----------------------------------------------------------------------------
double FindHighest(ChSystem* sys) {
    double highest = 0;
    for (auto body : sys->GetBodies()) {
        if (body->GetTag() >= tag_particles && body->GetPos().z() > highest)
            highest = body->GetPos().z();
    }
    return highest;
}

double FindLowest(ChSystem* sys) {
    double lowest = 1000;
    for (auto body : sys->GetBodies()) {
        if (body->GetTag() >= tag_particles && body->GetPos().z() < lowest)
            lowest = body->GetPos().z();
    }
    return lowest;
}

// -----------------------------------------------------------------------------
// Create the falling object such that its bottom point is at the specified height
// and its downward initial velocity has the specified magnitude.
// -----------------------------------------------------------------------------
std::shared_ptr<ChBody> CreatePenetrator(ChSystemMulticore* sys) {
    // Estimate object initial location and velocity
    double z = FindHighest(sys);
    double vz = std::sqrt(2 * gravity * h);
    double initLoc = RecalcPenetratorLocation(z);
    cout << "creating object at " << initLoc << " and velocity " << vz << endl;

// Create a material for the penetrator
#ifdef USE_SMC
    auto mat = chrono_types::make_shared<ChContactMaterialSMC>();
    mat->SetYoungModulus(Y_b);
    mat->SetFriction(mu_b);
    mat->SetRestitution(cr_b);
#else
    auto mat = chrono_types::make_shared<ChContactMaterialNSC>();
    mat->SetFriction(mu_b);
#endif

    // Create the falling object
    auto obj = chrono_types::make_shared<ChBody>();

    double mass;
    ChVector3d inertia;
    CalculatePenetratorInertia(mass, inertia);
    obj->SetMass(mass);
    obj->SetInertiaXX(inertia);
    obj->SetPos(ChVector3d(0, 0, initLoc));
    obj->SetRot(QuatFromAngleAxis(-CH_PI / 2, VECT_X));
    obj->SetLinVel(ChVector3d(0, 0, -vz));
    obj->EnableCollision(true);
    obj->SetFixed(false);

    CreatePenetratorGeometry(obj, mat);

    sys->AddBody(obj);
    return obj;
}

// -----------------------------------------------------------------------------
// Return true if all bodies in the granular mix have a linear velocity whose
// magnitude is below the specified value.
// -----------------------------------------------------------------------------
bool CheckSettled(ChSystem* sys, double threshold) {
    double t2 = threshold * threshold;

    for (auto body : sys->GetBodies()) {
        if (body->GetTag() >= tag_particles) {
            double vel2 = body->GetLinVel().Length2();
            if (vel2 > t2)
                return false;
        }
    }

    return true;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void SetArgumentsForMbdFromInput(int argc, char* argv[]) {
    if (argc > 1) {
        const char* text = argv[1];
        rho_b = atof(text);
    }
    int pType = 0;
    if (argc > 2) {
        const char* text = argv[2];
        pType = atoi(text);
        switch (pType) {
            case 0:
                penetGeom = P_SPHERE;
                break;
            case 1:
                penetGeom = P_CONE1;
                break;
            case 2:
                penetGeom = P_CONE2;
                break;
        }
    }

    const std::string simSettings = out_dir + "/simSetting.txt";
    std::ofstream sFile(simSettings.c_str());
    sFile << "density: " << rho_b << " penetrometer shape: " << pType << " (0: sphere, 1: cone1, 2: cone2) \n";
}

// -----------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    // Create output directories.
    if (!filesystem::create_directory(filesystem::path(out_dir))) {
        cout << "Error creating directory " << out_dir << endl;
        return 1;
    }
    if (!filesystem::create_directory(filesystem::path(pov_dir))) {
        cout << "Error creating directory " << pov_dir << endl;
        return 1;
    }
    
    // Get problem parameters from arguments
    SetArgumentsForMbdFromInput(argc, argv);

// Create system
#ifdef USE_SMC
    cout << "Create SMC system" << endl;
    ChSystemMulticoreSMC* sys = new ChSystemMulticoreSMC();
#else
    cout << "Create NSC system" << endl;
    ChSystemMulticoreNSC* sys = new ChSystemMulticoreNSC();
#endif

    sys->SetCollisionSystemType(ChCollisionSystem::Type::MULTICORE);

    // Debug log messages.
    ////sys->SetLoggingLevel(LOG_INFO, true);
    ////sys->SetLoggingLevel(LOG_TRACE, true);

    // Set number of threads.
    int max_threads = omp_get_num_procs();
    if (threads > max_threads)
        threads = max_threads;
    sys->SetNumThreads(threads);
    cout << "Using " << threads << " threads" << endl;

    // Set gravitational acceleration
    sys->SetGravitationalAcceleration(ChVector3d(0, 0, -gravity));

    // Edit system settings
    sys->GetSettings()->solver.use_full_inertia_tensor = false;
    sys->GetSettings()->solver.tolerance = tolerance;

#ifdef USE_SMC
    sys->GetSettings()->collision.narrowphase_algorithm = ChNarrowphase::Algorithm::HYBRID;

    sys->GetSettings()->solver.contact_force_model = contact_force_model;
    sys->GetSettings()->solver.tangential_displ_mode = tangential_displ_mode;
#else
    sys->GetSettings()->solver.solver_mode = SolverMode::SLIDING;
    sys->GetSettings()->solver.max_iteration_normal = max_iteration_normal;
    sys->GetSettings()->solver.max_iteration_sliding = max_iteration_sliding;
    sys->GetSettings()->solver.max_iteration_spinning = max_iteration_spinning;
    sys->GetSettings()->solver.alpha = 0;
    sys->GetSettings()->solver.contact_recovery_speed = contact_recovery_speed;
    sys->ChangeSolverType(SolverType::APGD);

    sys->GetSettings()->collision.collision_envelope = 0.05 * r_g;
#endif

    sys->GetSettings()->collision.bins_per_axis = vec3(20, 20, 20);

    // Depending on problem type:
    // - Select end simulation time
    // - Select output FPS
    // - Create granular material and container
    // - Create falling object
    double time_end;
    int out_fps;
    std::shared_ptr<ChBody> obj;

    if (problem == SETTLING) {
        time_end = time_settling_max;
        out_fps = out_fps_settling;

        cout << "Create granular material" << endl;
        CreateObjects(sys);
    }

    if (problem == DROPPING) {
        time_end = time_dropping;
        out_fps = out_fps_dropping;

        // Create the granular material and the container from the checkpoint file.
        cout << "Read checkpoint data from " << checkpoint_file;
        utils::ReadCheckpoint(sys, checkpoint_file);
        cout << "  done.  Read " << sys->GetBodies().size() << " bodies." << endl;
        obj = CreatePenetrator(sys);
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
    std::ofstream sfile(stats_file.c_str());
    std::ofstream hfile(height_file.c_str());

#ifdef CHRONO_OPENGL
    opengl::ChVisualSystemOpenGL vis;
    vis.AttachSystem(sys);
    vis.SetWindowTitle("Penetrator Test");
    vis.SetWindowSize(1280, 720);
    vis.SetRenderMode(opengl::WIREFRAME);
    vis.Initialize();
    vis.AddCamera(ChVector3d(0, -10 * hDimY, hDimZ), ChVector3d(0, 0, hDimZ));
    vis.SetCameraVertical(CameraVerticalDir::Z);
#endif

    while (time < time_end) {
        if (sim_frame == next_out_frame) {
            cout << endl;
            cout << "---- Frame:          " << out_frame << endl;
            cout << "     Sim frame:      " << sim_frame << endl;
            cout << "     Time:           " << time << endl;
            cout << "     Lowest point:   " << FindLowest(sys) << endl;
            cout << "     Avg. contacts:  " << num_contacts / out_steps << endl;
            cout << "     Execution time: " << exec_time << endl;

            sfile << time << "  " << exec_time << "  " << num_contacts / out_steps << "\n";

            // If enabled, output data for PovRay postprocessing.
            if (povray_output) {
                char filename[100];
                sprintf(filename, "%s/data_%03d.dat", pov_dir.c_str(), out_frame + 1);
                utils::WriteVisualizationAssets(sys, filename, false);
            }

            // Create a checkpoint from the current state.
            if (problem == SETTLING) {
                cout << "     Write checkpoint data " << flush;
//                utils::WriteCheckpoint(sys, checkpoint_file);
                cout << sys->GetBodies().size() << " bodies" << endl;
            }

            // Save current penetrator height.
            if (problem == DROPPING) {
                hfile << time << "  " << obj->GetPos().z() << "\n";
                hfile.flush();
                cout << "     Penetrator height:    " << obj->GetPos().z() << endl;
            }

            out_frame++;
            next_out_frame += out_steps;
            num_contacts = 0;
        }

        if (problem == SETTLING && time > time_settling_min && CheckSettled(sys, zero_v)) {
            cout << "Granular material settled...  time = " << time << endl;
            break;
        }

// Advance simulation by one step
#ifdef CHRONO_OPENGL
        if (vis.Run()) {
            sys->DoStepDynamics(time_step);
            vis.Render();
        } else
            break;
#else
        sys->DoStepDynamics(time_step);
#endif

        ////progressbar(out_steps + sim_frame - next_out_frame + 1, out_steps);
        //TimingOutput(sys);

        time += time_step;
        sim_frame++;
        exec_time += sys->GetTimerStep();
        num_contacts += sys->GetNumContacts();

        // If requested, output detailed timing information for this step
        if (sim_frame == timing_frame)
            sys->PrintStepStats();
    }

    // Create a checkpoint from the last state
    if (problem == SETTLING) {
        cout << "Write checkpoint data to " << checkpoint_file;
        utils::WriteCheckpoint(sys, checkpoint_file);
        cout << "  done.  Wrote " << sys->GetBodies().size() << " bodies." << endl;
    }

    // Final stats
    cout << "==================================" << endl;
    cout << "Number of bodies:  " << sys->GetBodies().size() << endl;
    cout << "Lowest position:   " << FindLowest(sys) << endl;
    cout << "Simulation time:   " << exec_time << endl;
    cout << "Number of threads: " << threads << endl;

    return 0;
}

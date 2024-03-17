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
// Authors: Radu Serban
// =============================================================================
//
// Demo program for simulating a rollover testing rig.
//
// =============================================================================

#include <cstdio>
#include <vector>
#include <cmath>

#include "chrono/ChConfig.h"
#include "chrono/assets/ChVisualShapeBox.h"
#include "chrono/assets/ChVisualShapeCapsule.h"
#include "chrono/utils/ChUtilsGeometry.h"
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

using namespace chrono;

using std::cout;
using std::endl;
using std::flush;

// -----------------------------------------------------------------------------
// Problem setup
// -----------------------------------------------------------------------------

// Specify solution method (comment next line for NSC)
#define USE_SMC

// Simulation phase
enum ProblemType { SETTLING, PUSHING, TESTING };

ProblemType problem = TESTING;

// Wheel contact shape (one of collision::CYLINDER or collision::ROUNDEDCYL)
ChCollisionShape::Type wheel_shape = ChCollisionShape::Type::CYLINDER;

// -----------------------------------------------------------------------------
// Simulation parameters
// -----------------------------------------------------------------------------

// Desired number of OpenMP threads (will be clamped to maximum available)
int threads = 20;

// Simulation duration.
double time_settling = 5;
double time_pushing = 2;

// Solver parameters
bool clamp_bilaterals = false;
double bilateral_clamp_speed = 1000;
double tolerance = 1e-3;

#ifdef USE_SMC
double time_step = 1e-4;
int max_iteration_bilateral = 100;
#else
double time_step = 5e-4;
int max_iteration_normal = 0;
int max_iteration_sliding = 1000;
int max_iteration_spinning = 0;
int max_iteration_bilateral = 0;
float contact_recovery_speed = 1e4;
#endif

// Output
#ifdef USE_SMC
const std::string out_dir = "../SCHMUTZ_SMC";
#else
const std::string out_dir = "../SCHMUTZ_NSC";
#endif
const std::string pov_dir = out_dir + "/POVRAY";
const std::string checkpoint_file = out_dir + "/settled.dat";
const std::string stats_file = out_dir + "/stats.dat";
const std::string results_file = out_dir + "/results.dat";

int out_fps_settling = 30;
int out_fps_pushing = 60;

// -----------------------------------------------------------------------------
// Parameters for the granular material (identical spheres)
// -----------------------------------------------------------------------------
double r_g = 0.01;
double rho_g = 2700;
float Y_g = 5e7;
float cr_g = 0.1f;
float mu_g = 0.4f;
unsigned int desired_num_particles = 10000;

// -----------------------------------------------------------------------------
// Parameters for the test rig
// -----------------------------------------------------------------------------
double mass1 = 550;
double mass_wheel = 351;

ChVector3d inertia_sled(100, 100, 100);
ChVector3d inertia_wheel(50, 138, 138);

double a = 0.951;
double b = 0.169;
double c = 0.8;

double e = 1;

double r_w = 0.290;
double w_w = 0.200;
double s_w = 0.058;

double L = 4.0;
double W = 0.8;
double H = 2.0;

double d = L / 2 - 1.2 * w_w;

double init_vel = 5;
double init_angle = (CH_C_PI / 180) * 4;

// =============================================================================
// This class encapsulates the rig's mechanism
// =============================================================================
class Mechanism {
  public:
    Mechanism(ChSystemMulticore* system, double h);

    const ChVector3d& GetSledVelocity() const { return m_sled->GetLinVel(); }
    const ChVector3d& GetWheelVelocity() const { return m_wheel->GetLinVel(); }

    void WriteResults(std::ofstream& f, double time);

  private:
    ChVector3d calcLocationWheel(double h) {
        double ca = std::cos(init_angle);
        double sa = std::sin(init_angle);

        return ChVector3d(-d - ca * (c + w_w / 2) + sa * (b + r_w), 0, h + sa * (c + w_w / 2) + ca * (b + r_w));
    }

    ChVector3d calcLocationRevolute(double h) {
        double ca = std::cos(init_angle);
        double sa = std::sin(init_angle);

        return ChVector3d(-d - ca * (a + c + w_w / 2) + sa * r_w, 0, h + sa * (a + c + w_w / 2) + ca * r_w);
    }

    std::shared_ptr<ChBody> m_ground;
    std::shared_ptr<ChBody> m_sled;
    std::shared_ptr<ChBody> m_wheel;

    std::shared_ptr<ChLinkLockPrismatic> m_prismatic;
    std::shared_ptr<ChLinkLockRevolute> m_revolute;
};

Mechanism::Mechanism(ChSystemMulticore* system, double h) {
    // Calculate hardpoint locations at initial configuration (expressed in the
    // global frame)
    ChVector3d loc_revolute = calcLocationRevolute(h);
    ChVector3d loc_wheel = calcLocationWheel(h);
    ChVector3d loc_sled = loc_revolute - ChVector3d(e, 0, 0);
    ChVector3d loc_prismatic = loc_sled - ChVector3d(0, 0, e / 4);

    // Create the ground body
    m_ground = chrono_types::make_shared<ChBody>();
    m_ground->SetIdentifier(-1);
    m_ground->SetFixed(true);
    m_ground->EnableCollision(false);

    system->AddBody(m_ground);

    // Create the sled body
    m_sled = chrono_types::make_shared<ChBody>();
    m_sled->SetIdentifier(1);
    m_sled->SetMass(mass1);
    m_sled->SetInertiaXX(inertia_sled);
    m_sled->SetPos(loc_sled);
    m_sled->SetLinVel(ChVector3d(init_vel, 0, 0));
    m_sled->SetFixed(false);
    m_sled->EnableCollision(false);

    auto box_sled = chrono_types::make_shared<ChVisualShapeBox>(2.0 * e, (2.0 / 3) * e, (2.0 / 3) * e);
    box_sled->SetColor(ChColor(0.7f, 0.3f, 0.3f));
    m_sled->AddVisualShape(box_sled);

    system->AddBody(m_sled);

    // Create a material for the wheel body
#ifdef USE_SMC
    auto mat_w = chrono_types::make_shared<ChContactMaterialSMC>();
    mat_w->SetYoungModulus(2e6f);
    mat_w->SetFriction(0.4f);
    mat_w->SetRestitution(0.1f);
#else
    auto mat_w = chrono_types::make_shared<ChContactMaterialNSC>();
    mat_w->SetFriction(0.4f);
#endif

    // Create the wheel body
    m_wheel = chrono_types::make_shared<ChBody>();
    m_wheel->SetIdentifier(2);
    m_wheel->SetMass(mass_wheel);
    m_wheel->SetInertiaXX(inertia_wheel);
    m_wheel->SetPos(loc_wheel);
    m_wheel->SetRot(QuatFromAngleY(init_angle));
    m_wheel->SetLinVel(ChVector3d(init_vel, 0, 0));
    m_wheel->SetFixed(false);
    m_wheel->EnableCollision(true);

    switch (wheel_shape) {
        case ChCollisionShape::Type::CYLINDER:
            utils::AddCylinderGeometry(m_wheel.get(), mat_w, r_w, w_w / 2, ChVector3d(c, 0, -b), QuatFromAngleY(CH_C_PI_2));
            break;
        case ChCollisionShape::Type::ROUNDEDCYL:
            utils::AddRoundedCylinderGeometry(m_wheel.get(), mat_w, r_w - s_w, w_w / 2 - s_w, s_w, ChVector3d(c, 0, -b),
                                              QuatFromAngleY(CH_C_PI_2));
            break;
    }

    auto cap_wheel = chrono_types::make_shared<ChVisualShapeCapsule>(w_w / 4, a + c - w_w / 2);
    cap_wheel->SetColor(ChColor(0.3f, 0.3f, 0.7f));
    m_wheel->AddVisualShape(cap_wheel, ChFrame<>(ChVector3d((c - a) / 2, 0, -b), QuatFromAngleY(CH_C_PI_2)));

    system->AddBody(m_wheel);

    // Create and initialize translational joint ground - sled
    m_prismatic = chrono_types::make_shared<ChLinkLockPrismatic>();
    m_prismatic->Initialize(m_ground, m_sled, ChFrame<>(loc_prismatic, QuatFromAngleY(CH_C_PI_2)));
    system->AddLink(m_prismatic);

    // Create and initialize revolute joint sled - wheel
    m_revolute = chrono_types::make_shared<ChLinkLockRevolute>();
    m_revolute->Initialize(m_wheel, m_sled, ChFrame<>(loc_revolute, QuatFromAngleX(CH_C_PI_2)));
    system->AddLink(m_revolute);
}

void Mechanism::WriteResults(std::ofstream& f, double time) {
    // Velocity of sled body (in absolute frame)
    ChVector3d sled_vel = m_sled->GetLinVel();

    // Velocity of wheel body (in absolute frame)
    ChVector3d wheel_vel = m_wheel->GetLinVel();

    // Coordinate system of the revolute joint (relative to the frame of body2, in
    // this case the sled)
    auto revFrame = m_revolute->GetFrame2Rel();

    // Reaction force in revolute joint
    ChVector3d force_jointsys = m_revolute->GetReaction2().force;
    ChVector3d force_bodysys = revFrame.TransformDirectionLocalToParent(force_jointsys);
    ChVector3d force_abssys = m_sled->GetCoordsys().TransformDirectionLocalToParent(force_bodysys);

    f << time << "  " << sled_vel.x() << "  " << sled_vel.y() << "  " << sled_vel.z() << "      " << wheel_vel.x() << "  "
      << wheel_vel.y() << "  " << wheel_vel.z() << "      " << force_abssys.x() << "  " << force_abssys.y() << "  "
      << force_abssys.z() << "\n";
}

// =============================================================================
// Create container bin.
// =============================================================================
void CreateContainer(ChSystemMulticore* system) {
    int id_c = -200;
    double thickness = 0.2;

#ifdef USE_SMC
    auto mat_c = chrono_types::make_shared<ChContactMaterialSMC>();
    mat_c->SetYoungModulus(2e6f);
    mat_c->SetFriction(0.4f);
    mat_c->SetRestitution(0.1f);

    utils::CreateBoxContainer(system, id_c, mat_c, ChVector3d(L / 2, W / 2, H / 2), thickness / 2);
#else
    auto mat_c = chrono_types::make_shared<ChContactMaterialNSC>();
    mat_c->SetFriction(0.4f);

    utils::CreateBoxContainer(system, id_c, mat_c, ChVector3d(L / 2, W / 2, H / 2), thickness / 2);

#endif
}

// =============================================================================
// Create granular material.
// =============================================================================
void CreateParticles(ChSystemMulticore* system) {
// Create a material for the ball mixture.
#ifdef USE_SMC
    auto mat_g = chrono_types::make_shared<ChContactMaterialSMC>();
    mat_g->SetYoungModulus(Y_g);
    mat_g->SetFriction(mu_g);
    mat_g->SetRestitution(cr_g);
#else
    auto mat_g = chrono_types::make_shared<ChContactMaterialNSC>();
    mat_g->SetFriction(mu_g);
#endif

    // Create a mixture entirely made out of spheres.
    double r = 1.01 * r_g;
    utils::PDSampler<double> sampler(2 * r);
    utils::Generator gen(system);

    std::shared_ptr<utils::MixtureIngredient> m1 = gen.AddMixtureIngredient(utils::MixtureType::SPHERE, 1.0);
    m1->setDefaultMaterial(mat_g);
    m1->setDefaultDensity(rho_g);
    m1->setDefaultSize(r_g);

    // Create particles, one layer at a time, until the desired number is reached.
    gen.setBodyIdentifier(100);

    ChVector3d hdims(L / 2 - r, W / 2 - r, 0);
    ChVector3d center(0, 0, 2 * r);

    int layer = 1;
    while (gen.getTotalNumBodies() < desired_num_particles) {
        gen.CreateObjectsBox(sampler, center, hdims);
        cout << "layer " << layer << "    total particles: " << gen.getTotalNumBodies() << endl;
        center.z() += 2 * r;
        layer++;
    }
}

// =============================================================================
// Find the height of the highest and lowest sphere in the granular mix.
// We only look at bodies whith identifiers larger than 100 (to exclude all
// other bodies).
// =============================================================================
void FindRange(ChSystem* sys, double& lowest, double& highest) {
    highest = -1000;
    lowest = 1000;
    for (auto body : sys->GetBodies()) {
        if (body->GetIdentifier() < 100)
            continue;
        double h = body->GetPos().z();
        if (h < lowest)
            lowest = h;
        if (h > highest)
            highest = h;
    }
}

// =============================================================================
// =============================================================================
int main(int argc, char* argv[]) {
    // --------------------------
    // Create output directories.
    // --------------------------

    if (!filesystem::create_directory(filesystem::path(out_dir))) {
        cout << "Error creating directory " << out_dir << endl;
        return 1;
    }
    if (!filesystem::create_directory(filesystem::path(pov_dir))) {
        cout << "Error creating directory " << pov_dir << endl;
        return 1;
    }

// --------------
// Create system.
// --------------

#ifdef USE_SMC
    cout << "Create SMC system" << endl;
    ChSystemMulticoreSMC* sys = new ChSystemMulticoreSMC();
#else
    cout << "Create NSC system" << endl;
    ChSystemMulticoreNSC* sys = new ChSystemMulticoreNSC();
#endif

    sys->SetCollisionSystemType(ChCollisionSystem::Type::MULTICORE);

    sys->SetGravitationalAcceleration(ChVector3d(0, 0, -9.81));

    // ----------------------
    // Set number of threads.
    // ----------------------

    int max_threads = omp_get_num_procs();
    if (threads > max_threads)
        threads = max_threads;
    sys->SetNumThreads(threads);
    cout << "Using " << threads << " threads" << endl;

    // ---------------------
    // Edit system settings.
    // ---------------------

    sys->GetSettings()->solver.tolerance = tolerance;
    sys->GetSettings()->solver.max_iteration_bilateral = max_iteration_bilateral;
    sys->GetSettings()->solver.clamp_bilaterals = clamp_bilaterals;
    sys->GetSettings()->solver.bilateral_clamp_speed = bilateral_clamp_speed;

#ifdef USE_SMC
    sys->GetSettings()->collision.narrowphase_algorithm = ChNarrowphase::Algorithm::PRIMS;
#else
    sys->GetSettings()->solver.solver_mode = SolverMode::SLIDING;
    sys->GetSettings()->solver.max_iteration_normal = max_iteration_normal;
    sys->GetSettings()->solver.max_iteration_sliding = max_iteration_sliding;
    sys->GetSettings()->solver.max_iteration_spinning = max_iteration_spinning;
    sys->GetSettings()->solver.alpha = 0;
    sys->GetSettings()->solver.contact_recovery_speed = contact_recovery_speed;
    sys->ChangeSolverType(SolverType::APGDREF);

    sys->GetSettings()->collision.collision_envelope = 0.05 * r_g;
#endif

    sys->GetSettings()->collision.bins_per_axis = vec3(10, 10, 10);

    // ----------------------------------------
    // Depending on problem type:
    // - Select end simulation time
    // - Create granular material and container
    // - Create mechanism
    // ----------------------------------------

    double time_end;
    int out_fps;
    Mechanism* mech = NULL;

    switch (problem) {
        case SETTLING:
            time_end = time_settling;
            out_fps = out_fps_settling;

            cout << "Create granular material" << endl;
            CreateContainer(sys);
            CreateParticles(sys);

            break;

        case PUSHING: {
            time_end = time_pushing;
            out_fps = out_fps_pushing;

            // Create the granular material and the container from the checkpoint file.
            cout << "Read checkpoint data from " << checkpoint_file;
            utils::ReadCheckpoint(sys, checkpoint_file);
            cout << "  done.  Read " << sys->GetBodies().size() << " bodies." << endl;

            // Create the mechanism with the wheel just above the granular material.
            double lowest, highest;
            FindRange(sys, lowest, highest);
            cout << "Create mechanism above height " << highest + r_g << endl;
            mech = new Mechanism(sys, highest + r_g);
        }

        break;

        case TESTING:
            time_end = time_pushing;
            out_fps = out_fps_pushing;

            mech = new Mechanism(sys, 0.9 * H);

            break;
    }

    // Number of steps.
    int num_steps = (int)std::ceil(time_end / time_step);
    int out_steps = (int)std::ceil((1.0 / time_step) / out_fps);

    // -----------------------
    // Perform the simulation.
    // -----------------------
    double time = 0;
    int sim_frame = 0;
    int out_frame = 0;
    int next_out_frame = 0;
    double exec_time = 0;
    int num_contacts = 0;
    std::ofstream sfile(stats_file.c_str());
    std::ofstream rfile(results_file.c_str());

    // Disable buffering on output file streams.
    sfile.rdbuf()->pubsetbuf(0, 0);
    rfile.rdbuf()->pubsetbuf(0, 0);

#ifdef CHRONO_OPENGL
    opengl::ChVisualSystemOpenGL vis;
    vis.AttachSystem(sys);
    vis.SetWindowTitle("Pressure Sinkage Test");
    vis.SetWindowSize(1280, 720);
    vis.SetRenderMode(opengl::WIREFRAME);
    vis.Initialize();
    vis.AddCamera(ChVector3d(0, -8, 0), ChVector3d(0, 0, 0));
    vis.SetCameraVertical(CameraVerticalDir::Z);
#endif

    while (time < time_end) {
        if (sim_frame == next_out_frame) {
            char filename[100];
            sprintf(filename, "%s/data_%03d.dat", pov_dir.c_str(), out_frame + 1);
            utils::WriteVisualizationAssets(sys, filename, problem == TESTING);
            double highest, lowest;
            FindRange(sys, lowest, highest);
            cout << "------------ Output frame:   " << out_frame + 1 << endl;
            cout << "             Sim frame:      " << sim_frame << endl;
            cout << "             Time:           " << time << endl;
            cout << "             Lowest point:   " << lowest << endl;
            cout << "             Highest point:  " << highest << endl;
            if (problem != SETTLING) {
                cout << "             Sled X vel.   : " << mech->GetSledVelocity().x() << endl;
                cout << "             Wheel X vel.  : " << mech->GetWheelVelocity().x() << endl;
            }
            cout << "             Avg. contacts:  " << num_contacts / out_steps << endl;
            cout << "             Execution time: " << exec_time << endl;

            sfile << time << "  " << exec_time << "  " << num_contacts / out_steps << "\n";
            sfile.flush();

            // Create a checkpoint from the current state.
            if (problem == SETTLING) {
                cout << "             Write checkpoint data " << flush;
                utils::WriteCheckpoint(sys, checkpoint_file);
                cout << sys->GetBodies().size() << " bodies" << endl;
            }

            out_frame++;
            next_out_frame += out_steps;
            num_contacts = 0;
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

        // Increment counters
        time += time_step;
        sim_frame++;
        exec_time += sys->GetTimerStep();
        num_contacts += sys->GetNumContacts();

        // Save results
        if (problem != SETTLING) {
            mech->WriteResults(rfile, time);
            rfile.flush();
        }
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
    cout << "Simulation time:   " << exec_time << endl;
    cout << "Number of threads: " << threads << endl;

    return 0;
}

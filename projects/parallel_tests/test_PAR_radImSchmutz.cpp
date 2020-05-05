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
#include "chrono/core/ChStream.h"
#include "chrono/assets/ChBoxShape.h"
#include "chrono/assets/ChCapsuleShape.h"
#include "chrono/assets/ChColorAsset.h"
#include "chrono/utils/ChUtilsGeometry.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsGenerators.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_parallel/solver/ChSystemDescriptorParallel.h"
#include "chrono_parallel/collision/ChNarrowphaseRUtils.h"

#include "chrono_thirdparty/filesystem/path.h"

// Control use of OpenGL run-time rendering
// Note: CHRONO_OPENGL is defined in ChConfig.h
//#undef CHRONO_OPENGL

#ifdef CHRONO_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif

using namespace chrono;
using namespace chrono::collision;

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

// Perform dynamic tuning of number of threads?
bool thread_tuning = false;

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

ChVector<> inertia_sled(100, 100, 100);
ChVector<> inertia_wheel(50, 138, 138);

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
    Mechanism(ChSystemParallel* system, double h);

    const ChVector<>& GetSledVelocity() const { return m_sled->GetPos_dt(); }
    const ChVector<>& GetWheelVelocity() const { return m_wheel->GetPos_dt(); }

    void WriteResults(ChStreamOutAsciiFile& f, double time);

  private:
    ChVector<> calcLocationWheel(double h) {
        double ca = std::cos(init_angle);
        double sa = std::sin(init_angle);

        return ChVector<>(-d - ca * (c + w_w / 2) + sa * (b + r_w), 0, h + sa * (c + w_w / 2) + ca * (b + r_w));
    }

    ChVector<> calcLocationRevolute(double h) {
        double ca = std::cos(init_angle);
        double sa = std::sin(init_angle);

        return ChVector<>(-d - ca * (a + c + w_w / 2) + sa * r_w, 0, h + sa * (a + c + w_w / 2) + ca * r_w);
    }

    std::shared_ptr<ChBody> m_ground;
    std::shared_ptr<ChBody> m_sled;
    std::shared_ptr<ChBody> m_wheel;

    std::shared_ptr<ChLinkLockPrismatic> m_prismatic;
    std::shared_ptr<ChLinkLockRevolute> m_revolute;
};

Mechanism::Mechanism(ChSystemParallel* system, double h) {
    // Calculate hardpoint locations at initial configuration (expressed in the
    // global frame)
    ChVector<> loc_revolute = calcLocationRevolute(h);
    ChVector<> loc_wheel = calcLocationWheel(h);
    ChVector<> loc_sled = loc_revolute - ChVector<>(e, 0, 0);
    ChVector<> loc_prismatic = loc_sled - ChVector<>(0, 0, e / 4);

    // Create the ground body
    m_ground = chrono_types::make_shared<ChBody>(chrono_types::make_shared<ChCollisionModelParallel>());
    m_ground->SetIdentifier(-1);
    m_ground->SetBodyFixed(true);
    m_ground->SetCollide(false);

    system->AddBody(m_ground);

    // Create the sled body
    m_sled = chrono_types::make_shared<ChBody>(chrono_types::make_shared<ChCollisionModelParallel>());
    m_sled->SetIdentifier(1);
    m_sled->SetMass(mass1);
    m_sled->SetInertiaXX(inertia_sled);
    m_sled->SetPos(loc_sled);
    m_sled->SetPos_dt(ChVector<>(init_vel, 0, 0));
    m_sled->SetBodyFixed(false);
    m_sled->SetCollide(false);

    auto box_sled = chrono_types::make_shared<ChBoxShape>();
    box_sled->GetBoxGeometry().Size = ChVector<>(e, e / 3, e / 3);
    box_sled->Pos = ChVector<>(0, 0, 0);
    box_sled->Rot = ChQuaternion<>(1, 0, 0, 0);
    m_sled->AddAsset(box_sled);

    auto col_sled = chrono_types::make_shared<ChColorAsset>();
    col_sled->SetColor(ChColor(0.7f, 0.3f, 0.3f));
    m_sled->AddAsset(col_sled);

    system->AddBody(m_sled);

    // Create a material for the wheel body
#ifdef USE_SMC
    auto mat_w = chrono_types::make_shared<ChMaterialSurfaceSMC>();
    mat_w->SetYoungModulus(2e6f);
    mat_w->SetFriction(0.4f);
    mat_w->SetRestitution(0.1f);
#else
    auto mat_w = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    mat_w->SetFriction(0.4f);
#endif

    // Create the wheel body
    m_wheel = chrono_types::make_shared<ChBody>(chrono_types::make_shared<ChCollisionModelParallel>());
    m_wheel->SetIdentifier(2);
    m_wheel->SetMass(mass_wheel);
    m_wheel->SetInertiaXX(inertia_wheel);
    m_wheel->SetPos(loc_wheel);
    m_wheel->SetRot(Q_from_AngY(init_angle));
    m_wheel->SetPos_dt(ChVector<>(init_vel, 0, 0));
    m_wheel->SetBodyFixed(false);
    m_wheel->SetCollide(true);

    m_wheel->GetCollisionModel()->ClearModel();
    switch (wheel_shape) {
        case ChCollisionShape::Type::CYLINDER:
            utils::AddCylinderGeometry(m_wheel.get(), mat_w, r_w, w_w / 2, ChVector<>(c, 0, -b), Q_from_AngZ(CH_C_PI_2));
            break;
        case ChCollisionShape::Type::ROUNDEDCYL:
            utils::AddRoundedCylinderGeometry(m_wheel.get(), mat_w, r_w - s_w, w_w / 2 - s_w, s_w, ChVector<>(c, 0, -b),
                                              Q_from_AngZ(CH_C_PI_2));
            break;
    }
    m_wheel->GetCollisionModel()->BuildModel();

    auto cap_wheel = chrono_types::make_shared<ChCapsuleShape>();
    cap_wheel->GetCapsuleGeometry().hlen = (a + c) / 2 - w_w / 4;
    cap_wheel->GetCapsuleGeometry().rad = w_w / 4;
    cap_wheel->Pos = ChVector<>((c - a) / 2, 0, -b);
    cap_wheel->Rot = Q_from_AngZ(CH_C_PI_2);
    m_wheel->AddAsset(cap_wheel);

    auto col_wheel = chrono_types::make_shared<ChColorAsset>();
    col_wheel->SetColor(ChColor(0.3f, 0.3f, 0.7f));
    m_wheel->AddAsset(col_wheel);

    system->AddBody(m_wheel);

    // Create and initialize translational joint ground - sled
    m_prismatic = chrono_types::make_shared<ChLinkLockPrismatic>();
    m_prismatic->Initialize(m_ground, m_sled, ChCoordsys<>(loc_prismatic, Q_from_AngY(CH_C_PI_2)));
    system->AddLink(m_prismatic);

    // Create and initialize revolute joint sled - wheel
    m_revolute = chrono_types::make_shared<ChLinkLockRevolute>();
    m_revolute->Initialize(m_wheel, m_sled, ChCoordsys<>(loc_revolute, Q_from_AngX(CH_C_PI_2)));
    system->AddLink(m_revolute);
}

void Mechanism::WriteResults(ChStreamOutAsciiFile& f, double time) {
    // Velocity of sled body (in absolute frame)
    ChVector<> sled_vel = m_sled->GetPos_dt();

    // Velocity of wheel body (in absolute frame)
    ChVector<> wheel_vel = m_wheel->GetPos_dt();

    // Coordinate system of the revolute joint (relative to the frame of body2, in
    // this case the sled)
    ChCoordsys<> revCoordsys = m_revolute->GetLinkRelativeCoords();

    // Reaction force in revolute joint
    ChVector<> force_jointsys = m_revolute->Get_react_force();
    ChVector<> force_bodysys = revCoordsys.TransformDirectionLocalToParent(force_jointsys);
    ChVector<> force_abssys = m_sled->GetCoord().TransformDirectionLocalToParent(force_bodysys);

    f << time << "  " << sled_vel.x() << "  " << sled_vel.y() << "  " << sled_vel.z() << "      " << wheel_vel.x() << "  "
      << wheel_vel.y() << "  " << wheel_vel.z() << "      " << force_abssys.x() << "  " << force_abssys.y() << "  "
      << force_abssys.z() << "\n";
}

// =============================================================================
// Create container bin.
// =============================================================================
void CreateContainer(ChSystemParallel* system) {
    int id_c = -200;
    double thickness = 0.2;

#ifdef USE_SMC
    auto mat_c = chrono_types::make_shared<ChMaterialSurfaceSMC>();
    mat_c->SetYoungModulus(2e6f);
    mat_c->SetFriction(0.4f);
    mat_c->SetRestitution(0.1f);

    utils::CreateBoxContainer(system, id_c, mat_c, ChVector<>(L / 2, W / 2, H / 2), thickness / 2);
#else
    auto mat_c = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    mat_c->SetFriction(0.4f);

    utils::CreateBoxContainer(system, id_c, mat_c, ChVector<>(L / 2, W / 2, H / 2), thickness / 2);

#endif
}

// =============================================================================
// Create granular material.
// =============================================================================
void CreateParticles(ChSystemParallel* system) {
// Create a material for the ball mixture.
#ifdef USE_SMC
    auto mat_g = chrono_types::make_shared<ChMaterialSurfaceSMC>();
    mat_g->SetYoungModulus(Y_g);
    mat_g->SetFriction(mu_g);
    mat_g->SetRestitution(cr_g);
#else
    auto mat_g = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    mat_g->SetFriction(mu_g);
#endif

    // Create a mixture entirely made out of spheres.
    utils::Generator gen(system);

    std::shared_ptr<utils::MixtureIngredient> m1 = gen.AddMixtureIngredient(utils::MixtureType::SPHERE, 1.0);
    m1->setDefaultMaterial(mat_g);
    m1->setDefaultDensity(rho_g);
    m1->setDefaultSize(r_g);

    // Create particles, one layer at a time, until the desired number is reached.
    gen.setBodyIdentifier(100);

    double r = 1.01 * r_g;
    ChVector<> hdims(L / 2 - r, W / 2 - r, 0);
    ChVector<> center(0, 0, 2 * r);

    int layer = 1;
    while (gen.getTotalNumBodies() < desired_num_particles) {
        gen.createObjectsBox(utils::SamplingType::POISSON_DISK, 2 * r, center, hdims);
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
    for (auto body : sys->Get_bodylist()) {
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
    ChSystemParallelSMC* msystem = new ChSystemParallelSMC();
#else
    cout << "Create NSC system" << endl;
    ChSystemParallelNSC* msystem = new ChSystemParallelNSC();
#endif

    msystem->Set_G_acc(ChVector<>(0, 0, -9.81));

    // ----------------------
    // Set number of threads.
    // ----------------------

    int max_threads = omp_get_num_procs();
    if (threads > max_threads)
        threads = max_threads;
    omp_set_num_threads(threads);
    cout << "Using " << threads << " threads" << endl;

    msystem->GetSettings()->max_threads = threads;
    msystem->GetSettings()->perform_thread_tuning = thread_tuning;

    // ---------------------
    // Edit system settings.
    // ---------------------

    msystem->GetSettings()->solver.tolerance = tolerance;
    msystem->GetSettings()->solver.max_iteration_bilateral = max_iteration_bilateral;
    msystem->GetSettings()->solver.clamp_bilaterals = clamp_bilaterals;
    msystem->GetSettings()->solver.bilateral_clamp_speed = bilateral_clamp_speed;

#ifdef USE_SMC
    msystem->GetSettings()->collision.narrowphase_algorithm = NarrowPhaseType::NARROWPHASE_R;
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

    msystem->GetSettings()->collision.bins_per_axis = vec3(10, 10, 10);

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
            CreateContainer(msystem);
            CreateParticles(msystem);

            break;

        case PUSHING: {
            time_end = time_pushing;
            out_fps = out_fps_pushing;

            // Create the granular material and the container from the checkpoint file.
            cout << "Read checkpoint data from " << checkpoint_file;
            utils::ReadCheckpoint(msystem, checkpoint_file);
            cout << "  done.  Read " << msystem->Get_bodylist().size() << " bodies." << endl;

            // Create the mechanism with the wheel just above the granular material.
            double lowest, highest;
            FindRange(msystem, lowest, highest);
            cout << "Create mechanism above height " << highest + r_g << endl;
            mech = new Mechanism(msystem, highest + r_g);
        }

        break;

        case TESTING:
            time_end = time_pushing;
            out_fps = out_fps_pushing;

            mech = new Mechanism(msystem, 0.9 * H);

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
    ChStreamOutAsciiFile sfile(stats_file.c_str());
    ChStreamOutAsciiFile rfile(results_file.c_str());

    // Disable buffering on output file streams.
    sfile.GetFstream().rdbuf()->pubsetbuf(0, 0);
    rfile.GetFstream().rdbuf()->pubsetbuf(0, 0);

#ifdef CHRONO_OPENGL
    opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
    gl_window.Initialize(1280, 720, "Pressure Sinkage Test", msystem);
    gl_window.SetCamera(ChVector<>(0, -8, 0), ChVector<>(0, 0, 0), ChVector<>(0, 0, 1));
    gl_window.SetRenderMode(opengl::WIREFRAME);
#endif

    while (time < time_end) {
        if (sim_frame == next_out_frame) {
            char filename[100];
            sprintf(filename, "%s/data_%03d.dat", pov_dir.c_str(), out_frame + 1);
            utils::WriteShapesPovray(msystem, filename, problem == TESTING);
            double highest, lowest;
            FindRange(msystem, lowest, highest);
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
            sfile.GetFstream().flush();

            // Create a checkpoint from the current state.
            if (problem == SETTLING) {
                cout << "             Write checkpoint data " << flush;
                utils::WriteCheckpoint(msystem, checkpoint_file);
                cout << msystem->Get_bodylist().size() << " bodies" << endl;
            }

            out_frame++;
            next_out_frame += out_steps;
            num_contacts = 0;
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

        // Increment counters
        time += time_step;
        sim_frame++;
        exec_time += msystem->GetTimerStep();
        num_contacts += msystem->GetNcontacts();

        // Save results
        if (problem != SETTLING) {
            mech->WriteResults(rfile, time);
            rfile.GetFstream().flush();
        }
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
    cout << "Simulation time:   " << exec_time << endl;
    cout << "Number of threads: " << threads << endl;

    return 0;
}

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
// Chrono::Vehicle + ChronoParallel demo program for simulating a HMMWV vehicle
// traversing a gravel ditch.
//
// The vehicle model uses the utility class ChWheeledVehicleAssembly and is
// based on JSON specification files from the Chrono data directory.
//
// Contact uses the DEM-P (penalty) formulation.
//
// The global reference frame has Z up.
// All units SI.
// =============================================================================

#include <stdio.h>
#include <vector>
#include <cmath>

#include "chrono/ChConfig.h"
#include "chrono/core/ChFileutils.h"
#include "chrono/core/ChStream.h"
#include "chrono/utils/ChUtilsGeometry.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsGenerators.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#include "chrono_vehicle/ChVehicleModelData.h"
#include "chrono_vehicle/wheeled_vehicle/utils/ChWheeledVehicleAssembly.h"

#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_parallel/lcp/ChLcpSystemDescriptorParallel.h"
#include "chrono_parallel/collision/ChCNarrowphaseRUtils.h"

// Control use of OpenGL run-time rendering
// Note: CHRONO_OPENGL is defined in ChConfig.h
//#undef CHRONO_OPENGL

#ifdef CHRONO_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif

#include "utils/demo_utils.h"

using namespace chrono;
using namespace chrono::collision;
using namespace chrono::vehicle;

using std::cout;
using std::endl;

// =============================================================================

// -----------------------------------------------------------------------------
// Specification of the terrain
// -----------------------------------------------------------------------------

// Control visibility of containing bin walls
bool visible_walls = false;

// Dimensions
double hdimX = 2.5;
double hdimY = 1.75;
double hdimZ = 0.5;
double hthick = 0.25;
double hlen = 2.5;

// Parameters for granular material
int Id_g = 1;
double r_g = 0.02;
double rho_g = 2500;
double vol_g = (4.0 / 3) * CH_C_PI * r_g * r_g * r_g;
double mass_g = rho_g * vol_g;
ChVector<> inertia_g = 0.4 * mass_g * r_g * r_g * ChVector<>(1, 1, 1);

float cohesion_g = 200;

float Y_g = 1e8;
float cr_g = 0.1;
float mu_g = 0.8;

int coll_fam_g = 1;

int num_particles = 100000;

// -----------------------------------------------------------------------------
// Specification of the vehicle model
// -----------------------------------------------------------------------------

ChWheeledVehicleAssembly* vehicle_assembly = NULL;
ChDriverInputsCallback* driver_cb = NULL;
ChTireContactCallback* tire_cb = NULL;

enum WheelType { CYLINDRICAL, LUGGED };

// Type of wheel/tire (controls both contact and visualization)
WheelType wheel_type = CYLINDRICAL;

// JSON files for vehicle model (using different wheel visualization meshes)
std::string vehicle_file_cyl("hmmwv/vehicle/HMMWV_Vehicle_simple.json");
std::string vehicle_file_lug("hmmwv/vehicle/HMMWV_Vehicle_simple_lugged.json");

// JSON files for powertrain (simple)
std::string simplepowertrain_file("hmmwv/powertrain/HMMWV_SimplePowertrain.json");

// Initial vehicle position and orientation
ChVector<> initLoc(-hdimX - hlen, 0, 0.6);
ChQuaternion<> initRot(1, 0, 0, 0);

// Contact material properties for rigid tires
float Y_t = 1e8;
float cr_t = 0.1;
float mu_t = 0.8;

int coll_fam_t = 2;

// -----------------------------------------------------------------------------
// Simulation parameters
// -----------------------------------------------------------------------------

// Desired number of OpenMP threads (will be clamped to maximum available)
int threads = 20;

// Perform dynamic tuning of number of threads?
bool thread_tuning = false;

// Total simulation duration.
double time_end = 7;

// Duration of the "hold time" (vehicle chassis fixed and no driver inputs).
// This can be used to allow the granular material to settle.
double time_hold = 0.3;

// Solver parameters
double time_step = 2e-5;

double tolerance = 0.1;

int max_iteration_bilateral = 1000;

// Contact force model
ChSystemDEM::ContactForceModel contact_force_model = ChSystemDEM::ContactForceModel::Hooke;
ChSystemDEM::TangentialDisplacementModel tangential_displ_mode = ChSystemDEM::TangentialDisplacementModel::OneStep;

// Periodically monitor maximum bilateral constraint violation
bool monitor_bilaterals = false;
int bilateral_frame_interval = 100;

// Output
bool povray_output = false;

const std::string out_dir = "../HMMWV_DEM_DITCH";
const std::string pov_dir = out_dir + "/POVRAY";

int out_fps = 60;

// =============================================================================

// Callback class for providing driver inputs.
class MyDriverInputs : public ChDriverInputsCallback {
  public:
    MyDriverInputs(double delay) : m_delay(delay) {}

    virtual void onCallback(double time, double& throttle, double& steering, double& braking) {
        throttle = 0;
        steering = 0;
        braking = 0;

        double eff_time = time - m_delay;

        // Do not generate any driver inputs for a duration equal to m_delay.
        if (eff_time < 0)
            return;

        if (eff_time > 0.2)
            throttle = 1.0;
        else if (eff_time > 0.1)
            throttle = 10 * (eff_time - 0.1);

        ////cout << "     t = " << time << "  eff_t = " << eff_time << " Throttle = " << throttle << endl;
    }

  private:
    double m_delay;
};

// Callback class for specifying rigid tire contact model.
// This version uses cylindrical contact shapes.
class MyCylindricalTire : public ChTireContactCallback {
  public:
    virtual void onCallback(std::shared_ptr<ChBody> wheelBody, double radius, double width) {
        wheelBody->ChangeCollisionModel(new collision::ChCollisionModelParallel);

        wheelBody->GetCollisionModel()->ClearModel();
        wheelBody->GetCollisionModel()->AddCylinder(0.46, 0.46, width / 2);
        wheelBody->GetCollisionModel()->BuildModel();

        wheelBody->GetCollisionModel()->SetFamily(coll_fam_t);

        wheelBody->GetMaterialSurfaceDEM()->SetFriction(mu_t);
        wheelBody->GetMaterialSurfaceDEM()->SetYoungModulus(Y_t);
        wheelBody->GetMaterialSurfaceDEM()->SetRestitution(cr_t);
    }
};

// Callback class for specifying rigid tire contact model.
// This version uses a collection of convex contact shapes (meshes).
class MyLuggedTire : public ChTireContactCallback {
  public:
    MyLuggedTire() {
        std::string lugged_file("hmmwv/lugged_wheel_section.obj");
        geometry::ChTriangleMeshConnected lugged_mesh;
        utils::LoadConvexMesh(vehicle::GetDataFile(lugged_file), lugged_mesh, lugged_convex);
        num_hulls = lugged_convex.GetHullCount();
    }

    virtual void onCallback(std::shared_ptr<ChBody> wheelBody, double radius, double width) {
        wheelBody->ChangeCollisionModel(new collision::ChCollisionModelParallel);

        ChCollisionModelParallel* coll_model = (ChCollisionModelParallel*)wheelBody->GetCollisionModel();
        coll_model->ClearModel();

        // Assemble the tire contact from 15 segments, properly offset.
        // Each segment is further decomposed in convex hulls.
        for (int iseg = 0; iseg < 15; iseg++) {
            ChQuaternion<> rot = Q_from_AngAxis(iseg * 24 * CH_C_DEG_TO_RAD, VECT_Y);
            for (int ihull = 0; ihull < num_hulls; ihull++) {
                std::vector<ChVector<> > convexhull;
                lugged_convex.GetConvexHullResult(ihull, convexhull);
                coll_model->AddConvexHull(convexhull, VNULL, rot);
            }
        }

        // Add a cylinder to represent the wheel hub.
        coll_model->AddCylinder(0.223, 0.223, 0.126);

        coll_model->BuildModel();

        coll_model->SetFamily(coll_fam_t);

        wheelBody->GetMaterialSurfaceDEM()->SetFriction(mu_t);
        wheelBody->GetMaterialSurfaceDEM()->SetYoungModulus(Y_t);
        wheelBody->GetMaterialSurfaceDEM()->SetRestitution(cr_t);
    }

  private:
    ChConvexDecompositionHACDv2 lugged_convex;
    int num_hulls;
};

// =============================================================================

void CreateGroundGeometry(std::shared_ptr<ChBody> ground) {
    ground->GetCollisionModel()->ClearModel();

    // Bottom box
    utils::AddBoxGeometry(ground.get(), ChVector<>(hdimX, hdimY, hthick), ChVector<>(0, 0, -hthick),
                          ChQuaternion<>(1, 0, 0, 0), true);
    // Left box
    utils::AddBoxGeometry(ground.get(), ChVector<>(hdimX, hthick, hdimZ + hthick),
                          ChVector<>(0, hdimY + hthick, hdimZ - hthick), ChQuaternion<>(1, 0, 0, 0), visible_walls);
    // Right box
    utils::AddBoxGeometry(ground.get(), ChVector<>(hdimX, hthick, hdimZ + hthick),
                          ChVector<>(0, -hdimY - hthick, hdimZ - hthick), ChQuaternion<>(1, 0, 0, 0), visible_walls);

    // Front box
    utils::AddBoxGeometry(ground.get(), ChVector<>(hthick, hdimY, hdimZ + hthick),
                          ChVector<>(hdimX + hthick, 0, hdimZ - hthick), ChQuaternion<>(1, 0, 0, 0), visible_walls);
    // Rear box
    utils::AddBoxGeometry(ground.get(), ChVector<>(hthick, hdimY, hdimZ + hthick),
                          ChVector<>(-hdimX - hthick, 0, hdimZ - hthick), ChQuaternion<>(1, 0, 0, 0), visible_walls);

    ground->GetCollisionModel()->BuildModel();

    ground->GetCollisionModel()->SetFamily(coll_fam_g);
    ground->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(coll_fam_t);
}

void AdjustGroundGeometry(std::shared_ptr<ChBody> ground, double platform_height) {
    ground->GetAssets().clear();
    ground->GetCollisionModel()->ClearModel();

    // Bottom box
    utils::AddBoxGeometry(ground.get(), ChVector<>(hdimX, hdimY, hthick), ChVector<>(0, 0, -hthick),
                          ChQuaternion<>(1, 0, 0, 0), true);
    // Left box
    utils::AddBoxGeometry(ground.get(), ChVector<>(hdimX, hthick, hdimZ + hthick),
                          ChVector<>(0, hdimY + hthick, hdimZ - hthick), ChQuaternion<>(1, 0, 0, 0), visible_walls);
    // Right box
    utils::AddBoxGeometry(ground.get(), ChVector<>(hdimX, hthick, hdimZ + hthick),
                          ChVector<>(0, -hdimY - hthick, hdimZ - hthick), ChQuaternion<>(1, 0, 0, 0), visible_walls);

    // Front platform
    utils::AddBoxGeometry(ground.get(), ChVector<>(hlen, hdimY, platform_height / 2 + hthick),
                          ChVector<>(hdimX + hlen, 0, platform_height / 2 - hthick), ChQuaternion<>(1, 0, 0, 0), true);
    // Rear platform
    utils::AddBoxGeometry(ground.get(), ChVector<>(hlen, hdimY, platform_height / 2 + hthick),
                          ChVector<>(-hdimX - hlen, 0, platform_height / 2 - hthick), ChQuaternion<>(1, 0, 0, 0), true);

    ground->GetCollisionModel()->BuildModel();
}

void CreatePlatform(ChSystem* system, double platform_height) {
    auto platform = std::make_shared<ChBody>(new collision::ChCollisionModelParallel, ChMaterialSurfaceBase::DEM);
    platform->SetIdentifier(-2);
    platform->SetMass(1000);
    platform->SetBodyFixed(true);
    platform->SetCollide(true);

    platform->GetMaterialSurfaceDEM()->SetFriction(mu_g);
    platform->GetMaterialSurfaceDEM()->SetYoungModulus(Y_g);
    platform->GetMaterialSurfaceDEM()->SetRestitution(cr_g);
    platform->GetMaterialSurfaceDEM()->SetAdhesion(cohesion_g);

    platform->GetCollisionModel()->ClearModel();

    // Front platform
    utils::AddBoxGeometry(platform.get(), ChVector<>(hlen, hdimY, platform_height / 2 + hthick),
                          ChVector<>(hdimX + hlen, 0, platform_height / 2 - hthick), ChQuaternion<>(1, 0, 0, 0), true);
    // Rear platform
    utils::AddBoxGeometry(platform.get(), ChVector<>(hlen, hdimY, platform_height / 2 + hthick),
                          ChVector<>(-hdimX - hlen, 0, platform_height / 2 - hthick), ChQuaternion<>(1, 0, 0, 0), true);

    platform->GetCollisionModel()->BuildModel();

    system->AddBody(platform);
}

// =============================================================================

void CreateVehicleAssembly(ChSystem* system, double vertical_offset) {
    // Create the vehicle assembly and the callback object for tire contact
    // according to the specified type of tire/wheel.
    switch (wheel_type) {
        case CYLINDRICAL: {
            vehicle_assembly = new ChWheeledVehicleAssembly(system, vehicle_file_cyl, simplepowertrain_file);
            tire_cb = new MyCylindricalTire();
        } break;
        case LUGGED: {
            vehicle_assembly = new ChWheeledVehicleAssembly(system, vehicle_file_lug, simplepowertrain_file);
            tire_cb = new MyLuggedTire();
        } break;
    }

    vehicle_assembly->SetTireContactCallback(tire_cb);

    // Set the callback object for driver inputs. Pass the hold time as a delay in
    // generating driver inputs.
    driver_cb = new MyDriverInputs(time_hold);
    vehicle_assembly->SetDriverInputsCallback(driver_cb);

    // Initialize the vehicle at the specified location.
    vehicle_assembly->Initialize(initLoc + ChVector<>(0, 0, vertical_offset), initRot);
}

// =============================================================================

int CreateParticles(ChSystem* system) {
    // Create a material
    auto mat_g = std::make_shared<ChMaterialSurfaceDEM>();
    mat_g->SetYoungModulus(Y_g);
    mat_g->SetFriction(mu_g);
    mat_g->SetRestitution(cr_g);
    mat_g->SetAdhesion(cohesion_g);

    // Create a particle generator and a mixture entirely made out of spheres
    utils::Generator gen(system);
    std::shared_ptr<utils::MixtureIngredient>& m1 = gen.AddMixtureIngredient(utils::SPHERE, 1.0);
    m1->setDefaultMaterialDEM(mat_g);
    m1->setDefaultDensity(rho_g);
    m1->setDefaultSize(r_g);

    // Set starting value for body identifiers
    gen.setBodyIdentifier(Id_g);

    // Create particles in layers until reaching the desired number of particles
    double r = 1.01 * r_g;
    ChVector<> hdims(hdimX - r, hdimY - r, 0);
    ChVector<> center(0, 0, 2 * r);

    while (gen.getTotalNumBodies() < num_particles) {
        gen.createObjectsBox(utils::POISSON_DISK, 2 * r, center, hdims);
        center.z += 2 * r;
    }

    return gen.getTotalNumBodies();
}

double FindHighestParticle(ChSystem* system) {
    double highest = 0;
    for (size_t i = 0; i < system->Get_bodylist()->size(); ++i) {
        auto body = (*system->Get_bodylist())[i];
        if (body->GetIdentifier() > 0 && body->GetPos().z > highest)
            highest = body->GetPos().z;
    }
    return highest;
}

// =============================================================================

int main(int argc, char* argv[]) {
    // Set path to Chrono and Chrono::Vehicle data directories
    SetChronoDataPath(CHRONO_DATA_DIR);
    vehicle::SetDataPath(CHRONO_VEHICLE_DATA_DIR);

    // --------------------------
    // Create output directories.
    // --------------------------

    if (povray_output) {
        if (ChFileutils::MakeDirectory(out_dir.c_str()) < 0) {
            cout << "Error creating directory " << out_dir << endl;
            return 1;
        }
        if (ChFileutils::MakeDirectory(pov_dir.c_str()) < 0) {
            cout << "Error creating directory " << pov_dir << endl;
            return 1;
        }
    }

    // --------------
    // Create system.
    // --------------

    ChSystemParallelDEM* system = new ChSystemParallelDEM();

    system->Set_G_acc(ChVector<>(0, 0, -9.81));

    // ----------------------
    // Enable debug log
    // ----------------------

    ////system->SetLoggingLevel(LOG_INFO, true);
    ////system->SetLoggingLevel(LOG_TRACE, true);

    // ----------------------
    // Set number of threads.
    // ----------------------
    int max_threads = omp_get_num_procs();
    if (threads > max_threads)
        threads = max_threads;
    system->SetParallelThreadNumber(threads);
    omp_set_num_threads(threads);
    cout << "Using " << threads << " threads" << endl;

    system->GetSettings()->perform_thread_tuning = thread_tuning;

    // ---------------------
    // Edit system settings.
    // ---------------------

    system->GetSettings()->solver.use_full_inertia_tensor = false;

    system->GetSettings()->solver.tolerance = tolerance;
    system->GetSettings()->solver.max_iteration_bilateral = max_iteration_bilateral;

    system->GetSettings()->solver.adhesion_force_model = ChSystemDEM::AdhesionForceModel::Constant;
    system->GetSettings()->solver.contact_force_model = contact_force_model;
    system->GetSettings()->solver.tangential_displ_mode = tangential_displ_mode;

    system->GetSettings()->collision.narrowphase_algorithm = NARROWPHASE_HYBRID_MPR;
    system->GetSettings()->collision.bins_per_axis = I3(20, 20, 10);

    // -------------------
    // Create the terrain.
    // -------------------

    // Ground body
    auto ground = std::make_shared<ChBody>(new collision::ChCollisionModelParallel, ChMaterialSurfaceBase::DEM);
    ground->SetIdentifier(-1);
    ground->SetMass(1000);
    ground->SetBodyFixed(true);
    ground->SetCollide(true);

    ground->GetMaterialSurfaceDEM()->SetFriction(mu_g);
    ground->GetMaterialSurfaceDEM()->SetYoungModulus(Y_g);
    ground->GetMaterialSurfaceDEM()->SetRestitution(cr_g);
    ground->GetMaterialSurfaceDEM()->SetAdhesion(cohesion_g);

    CreateGroundGeometry(ground);

    system->AddBody(ground);

    // Create the granular material.
    int num_particles = CreateParticles(system);
    cout << "Created " << num_particles << " particles." << endl;

    // -------------------
    // Specify active box.
    // -------------------

    system->GetSettings()->collision.use_aabb_active = true;
    system->GetSettings()->collision.aabb_min = R3(-hdimX - 2 * hlen, -hdimY, 0);
    system->GetSettings()->collision.aabb_max = R3(hdimX + 2 * hlen, hdimY, 2 * hdimZ);

// -----------------------
// Start the simulation.
// -----------------------

#ifdef CHRONO_OPENGL
    // Initialize OpenGL
    opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
    gl_window.Initialize(1280, 720, "HMMWV ditch (DEM contact)", system);
    gl_window.SetCamera(ChVector<>(0, -10, 0), ChVector<>(0, 0, 0), ChVector<>(0, 0, 1));
    gl_window.SetRenderMode(opengl::WIREFRAME);
#endif

    // Run simulation for specified time.
    int out_steps = std::ceil((1.0 / time_step) / out_fps);

    double time = 0;
    int sim_frame = 0;
    int out_frame = 0;
    int next_out_frame = 0;
    double exec_time = 0;
    int num_contacts = 0;

    while (time < time_end) {
        // If enabled, output data for PovRay postprocessing.
        if (sim_frame == next_out_frame) {
            double speed = vehicle_assembly ? vehicle_assembly->GetVehicle()->GetVehicleSpeed() : 0;
            cout << endl;
            cout << "---- Frame:          " << out_frame + 1 << endl;
            cout << "     Sim frame:      " << sim_frame << endl;
            cout << "     Time:           " << time << endl;
            cout << "     Speed:          " << speed << endl;
            cout << "     Avg. contacts:  " << num_contacts / out_steps << endl;
            cout << "     Execution time: " << exec_time << endl;

            if (povray_output) {
                char filename[100];
                sprintf(filename, "%s/data_%03d.dat", pov_dir.c_str(), out_frame + 1);
                utils::WriteShapesPovray(system, filename);
            }

            out_frame++;
            next_out_frame += out_steps;
            num_contacts = 0;
        }

        // At the end of the hold time, adjust ground geometry and create the vehicle.
        if (!vehicle_assembly && time > time_hold) {
            cout << endl << "Create vehicle at t = " << time << endl;
            double max_height = FindHighestParticle(system);
            ////AdjustGroundGeometry(ground, max_height);
            CreatePlatform(system, max_height);
            CreateVehicleAssembly(system, max_height);
        }

        // Update vehicle
        if (vehicle_assembly)
            vehicle_assembly->Synchronize(time);

// Advance dynamics.
#ifdef CHRONO_OPENGL
        if (gl_window.Active()) {
            gl_window.DoStepDynamics(time_step);
            gl_window.Render();
        } else
            break;
#else
        system->DoStepDynamics(time_step);
#endif

        ////progressbar(out_steps + sim_frame - next_out_frame + 1, out_steps);
        TimingOutput(system);

        // Periodically display maximum constraint violation
        if (monitor_bilaterals && sim_frame % bilateral_frame_interval == 0) {
            std::vector<double> cvec;
            cout << "  Max. violation = " << system->CalculateConstraintViolation(cvec) << endl;
        }

        // Update counters.
        time += time_step;
        sim_frame++;
        exec_time += system->GetTimerStep();
        num_contacts += system->GetNcontacts();
    }

    // Final stats
    cout << "==================================" << endl;
    cout << "Simulation time:   " << exec_time << endl;
    cout << "Number of threads: " << threads << endl;

    delete vehicle_assembly;
    delete driver_cb;
    delete tire_cb;

    return 0;
}

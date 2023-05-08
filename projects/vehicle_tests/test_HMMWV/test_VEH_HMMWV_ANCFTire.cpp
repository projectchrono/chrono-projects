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
// HMMWV full model using ANCF, RIGID, or RIGID_MESH tires on rigid terrain.
//
// The vehicle reference frame has Z up, X towards the front of the vehicle, and
// Y pointing to the left.
//
// =============================================================================

////#include <float.h>
////unsigned int fp_control_state = _controlfp(_EM_INEXACT, _MCW_EM);

#include "chrono/core/ChStream.h"
#include "chrono/physics/ChSystem.h"
#include "chrono/utils/ChOpenMP.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#ifdef CHRONO_PARDISO_MKL
#include "chrono_pardisomkl/ChSolverPardisoMKL.h"
#endif

#include "chrono_vehicle/ChConfigVehicle.h"
#include "chrono_vehicle/ChVehicleModelData.h"
#include "chrono_vehicle/driver/ChDataDriver.h"
#include "chrono_vehicle/driver/ChInteractiveDriverIRR.h"
#include "chrono_vehicle/terrain/RigidTerrain.h"
#include "chrono_vehicle/wheeled_vehicle/ChWheeledVehicleVisualSystemIrrlicht.h"

#include "chrono_models/vehicle/hmmwv/HMMWV.h"
#include "chrono_models/vehicle/hmmwv/tire/HMMWV_ANCFTire.h"
#include "chrono_models/vehicle/hmmwv/tire/HMMWV_RigidTire.h"

#include "chrono_thirdparty/filesystem/path.h"

using namespace chrono;
using namespace chrono::vehicle;
using namespace chrono::vehicle::hmmwv;

// =============================================================================

// Number of OpenMP threads
int num_threads = 4;

// Initial vehicle location and orientation
ChVector<> initLoc(0, 0, 1.2);
ChQuaternion<> initRot(1, 0, 0, 0);
////ChQuaternion<> initRot(0.866025, 0, 0, 0.5);
////ChQuaternion<> initRot(0.7071068, 0, 0, 0.7071068);
////ChQuaternion<> initRot(0.25882, 0, 0, 0.965926);
////ChQuaternion<> initRot(0, 0, 0, 1);

// Visualization type for chassis (PRIMITIVES, MESH, or NONE)
VisualizationType chassis_vis_type = VisualizationType::PRIMITIVES;

// Type of tire type (ANCF, RIGID, RIGID_MESH)
TireModelType tire_model = TireModelType::ANCF;

// Type of powertrain models (SHAFTS, SIMPLE)
EngineModelType engine_model = EngineModelType::SHAFTS;
TransmissionModelType transmission_model = TransmissionModelType::SHAFTS;

// Drive type (FWD, RWD, or AWD)
DrivelineTypeWV drive_type = DrivelineTypeWV::AWD;

// Rigid terrain (RigidTerrain::PatchType::FLAT, RigidTerrain::PatchType::HEIGHT_MAP, RigidTerrain::PatchType::MESH)
RigidTerrain::PatchType terrain_model = RigidTerrain::PatchType::BOX;

// Use material properties for SMC contact method?
bool use_mat_properties = true;

// Terrain dimensions (for FLAT terrain)
double terrainHeight = 0;      // terrain height (FLAT terrain only)
double terrainLength = 100.0;  // size in X direction
double terrainWidth = 100.0;   // size in Y direction

// Point on chassis tracked by the camera (for Irrlicht visualization)
ChVector<> trackPoint(0.0, 0.0, 1.75);

// Simulation step sizes
double step_size = 1e-4;
// Simulation end time
double t_end = 5;
// Verbose solver output
bool verbose = false;

// Time interval between two render frames (1/FPS)
double render_step_size = 1.0 / 50;

// Output directories
const std::string out_dir = "../HMMWV_ANCF";
const std::string pov_dir = out_dir + "/POVRAY";

// Debug logging
bool debug_output = false;
// Debug output frequency (1/FPS)
double debug_step_size = 1.0 / 2;

// POV-Ray output
bool povray_output = false;

// =============================================================================

class MyDriver : public ChDriver {
  public:
    MyDriver(ChVehicle& vehicle, double delay) : ChDriver(vehicle), m_delay(delay) {}
    ~MyDriver() {}

    virtual void Synchronize(double time) override {
        m_throttle = 0;
        m_steering = 0;
        m_braking = 0;

        double eff_time = time - m_delay;

        // Do not generate any driver inputs for a duration equal to m_delay.
        if (eff_time < 0)
            return;

        if (eff_time > 0.2)
            m_throttle = 0.8;
        else
            m_throttle = 4 * eff_time;
    }

  private:
    double m_delay;
};

// =============================================================================

int main(int argc, char* argv[]) {
    // Set path to Chrono and Chrono::Vehicle data directories
    SetChronoDataPath(CHRONO_DATA_DIR);
    vehicle::SetDataPath(CHRONO_VEHICLE_DATA_DIR);

    // ----------------------------------
    // Create the (sequential) SMC system
    // ----------------------------------

    ChSystemSMC* system = new ChSystemSMC(use_mat_properties);
    system->Set_G_acc(ChVector<>(0, 0, -9.81));

    // Set number threads
    system->SetNumThreads(num_threads);

#ifdef CHRONO_PARDISO_MKL
    // PardisoMKL solver settings
    auto mkl_solver = chrono_types::make_shared<ChSolverPardisoMKL>();
    mkl_solver->LockSparsityPattern(true);
    mkl_solver->SetVerbose(verbose);
    system->SetSolver(mkl_solver);
#else
    // Default solver settings
    system->SetSolverType(ChSolver::Type::PSOR);
    system->SetSolverMaxIterations(100);
    system->SetSolverTolerance(1e-10);
    system->SetSolverForceTolerance(1e-8);
#endif

    // Integrator settings
    system->SetTimestepperType(ChTimestepper::Type::HHT);
    auto integrator = std::static_pointer_cast<ChTimestepperHHT>(system->GetTimestepper());
    integrator->SetAlpha(-0.2);
    integrator->SetMaxiters(50);
    integrator->SetAbsTolerances(5e-05, 1.8);
    integrator->SetMode(ChTimestepperHHT::POSITION);
    integrator->SetStepControl(true);
    integrator->SetModifiedNewton(false);
    integrator->SetScaling(true);
    integrator->SetVerbose(verbose);
    integrator->SetMaxItersSuccess(5);

    // --------------
    // Create systems
    // --------------

    // Create the HMMWV vehicle, set parameters, and initialize
    HMMWV_Full my_hmmwv(system);
    my_hmmwv.SetChassisFixed(false);
    my_hmmwv.SetInitPosition(ChCoordsys<>(initLoc, initRot));
    my_hmmwv.SetEngineType(engine_model);
    my_hmmwv.SetTransmissionType(transmission_model);
    my_hmmwv.SetDriveType(drive_type);
    my_hmmwv.SetTireType(tire_model);
    my_hmmwv.Initialize();

    my_hmmwv.SetChassisVisualizationType(chassis_vis_type);
    my_hmmwv.SetSuspensionVisualizationType(VisualizationType::PRIMITIVES);
    my_hmmwv.SetSteeringVisualizationType(VisualizationType::PRIMITIVES);
    my_hmmwv.SetWheelVisualizationType(VisualizationType::NONE);
    my_hmmwv.SetTireVisualizationType(VisualizationType::MESH);

    // Downcast tires (if needed)
    switch (tire_model) {
        case TireModelType::ANCF: {
            my_hmmwv.GetVehicle().GetWheel(0, LEFT)->GetTire();
            auto tire_FL = std::static_pointer_cast<HMMWV_ANCFTire>(my_hmmwv.GetVehicle().GetWheel(0, LEFT)->GetTire());
            auto tire_FR = std::static_pointer_cast<HMMWV_ANCFTire>(my_hmmwv.GetVehicle().GetWheel(0, LEFT)->GetTire());
            auto tire_RL = std::static_pointer_cast<HMMWV_ANCFTire>(my_hmmwv.GetVehicle().GetWheel(0, LEFT)->GetTire());
            auto tire_RR = std::static_pointer_cast<HMMWV_ANCFTire>(my_hmmwv.GetVehicle().GetWheel(0, LEFT)->GetTire());
            ////tire_FL->EnablePressure(false);
            ////tire_FR->EnablePressure(false);
            ////tire_RL->EnablePressure(false);
            ////tire_RR->EnablePressure(false);

            break;
        }
    }

    // Create the terrain
    auto patch_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    patch_mat->SetFriction(0.9f);
    patch_mat->SetRestitution(0.01f);

    RigidTerrain terrain(my_hmmwv.GetSystem());
    std::shared_ptr<RigidTerrain::Patch> patch;
    switch (terrain_model) {
        case RigidTerrain::PatchType::BOX:
            patch = terrain.AddPatch(patch_mat, ChCoordsys<>(ChVector<>(0, 0, terrainHeight), QUNIT), terrainLength,
                                     terrainWidth);
            patch->SetTexture(vehicle::GetDataFile("terrain/textures/tile4.jpg"), 200, 200);
            break;
        case RigidTerrain::PatchType::HEIGHT_MAP:
            patch = terrain.AddPatch(patch_mat, CSYSNORM, vehicle::GetDataFile("terrain/height_maps/test64.bmp"), 128,
                                     128, 0, 4);
            patch->SetTexture(vehicle::GetDataFile("terrain/textures/grass.jpg"), 16, 16);
            break;
        case RigidTerrain::PatchType::MESH:
            patch = terrain.AddPatch(patch_mat, CSYSNORM, vehicle::GetDataFile("terrain/meshes/test.obj"));
            patch->SetTexture(vehicle::GetDataFile("terrain/textures/grass.jpg"), 100, 100);
            break;
    }
    patch->SetColor(ChColor(0.8f, 0.8f, 0.5f));
    terrain.Initialize();

    // Create the vehicle Irrlicht interface
    auto vis = chrono_types::make_shared<ChWheeledVehicleVisualSystemIrrlicht>();
    vis->SetWindowTitle("HMMWV ANCF tires Test");
    vis->SetChaseCamera(trackPoint, 6.0, 0.5);
    vis->Initialize();
    vis->AddTypicalLights();
    vis->AddSkyBox();
    vis->AddLogo();
    vis->AttachVehicle(&my_hmmwv.GetVehicle());

    // -----------------
    // Initialize output
    // -----------------

    if (!filesystem::create_directory(filesystem::path(out_dir))) {
        std::cout << "Error creating directory " << out_dir << std::endl;
        return 1;
    }
    if (povray_output) {
        if (!filesystem::create_directory(filesystem::path(pov_dir))) {
            std::cout << "Error creating directory " << pov_dir << std::endl;
            return 1;
        }
        terrain.ExportMeshPovray(out_dir);
    }

    // ------------------------
    // Create the driver system
    // ------------------------

    MyDriver driver(my_hmmwv.GetVehicle(), 0.5);
    driver.Initialize();

    // ---------------
    // Simulation loop
    // ---------------

    if (debug_output) {
        GetLog() << "\n\n============ System Configuration ============\n";
        my_hmmwv.LogHardpointLocations();
    }

    // Number of simulation steps between miscellaneous events
    int render_steps = (int)std::ceil(render_step_size / step_size);
    int debug_steps = (int)std::ceil(debug_step_size / step_size);

    // Initialize simulation frame counter and simulation time
    int step_number = 0;
    int render_frame = 0;
    double time = 0;

    while (vis->Run()) {
        time = my_hmmwv.GetSystem()->GetChTime();

        // End simulation
        if (time >= t_end)
            break;

        vis->BeginScene();
        vis->Render();
        vis->EndScene();

        // Render scene and output POV-Ray data
        if (povray_output && step_number % render_steps == 0) {
            char filename[100];
            sprintf(filename, "%s/data_%03d.dat", pov_dir.c_str(), render_frame + 1);
            utils::WriteVisualizationAssets(my_hmmwv.GetSystem(), filename);
            render_frame++;
        }

        // Debug logging
        if (debug_output && step_number % debug_steps == 0) {
            GetLog() << "\n\n============ System Information ============\n";
            GetLog() << "Time = " << time << "\n\n";
            my_hmmwv.DebugLog(OUT_SPRINGS | OUT_SHOCKS | OUT_CONSTRAINTS);
        }

        // Driver inputs
        DriverInputs driver_inputs = driver.GetInputs();

        // Update modules (process inputs from other modules)
        driver.Synchronize(time);
        terrain.Synchronize(time);
        my_hmmwv.Synchronize(time, driver_inputs, terrain);
        vis->Synchronize(time, driver_inputs);

        // Advance simulation for one timestep for all modules
        driver.Advance(step_size);
        terrain.Advance(step_size);
        my_hmmwv.Advance(step_size);
        vis->Advance(step_size);

        // Increment frame number
        step_number++;
    }

    return 0;
}

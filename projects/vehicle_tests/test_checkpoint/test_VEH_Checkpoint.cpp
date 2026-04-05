// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2026 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Radu Serban
// =============================================================================
//
// Demonstration of using checkpoints with Chrono::Vehicle simulations.
//
// The vehicle reference frame has Z up, X towards the front of the vehicle, and
// Y pointing to the left.
//
// =============================================================================

#include <cmath>

#include "chrono/physics/ChSystemSMC.h"

#include "chrono/functions/ChFunction.h"
#include "chrono/utils/ChUtils.h"
#include "chrono/input_output/ChWriterCSV.h"
#include "chrono/solver/ChSolverBB.h"

#include "chrono_vehicle/ChWorldFrame.h"
#include "chrono_vehicle/terrain/RigidTerrain.h"
#include "chrono_vehicle/ChVehicleVisualSystem.h"

#ifdef CHRONO_VSG
    #include "chrono_vehicle/wheeled_vehicle/ChWheeledVehicleVisualSystemVSG.h"
#endif

#ifdef CHRONO_POSTPROCESS
    #include "chrono_postprocess/ChGnuPlot.h"
#endif

#include "chrono_thirdparty/filesystem/path.h"

#include "demos/vehicle/WheeledVehicleModels.h"

using std::cout;
using std::cerr;
using std::endl;

// =============================================================================

// Create a Chrono physical system with default settings for wheeled vehicle simulations.
std::unique_ptr<ChSystemSMC> CreateSystem() {
    auto sys = chrono_types::make_unique<ChSystemSMC>();

    sys->SetGravitationalAcceleration(-9.81 * ChWorldFrame::Vertical());
    sys->SetCollisionSystemType(ChCollisionSystem::Type::BULLET);

    auto solver = chrono_types::make_shared<ChSolverBB>();
    solver->SetMaxIterations(100);
    solver->SetOmega(0.8);
    solver->SetSharpnessLambda(1.0);
    sys->SetSolver(solver);

    sys->SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT_LINEARIZED);

    return std::move(sys);
}

// =============================================================================

// Create a rigid terrain patch at the specified height.
RigidTerrain CreateTerrain(ChSystem* sys, double height) {
    // Terrain patch dimensions
    double x_size = 200;
    double y_size = 200;

    // Create the terrain
    ChContactMaterialData minfo;
    minfo.mu = 0.9f;
    minfo.cr = 0.01f;
    minfo.Y = 2e7f;
    auto patch_mat = minfo.CreateMaterial(sys->GetContactMethod());

    RigidTerrain terrain(sys);
    auto patch = terrain.AddPatch(patch_mat, ChCoordsysd(ChVector3d(0, 0, height), QUNIT), x_size, y_size);
    patch->SetTexture(GetChronoDataFile("textures/checker2.png"), 20, 20);
    patch->SetColor(ChColor(0.8f, 0.8f, 0.5f));

    terrain.Initialize();

    return terrain;
}

// =============================================================================

// Create a vehicle in the given system from the specified model at the origin with no rotation.
ChWheeledVehicle& CreateVehicle(ChSystem* sys, std::shared_ptr<WheeledVehicleModel> vehicle_model) {
    vehicle_model->Create(sys, ChCoordsysd(ChVector3d(0, 0, 0.5), QUNIT), false);
    return vehicle_model->GetVehicle();
}

// Create a vehicle in the given system from the specified model, and initialize its state from checkpoint files.
// The position of all new bodies is offset by the specified amount.
ChWheeledVehicle& CreateVehicle(ChSystem* sys, std::shared_ptr<WheeledVehicleModel> vehicle_model, const std::string& out_dir, const ChVector2d& offset) {
    // Cache all bodies in the system
    std::set<ChBody*> pre_existing;
    for (auto& body : sys->GetBodies())
        pre_existing.insert(body.get());

    // Create the vehicle from the specified model
    auto& vehicle = CreateVehicle(sys, vehicle_model);

    // Initialize vehicle and tires from checkpoint files
    vehicle.ImportCheckpoint(ChCheckpoint::Format::ASCII, out_dir + "/vehicle_checkpoint.txt");
    int tire_id = 0;
    for (const auto& a : vehicle.GetAxles()) {
        for (const auto& w : a->GetWheels()) {
            if (w->GetTire()) {
                w->GetTire()->ImportCheckpoint(ChCheckpoint::Format::ASCII, out_dir + "/tire_" + std::to_string(tire_id++) + "_checkpoint.txt");
            }
        }
    }

    // Find new bodies in the system
    auto for_vehicle_bodies = [&](auto fn) {
        for (auto& body : sys->GetBodies())
            if (pre_existing.find(body.get()) == pre_existing.end())
                fn(body);
    };

    // Offset position of all new bodies
    ChVector3d old_pos = vehicle.GetPos();
    if (offset.Length() > 1e-6) {
        for_vehicle_bodies([&](auto& body) { body->SetPos(body->GetPos() + ChVector3d(offset.x(), offset.y(), 0)); });
    }

    return vehicle;
}

// =============================================================================

// Driver system that uses specified functions for throttle, braking, and steering inputs.
class FunctionDriver : public ChDriver {
  public:
    FunctionDriver(ChVehicle& vehicle) : ChDriver(vehicle), m_throttle_fun(nullptr), m_braking_fun(nullptr), m_steering_fun(nullptr) {}

    void SetFunctions(std::shared_ptr<ChFunction> throttle, std::shared_ptr<ChFunction> braking, std::shared_ptr<ChFunction> steering) {
        m_throttle_fun = throttle;
        m_braking_fun = braking;
        m_steering_fun = steering;
    }

    virtual void Synchronize(double time) override {
        DriverInputs driver_inputs = GetInputs();
        m_throttle = m_throttle_fun->GetVal(time);
        m_braking = m_braking_fun->GetVal(time);
        m_steering = m_steering_fun->GetVal(time);
    }

  private:
    std::shared_ptr<ChFunction> m_throttle_fun;
    std::shared_ptr<ChFunction> m_braking_fun;
    std::shared_ptr<ChFunction> m_steering_fun;
};

// Steering function that starts at zero and then follows a cosine wave after a specified delay.
class FunctionCosineSteering : public ChFunction {
  public:
    FunctionCosineSteering(double delay, double frequency) : delay(delay), frequency(frequency) {}
    FunctionCosineSteering(const FunctionCosineSteering& other) {
        delay = other.delay;
        frequency = other.frequency;
    }
    virtual FunctionCosineSteering* Clone() const override { return new FunctionCosineSteering(*this); }
    virtual double GetVal(double time) const override {
        if (time < delay)
            return 0;
        return 0.5 + 0.5 * std::cos(CH_2PI * frequency * (time - delay) + CH_PI);
    }
  private:
    double delay;
    double frequency;
};

// Create a driver that uses the specified functions for throttle, braking, and steering inputs.
FunctionDriver& CreateDriver(ChVehicle& vehicle, std::shared_ptr<ChFunction> throttle_fun, std::shared_ptr<ChFunction> braking_fun, std::shared_ptr<ChFunction> steering_fun) {
    auto* driver = new FunctionDriver(vehicle);
    driver->SetFunctions(throttle_fun, braking_fun, steering_fun);
    driver->Initialize();
    return *driver;
}

// Create a driver and initialize its state from a checkpoint file.
// The driver inputs are kept constant at the checkpoint values.
FunctionDriver& CreateDriver(ChVehicle& vehicle, const std::string& out_dir) {
    auto* driver = new FunctionDriver(vehicle);
    driver->ImportCheckpoint(ChCheckpoint::Format::ASCII, out_dir + "/driver_checkpoint.txt");

    // Keep driver inputs constant at the checkpoint values
    auto throttle_fun = chrono_types::make_shared<ChFunctionConst>(driver->GetThrottle());
    auto braking_fun = chrono_types::make_shared<ChFunctionConst>(driver->GetBraking());
    auto steering_fun = chrono_types::make_shared<ChFunctionConst>(driver->GetSteering());
    driver->SetFunctions(throttle_fun, braking_fun, steering_fun);

    driver->Initialize();
    return *driver;
}

// =============================================================================

#ifdef CHRONO_VSG

// Create a VSG visualization system for the specified vehicle and driver.
std::shared_ptr<ChWheeledVehicleVisualSystemVSG> CreateVisualization(ChVehicle& vehicle, ChDriver& driver, std::shared_ptr<WheeledVehicleModel> vehicle_model) {
    auto vis = chrono_types::make_shared<ChWheeledVehicleVisualSystemVSG>();
    vis->SetWindowTitle(vehicle_model->ModelName());
    vis->AttachVehicle(&vehicle);
    vis->AttachDriver(&driver);
    vis->SetChaseCamera(vehicle_model->TrackPoint(), vehicle_model->CameraDistance(), vehicle_model->CameraHeight());
    vis->SetWindowSize(1280, 800);
    vis->EnableSkyTexture(SkyMode::DOME);
    vis->SetCameraAngleDeg(40);
    vis->SetLightIntensity(1.0f);
    vis->SetLightDirection(1.8 * CH_PI_2, CH_PI_4);
    vis->EnableShadows();
    vis->Initialize();
    vis->ToggleAbsFrameVisibility();

    return vis;
}

// Render the visualization system at the specified time and render frame rate.
// Return false if the visualization system is closed.
bool RenderVisualization(std::shared_ptr<ChWheeledVehicleVisualSystemVSG> vis, double time, double render_fps) {
    static int render_frame = 0;

    if (vis) {
        if (!vis->Run())
            return false;
        if (time >= render_frame / render_fps) {
            vis->Render();
            render_frame++;
        }
    }

    return true;
}

#endif

// =============================================================================

// Simulate a single vehicle from the specified model for the specified duration, and save final checkpoint.
void SimulateSingle(std::shared_ptr<WheeledVehicleModel> vehicle_model, double time_end, const std::string& out_dir) {
    cout << "Simulate single vehicle: " << vehicle_model->ModelName() << endl;

    // Create the containing Chrono system
    auto sys = CreateSystem();
    double step_size = 2e-3;

    // Create the vehicle model
    auto& vehicle = CreateVehicle(sys.get(), vehicle_model);
    vehicle.EnableRealtime(true);

    // Create the terrain
    RigidTerrain terrain = CreateTerrain(sys.get(), 0.0);

    // Create the driver system
    auto throttle_fun = chrono_types::make_shared<ChFunctionConst>(0.5);
    auto braking_fun = chrono_types::make_shared<ChFunctionConst>(0.0);
    auto steering_fun = chrono_types::make_shared<FunctionCosineSteering>(1.0, 0.5);
    auto& driver = CreateDriver(vehicle, throttle_fun, braking_fun, steering_fun);

    // Create the vehicle run-time visualization
#ifdef CHRONO_VSG
    auto vis = CreateVisualization(vehicle, driver, vehicle_model);
#endif

    // Simulation loop
    double render_fps = 50;

    while (true) {
        double time = sys->GetChTime();

#ifdef CHRONO_VSG
        if (!RenderVisualization(vis, time, render_fps))
            break;
#endif
        if (time > time_end)
            break;

        // Synchronize subsystems
        DriverInputs driver_inputs = driver.GetInputs();
        driver.Synchronize(time);
        terrain.Synchronize(time);
        vehicle_model->Synchronize(time, driver_inputs, terrain);
#ifdef CHRONO_VSG
        vis->Synchronize(time, driver_inputs);
#endif

        // Advance simulation
        driver.Advance(step_size);
        terrain.Advance(step_size);
        vehicle_model->Advance(step_size);
#ifdef CHRONO_VSG
        vis->Advance(step_size);
#endif

        // Advance state of containing system
        sys->DoStepDynamics(step_size);
    }

    // Checkpoint final vehicle and driver state
    cout << "Output vehicle checkpoint file: " << out_dir + "/vehicle_checkpoint.txt" << endl;
    cout << "Output tire checkpoint file:  " << out_dir + "/tire_X_checkpoint.txt" << endl;
    cout << "Output driver checkpoint file:  " << out_dir + "/driver_checkpoint.txt" << endl;
    cout << endl;

    vehicle.ExportCheckpoint(ChCheckpoint::Format::ASCII, out_dir + "/vehicle_checkpoint.txt");
    driver.ExportCheckpoint(ChCheckpoint::Format::ASCII, out_dir + "/driver_checkpoint.txt");
    int tire_id = 0;
    for (const auto& a : vehicle.GetAxles()) {
        for (const auto& w : a->GetWheels()) {
            if (w->GetTire()) {
                w->GetTire()->ExportCheckpoint(ChCheckpoint::Format::ASCII, out_dir + "/tire_" + std::to_string(tire_id++) + "_checkpoint.txt");
            }
        }
    }
}

// =============================================================================

// Simulate multiple vehicles from the specified model, initialized from the same checkpoint, and offset in space by the specified amounts.
void SimulateMultiple(std::shared_ptr<WheeledVehicleModel> vehicle_model, const std::vector<ChVector2d>& offsets, const std::string& out_dir) {
    cout << "Simulate multiple vehicles: " << vehicle_model->ModelName() << endl;

    // Create the containing Chrono system
    double step_size = 2e-3;
    auto sys = CreateSystem();

    // Create driver input functions
    auto throttle_fun = chrono_types::make_shared<ChFunctionConst>(0.5);
    auto braking_fun = chrono_types::make_shared<ChFunctionConst>(0.0);
    auto steering_fun = chrono_types::make_shared<ChFunctionConst>(0.0);

    // Create the vehicle models and associated driver systems
    std::vector<ChWheeledVehicle*> vehicles;
    std::vector<ChDriver*> drivers;
    std::vector<DriverInputs> driver_inputs;
    std::vector<ChWriterCSV> csv_writers;
    for (auto offset : offsets) {
        // Create the vehicle, initialize from checkpoint, and apply offset
        auto& vehicle = CreateVehicle(sys.get(), vehicle_model, out_dir, offset);
        vehicle.EnableRealtime(true);

        // Create the driver system and initialize from checkpoint
        auto& driver = CreateDriver(vehicle, out_dir);

        vehicles.push_back(&vehicle);
        drivers.push_back(&driver);
        driver_inputs.push_back(DriverInputs());
        csv_writers.push_back(ChWriterCSV(" "));
    }

    // Create the terrain
    RigidTerrain terrain = CreateTerrain(sys.get(), 0.0);

    // Create the vehicle run-time visualization
#ifdef CHRONO_VSG
    auto vis = CreateVisualization(*vehicles[0], *drivers[0], vehicle_model);
#else
    double time_end = 50;
#endif

    // Simulation loop
    double render_fps = 50;

    while (true) {
        double time = sys->GetChTime();

        for (int i = 0; i < vehicles.size(); i++) {
            csv_writers[i] << time << vehicles[i]->GetPos();
            csv_writers[i] << vehicles[i]->GetSpindlePos(0, VehicleSide::LEFT) << vehicles[i]->GetSpindlePos(0, VehicleSide::RIGHT);
            csv_writers[i] << vehicles[i]->GetSpindlePos(1, VehicleSide::LEFT) << vehicles[i]->GetSpindlePos(1, VehicleSide::RIGHT);
            csv_writers[i] << std::endl; 
        }

#ifdef CHRONO_VSG
        if (!RenderVisualization(vis, time, render_fps))
            break;
#else
        if (time > time_end)
            break;
#endif

        // Synchronize subsystems
        for (int i = 0; i < vehicles.size(); i++) {
            driver_inputs[i] = drivers[i]->GetInputs();
            drivers[i]->Synchronize(time);
            vehicles[i]->Synchronize(time, driver_inputs[i], terrain);
        }
        terrain.Synchronize(time);
#ifdef CHRONO_VSG
        vis->Synchronize(time, driver_inputs[0]);
#endif

        // Advance simulation
        for (int i = 0; i < vehicles.size(); i++) {
            drivers[i]->Advance(step_size);
            vehicles[i]->Advance(step_size);
        }
        terrain.Advance(step_size);
#ifdef CHRONO_VSG
        vis->Advance(step_size);
#endif

        // Advance state of entire system (containing all vehicles)
        sys->DoStepDynamics(step_size);
    }

#ifdef CHRONO_POSTPROCESS
    for (int i = 0; i < vehicles.size(); i++) {
        csv_writers[i].WriteToFile(out_dir + "/vehicle_" + std::to_string(i) + ".csv");
    }

    // Plot trajectories of all vehicles
    {
        postprocess::ChGnuPlot gplot_traj(out_dir + "/vehicle_traj.gpl");
        for (int i = 0; i < vehicles.size(); i++) {
            gplot_traj.Plot(out_dir + "/vehicle_" + std::to_string(i) + ".csv", 3, 2, "vehicle " + std::to_string(i), "with lines lw 2");
        }
        gplot_traj.SetTitle("Vehicle Trajectories");
        gplot_traj.SetLabelX("Y [m]");
        gplot_traj.SetLabelY("X [m]");
        gplot_traj.SetAxesEqual(true);
    }

    // Plot front-left spindle vertical positions for all vehicles
    {
        postprocess::ChGnuPlot gplot_fl(out_dir + "/vehicle_fl_spindle.gpl");
        for (int i = 0; i < vehicles.size(); i++) {
            gplot_fl.Plot(out_dir + "/vehicle_" + std::to_string(i) + ".csv", 1, 7, "vehicle " + std::to_string(i), "with lines lw 2");
        }
        gplot_fl.SetTitle("Front-left spindle vertical position");
        gplot_fl.SetLabelX("t [s]");
        gplot_fl.SetLabelY("Z [m]");
        gplot_fl.SetAxesEqual(true);
    }

#endif
}

// =============================================================================

int main(int argc, char* argv[]) {
    // Set path to Chrono and Chrono::Vehicle data directories
    SetChronoDataPath(CHRONO_DATA_DIR);
    SetVehicleDataPath(CHRONO_VEHICLE_DATA_DIR);

    // Vehicle model
    auto v = chrono_types::make_shared<HMMWV_Model>();

    // Create output directories
    std::string out_dir = GetChronoOutputPath() + "VEHICLE_CHECKPOINT_1";
    if (!filesystem::create_directory(filesystem::path(out_dir))) {
        cout << "Error creating directory " << out_dir << endl;
        return 1;
    }

    // Simulate vehicle for specified duration
    SimulateSingle(v, 3.5, out_dir);

    // Create multiple vehicles initialized from the same checkpoint and offset in space
    SimulateMultiple(v, {ChVector2d(0, -2), ChVector2d(0, +2)}, out_dir);

    return 0;
}

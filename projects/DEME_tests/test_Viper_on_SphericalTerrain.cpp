// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2021 projectchrono.org
// All right reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Ruochun Zhang
// =============================================================================
//
// Demo to show Viper Rover operated on DEM terrain represented by spherical
// particles, which serves as a comparison against GRC terrain
//
// =============================================================================

#include "chrono_models/robot/viper/Viper.h"

#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#include "chrono_vehicle/ChVehicleModelData.h"

#include <DEM/API.h>
#include <DEM/HostSideHelpers.hpp>
#include <DEM/utils/Samplers.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <map>
#include <random>

using namespace deme;
using namespace std::filesystem;

using namespace chrono;
using namespace chrono::geometry;
using namespace chrono::viper;

// Define Viper rover wheel type
ViperWheelType wheel_type = ViperWheelType::RealWheel;

// ChVector to/from float3
inline float3 ChVec2Float(const ChVector<>& vec) {
    return make_float3(vec.x(), vec.y(), vec.z());
}
inline ChVector<> Float2ChVec(float3 f3) {
    return ChVector<>(f3.x, f3.y, f3.z);
}

inline float4 ChQ2Float(const ChQuaternion<>& Q) {
    float4 f4;
    f4.w = Q.e0();
    f4.x = Q.e1();
    f4.y = Q.e2();
    f4.z = Q.e3();
    return f4;
}

void SaveParaViewFiles(Viper& rover, path& rover_dir, unsigned int frame_number);

int main(int argc, char* argv[]) {
    GetLog() << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";

    // Global parameter for moving patch size:
    double wheel_range = 0.5;
    ////double body_range = 1.2;

    // Create a Chrono::Engine physical system
    ChSystemSMC sys;
    ChVector<double> G = ChVector<double>(0, 0, -9.81);
    sys.Set_G_acc(G);

    const int nW = 4;  // 4 wheels

    // Create the rover
    // auto driver = chrono_types::make_shared<ViperDCMotorControl>();
    auto driver = chrono_types::make_shared<ViperSpeedDriver>(0, 3.14159 / 2);
    Viper viper(&sys, wheel_type);
    viper.SetDriver(driver);

    // driver->SetMotorNoLoadSpeed(0.8, ViperWheelID::V_LF);
    // driver->SetMotorNoLoadSpeed(0.8, ViperWheelID::V_RF);
    // driver->SetMotorNoLoadSpeed(0.8, ViperWheelID::V_LB);
    // driver->SetMotorNoLoadSpeed(0.8, ViperWheelID::V_RB);
    // driver->SetMotorStallTorque(50.0, ViperWheelID::V_LF);
    // driver->SetMotorStallTorque(50.0, ViperWheelID::V_RF);
    // driver->SetMotorStallTorque(50.0, ViperWheelID::V_LB);
    // driver->SetMotorStallTorque(50.0, ViperWheelID::V_RB);

    viper.Initialize(ChFrame<>(ChVector<>(-0.5, -0.0, -0.12), QUNIT));

    // Get wheels and bodies to set up SCM patches
    std::vector<std::shared_ptr<ChBodyAuxRef>> Wheels;
    std::vector<ChVector<>> wheel_pos;
    Wheels.push_back(viper.GetWheel(ViperWheelID::V_LF)->GetBody());
    Wheels.push_back(viper.GetWheel(ViperWheelID::V_RF)->GetBody());
    Wheels.push_back(viper.GetWheel(ViperWheelID::V_LB)->GetBody());
    Wheels.push_back(viper.GetWheel(ViperWheelID::V_RB)->GetBody());

    auto Body_1 = viper.GetChassis()->GetBody();
    double total_mass = viper.GetRoverMass();
    double wheel_mass = viper.GetWheelMass();

    for (int i = 0; i < nW; i++) {
        wheel_pos.push_back(Wheels[i]->GetFrame_REF_to_abs().GetPos());
    }

    //////////////////////////////////////////////
    // Now step up DEME
    //////////////////////////////////////////////

    DEMSolver DEMSim;
    DEMSim.SetVerbosity(INFO);
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    // DEMSim.SetOutputContent(OUTPUT_CONTENT::FAMILY);
    DEMSim.SetOutputContent(OUTPUT_CONTENT::XYZ);

    srand(759);

    float kg_g_conv = 1;
    // Define materials
    auto mat_type_terrain =
        DEMSim.LoadMaterial({{"E", 1e9 * kg_g_conv}, {"nu", 0.3}, {"CoR", 0.3}, {"mu", 0.5}, {"Crr", 0.0}});
    auto mat_type_wheel =
        DEMSim.LoadMaterial({{"E", 1e9 * kg_g_conv}, {"nu", 0.3}, {"CoR", 0.3}, {"mu", 0.5}, {"Crr", 0.0}});

    // Define the simulation world
    double world_y_size = 2.0;
    DEMSim.InstructBoxDomainDimension(2. * world_y_size, world_y_size, world_y_size);
    float bottom = -0.5;
    DEMSim.AddBCPlane(make_float3(0, 0, bottom), make_float3(0, 0, 1), mat_type_terrain);
    DEMSim.AddBCPlane(make_float3(0, world_y_size / 2, 0), make_float3(0, -1, 0), mat_type_terrain);
    DEMSim.AddBCPlane(make_float3(0, -world_y_size / 2, 0), make_float3(0, 1, 0), mat_type_terrain);
    // X-dir bounding planes
    DEMSim.AddBCPlane(make_float3(-world_y_size * 2 / 2, 0, 0), make_float3(1, 0, 0), mat_type_terrain);
    DEMSim.AddBCPlane(make_float3(world_y_size * 2 / 2, 0, 0), make_float3(-1, 0, 0), mat_type_terrain);

    // Define the wheel geometry
    float wheel_rad = 0.25;
    float wheel_width = 0.25;
    wheel_mass *= kg_g_conv;  // in kg or g
    // Our shelf wheel geometry is lying flat on ground with z being the axial direction
    float wheel_IYY = wheel_mass * wheel_rad * wheel_rad / 2;
    float wheel_IXX = (wheel_mass / 12) * (3 * wheel_rad * wheel_rad + wheel_width * wheel_width);
    float3 wheel_MOI = make_float3(wheel_IXX, wheel_IYY, wheel_IXX);
    auto wheel_template =
        DEMSim.LoadClumpType(wheel_mass, wheel_MOI, "../data/clumps/ViperWheelSimple.csv", mat_type_wheel);
    // The file contains no wheel particles size info, so let's manually set them
    wheel_template->radii = std::vector<float>(wheel_template->nComp, 0.01);
    // This wheel template is `lying down', but our reported MOI info is assuming it's in a position to roll
    // along X direction. Let's make it clear its principal axes is not what we used to report its component
    // sphere relative positions.
    wheel_template->InformCentroidPrincipal(make_float3(0), make_float4(0.7071, 0, 0, 0.7071));

    // Then the ground particle template
    float sp_rad = 0.005;
    // float sp_rad = 0.0025; // This one is fine too, with 5e-6 step size
    auto ground_particle_template =
        DEMSim.LoadSphereType(4. / 3. * 3.14 * std::pow(sp_rad, 3) * 2.6e3, sp_rad, mat_type_terrain);
    ground_particle_template->SetVolume(4. / 3. * 3.14 * std::pow(sp_rad, 3));

    // Now we load part1 clump locations from a output file
    std::cout << "Making terrain..." << std::endl;
    auto in_xyz = DEMBoxHCPSampler(make_float3(0, 0, -0.45),
                                   make_float3(world_y_size - 2 * sp_rad, world_y_size / 2 - sp_rad, 0.05 - sp_rad),
                                   2.01 * sp_rad);

    DEMSim.AddClumps(std::vector<std::shared_ptr<DEMClumpTemplate>>(in_xyz.size(), ground_particle_template), in_xyz);

    ////////////////////
    // Add wheel in DEM
    ////////////////////

    // Instantiate this wheel
    std::cout << "Making wheels..." << std::endl;
    DEMSim.SetFamilyFixed(100);
    std::vector<std::shared_ptr<DEMTracker>> trackers;
    std::vector<std::shared_ptr<DEMClumpBatch>> DEM_Wheels;
    for (int i = 0; i < nW; i++) {
        DEM_Wheels.push_back(
            DEMSim.AddClumps(wheel_template, make_float3(wheel_pos[i].x(), wheel_pos[i].y(), wheel_pos[i].z())));
        // DEM_Wheel->SetOriQ(make_float4(0.7071, 0.7071, 0, 0));
        DEM_Wheels[i]->SetFamily(100);
        trackers.push_back(DEMSim.Track(DEM_Wheels[i]));
    }
    DEMSim.DisableFamilyOutput(100);  // no need outputting wheels

    //////
    // Make ready for DEM simulation
    ///////
    auto max_v_finder = DEMSim.CreateInspector("clump_max_absv");
    auto max_z_finder = DEMSim.CreateInspector("clump_max_z");
    auto void_ratio_finder =
        DEMSim.CreateInspector("clump_volume", "return (abs(X) <= 0.48) && (abs(Y) <= 0.48) && (Z <= -0.45);");
    float total_volume = 0.96 * 0.96 * 0.05;

    // Now add a plane to compress the `road'
    auto compressor = DEMSim.AddExternalObject();
    compressor->AddPlane(make_float3(0, 0, 0), make_float3(0, 0, -1), mat_type_terrain);
    compressor->SetFamily(90);
    DEMSim.SetFamilyFixed(90);
    DEMSim.DisableContactBetweenFamilies(90, 100);  // no contact between compressor and wheels
    auto compressor_tracker = DEMSim.Track(compressor);

    float base_step_size = 2e-5;
    float step_size = base_step_size;
    float base_vel = 0.4;
    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetGravitationalAcceleration(ChVec2Float(G));
    // If you want to use a large UpdateFreq then you have to expand spheres to ensure safety
    DEMSim.SetCDUpdateFreq(20);
    DEMSim.SetMaxVelocity(40.);
    DEMSim.SetInitBinSize(sp_rad * 6);
    DEMSim.SetIntegrator(TIME_INTEGRATOR::EXTENDED_TAYLOR);

    DEMSim.Initialize();
    for (const auto& tracker : trackers) {
        std::cout << "A tracker is tracking owner " << tracker->obj->ownerID << std::endl;
    }

    ///////////////////////////////////////////
    // Compress the road first
    ///////////////////////////////////////////

    float time_end = 8.0;
    unsigned int fps = 20;
    unsigned int report_freq = 5000;
    unsigned int param_update_freq = 5000;
    unsigned int out_steps = (unsigned int)(1.0 / (fps * step_size));
    float frame_accu_thres = 1.0 / fps;
    unsigned int report_steps = (unsigned int)(1.0 / (report_freq * step_size));
    unsigned int param_update_steps = (unsigned int)(1.0 / (param_update_freq * step_size));

    path out_dir = current_path();
    out_dir += "/Viper_on_GRC_flat";
    path rover_dir = out_dir / "./rover";
    create_directory(out_dir);
    create_directory(rover_dir);
    unsigned int currframe = 0;
    unsigned int curr_step = 0;

    {
        // Figure out the exact size of wheels
        float max_z = max_z_finder->GetValue();
        float wheel_center_z = ChVec2Float(Wheels[0]->GetFrame_REF_to_abs().GetPos()).z;
        std::cout << "Exact wheel radius is " << max_z - wheel_center_z << std::endl;
    }

    step_size = 5e-6;
    DEMSim.SetInitTimeStep(step_size);
    DEMSim.UpdateStepSize();

    double now_z = -0.37;
    compressor_tracker->SetPos(make_float3(0, 0, now_z));
    float compress_time = 0.4;
    double compressor_final_dist = 0.4 - 0.37;
    double compressor_v = compressor_final_dist / compress_time;
    for (float t = 0; t < compress_time; t += step_size, curr_step++) {
        if (curr_step % out_steps == 0) {
            std::cout << "Frame: " << currframe << std::endl;
            DEMSim.ShowThreadCollaborationStats();
            char filename[200];
            sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe);
            DEMSim.WriteSphereFile(std::string(filename));
            SaveParaViewFiles(viper, rover_dir, currframe);
            currframe++;
        }
        now_z -= compressor_v * step_size;
        compressor_tracker->SetPos(make_float3(0, 0, now_z));
        DEMSim.DoDynamics(step_size);
    }
    for (float t = 0; t < compress_time; t += step_size, curr_step++) {
        if (curr_step % out_steps == 0) {
            std::cout << "Frame: " << currframe << std::endl;
            DEMSim.ShowThreadCollaborationStats();
            char filename[200];
            sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe);
            DEMSim.WriteSphereFile(std::string(filename));
            SaveParaViewFiles(viper, rover_dir, currframe);
            currframe++;
        }
        now_z += compressor_v * step_size;
        compressor_tracker->SetPos(make_float3(0, 0, now_z));
        DEMSim.DoDynamics(step_size);
    }

    DEMSim.DisableContactBetweenFamilies(90, 0);
    DEMSim.DoDynamicsThenSync(0);
    step_size = 2e-5;
    DEMSim.SetInitTimeStep(step_size);
    DEMSim.UpdateStepSize();

    ///////////////////////////////////////////
    // Real simulation
    ///////////////////////////////////////////

    // Timers
    std::chrono::high_resolution_clock::time_point h_start, d_start;
    std::chrono::high_resolution_clock::time_point h_end, d_end;
    std::chrono::duration<double> h_total, d_total;
    h_total = std::chrono::duration<double>(0);
    d_total = std::chrono::duration<double>(0);

    std::vector<ChQuaternion<>> wheel_rot(4);
    float max_v;
    int change_step = 0;
    float frame_accu = frame_accu_thres;
    for (float t = 0; t < time_end; t += step_size, curr_step++) {
        for (int i = 0; i < nW; i++) {
            wheel_pos[i] = Wheels[i]->GetFrame_REF_to_abs().GetPos();
            trackers[i]->SetPos(ChVec2Float(wheel_pos[i]));
            wheel_rot[i] = Wheels[i]->GetFrame_REF_to_abs().GetRot();
            trackers[i]->SetOriQ(ChQ2Float(wheel_rot[i]));

            if (curr_step % report_steps == 0) {
                std::cout << "Wheel " << i << " position: " << wheel_pos[i].x() << ", " << wheel_pos[i].y() << ", "
                          << wheel_pos[i].z() << std::endl;

                std::cout << "Wheel " << i << " rotation: " << wheel_rot[i].e0() << ", " << wheel_rot[i].e1() << ", "
                          << wheel_rot[i].e2() << ", " << wheel_rot[i].e3() << std::endl;
            }
        }

        // if (curr_step % out_steps == 0) {
        if (frame_accu >= frame_accu_thres) {
            frame_accu = 0.;
            std::cout << "Frame: " << currframe << std::endl;
            DEMSim.ShowThreadCollaborationStats();
            char filename[200];
            sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe);
            DEMSim.WriteSphereFile(std::string(filename));
            SaveParaViewFiles(viper, rover_dir, currframe);
            currframe++;
        }
        // Run DEM first
        d_start = std::chrono::high_resolution_clock::now();
        DEMSim.DoDynamics(step_size);
        d_end = std::chrono::high_resolution_clock::now();
        d_total += d_end - d_start;

        // Then feed force
        for (int i = 0; i < nW; i++) {
            float3 F = trackers[i]->ContactAcc();
            F *= wheel_mass;
            float3 tor = trackers[i]->ContactAngAccLocal();
            tor = wheel_MOI * tor;
            Wheels[i]->Empty_forces_accumulators();
            Wheels[i]->Accumulate_force(Float2ChVec(F), wheel_pos[i], false);
            Wheels[i]->Accumulate_torque(Float2ChVec(tor), true);  // torque in SMUG is local
        }
        h_start = std::chrono::high_resolution_clock::now();
        sys.DoStepDynamics(step_size);
        viper.Update();
        h_end = std::chrono::high_resolution_clock::now();
        h_total += h_end - h_start;

        t += step_size;
        frame_accu += step_size;

        if (curr_step % report_steps == 0) {
            float3 body_pos = ChVec2Float(Body_1->GetFrame_REF_to_abs().GetPos());
            std::cout << "Rover body is at " << body_pos.x << ", " << body_pos.y << ", " << body_pos.z << std::endl;
            std::cout << "Time is " << t << std::endl;
            max_v = max_v_finder->GetValue();
            std::cout << "Max vel in simulation is " << max_v << std::endl;
            float matter_volume = void_ratio_finder->GetValue();
            std::cout << "Void ratio now: " << (total_volume - matter_volume) / matter_volume << std::endl;
            std::cout << "========================" << std::endl;
        }
    }
    std::cout << h_total.count() << " seconds spent on host" << std::endl;
    std::cout << d_total.count() << " seconds spent on device" << std::endl;

    DEMSim.ShowThreadCollaborationStats();
    DEMSim.ShowTimingStats();

    return 0;
}

//------------------------------------------------------------------
// Function to save the povray files of the MBD
//------------------------------------------------------------------
void SaveParaViewFiles(Viper& rover, path& rover_dir, unsigned int frame_number) {
    path filename;

    char f_name[20];
    sprintf(f_name, "%04d", frame_number);
    filename = rover_dir / ("./viper_" + std::string(f_name) + ".obj");

    std::vector<geometry::ChTriangleMeshConnected> meshes;

    // save the VIPER body to obj/vtk files
    for (int i = 0; i < 1; i++) {
        auto body = rover.GetChassis()->GetBody();
        ChFrame<> body_ref_frame = body->GetFrame_REF_to_abs();
        ChVector<> body_pos = body_ref_frame.GetPos();      // body->GetPos();
        ChQuaternion<> body_rot = body_ref_frame.GetRot();  // body->GetRot();

        auto mmesh = chrono_types::make_shared<ChTriangleMeshConnected>();
        std::string obj_path = (GetChronoDataFile("robot/viper/obj/viper_chassis.obj"));
        double scale_ratio = 1.0;
        mmesh->LoadWavefrontMesh(obj_path, false, true);
        mmesh->Transform(ChVector<>(0, 0, 0), ChMatrix33<>(scale_ratio));  // scale to a different size
        mmesh->RepairDuplicateVertexes(1e-9);                              // if meshes are not watertight

        double mmass;
        ChVector<> mcog;
        ChMatrix33<> minertia;
        mmesh->ComputeMassProperties(true, mmass, mcog, minertia);
        mmesh->Transform(body_pos, ChMatrix33<>(body_rot));  // rotate the mesh based on the orientation of body

        // filename = rover_dir / ("./body_" + std::string(f_name) + ".obj");
        // std::vector<geometry::ChTriangleMeshConnected> meshes = {*mmesh};
        // geometry::ChTriangleMeshConnected::WriteWavefront(filename.string(), meshes);
        meshes.push_back(*mmesh);
    }

    // save the wheels to obj/vtk files
    for (int i = 0; i < 4; i++) {
        std::shared_ptr<ChBodyAuxRef> body;
        if (i == 0) {
            body = rover.GetWheel(ViperWheelID::V_LF)->GetBody();
        }
        if (i == 1) {
            body = rover.GetWheel(ViperWheelID::V_RF)->GetBody();
        }
        if (i == 2) {
            body = rover.GetWheel(ViperWheelID::V_LB)->GetBody();
        }
        if (i == 3) {
            body = rover.GetWheel(ViperWheelID::V_RB)->GetBody();
        }

        ChFrame<> body_ref_frame = body->GetFrame_REF_to_abs();
        ChVector<> body_pos = body_ref_frame.GetPos();      // body->GetPos();
        ChQuaternion<> body_rot = body_ref_frame.GetRot();  // body->GetRot();
        if (i == 0 || i == 2) {
            body_rot.Cross(body_rot, Q_from_AngZ(CH_C_PI));
        }

        auto mmesh = chrono_types::make_shared<ChTriangleMeshConnected>();
        std::string obj_path = GetChronoDataFile("robot/viper/obj/viper_wheel.obj");
        double scale_ratio = 1.0;
        mmesh->LoadWavefrontMesh(obj_path, false, true);
        mmesh->Transform(ChVector<>(0, 0, 0), ChMatrix33<>(scale_ratio));  // scale to a different size
        mmesh->RepairDuplicateVertexes(1e-9);                              // if meshes are not watertight

        double mmass;
        ChVector<> mcog;
        ChMatrix33<> minertia;
        mmesh->ComputeMassProperties(true, mmass, mcog, minertia);
        mmesh->Transform(body_pos, ChMatrix33<>(body_rot));  // rotate the mesh based on the orientation of body

        // filename = rover_dir / ("./wheel_" + std::to_string(i + 1) + "_" + std::string(f_name) + ".obj");
        // std::vector<geometry::ChTriangleMeshConnected> meshes = {*mmesh};
        // geometry::ChTriangleMeshConnected::WriteWavefront(filename.string(), meshes);
        meshes.push_back(*mmesh);
    }

    // save the steering rod to obj/vtk files
    for (int i = 0; i < 4; i++) {
        std::shared_ptr<ChBodyAuxRef> body;
        if (i == 0) {
            body = rover.GetUpright(ViperWheelID::V_LF)->GetBody();
        }
        if (i == 1) {
            body = rover.GetUpright(ViperWheelID::V_RF)->GetBody();
        }
        if (i == 2) {
            body = rover.GetUpright(ViperWheelID::V_LB)->GetBody();
        }
        if (i == 3) {
            body = rover.GetUpright(ViperWheelID::V_RB)->GetBody();
        }
        ChFrame<> body_ref_frame = body->GetFrame_REF_to_abs();
        ChVector<> body_pos = body_ref_frame.GetPos();      // body->GetPos();
        ChQuaternion<> body_rot = body_ref_frame.GetRot();  // body->GetRot();

        auto mmesh = chrono_types::make_shared<ChTriangleMeshConnected>();
        std::string obj_path = "";
        if (i == 0 || i == 2) {
            obj_path = GetChronoDataFile("robot/viper/obj/viper_L_steer.obj");
        }
        if (i == 1 || i == 3) {
            obj_path = GetChronoDataFile("robot/viper/obj/viper_R_steer.obj");
        }
        double scale_ratio = 1.0;
        mmesh->LoadWavefrontMesh(obj_path, false, true);
        mmesh->Transform(ChVector<>(0, 0, 0), ChMatrix33<>(scale_ratio));  // scale to a different size
        mmesh->RepairDuplicateVertexes(1e-9);                              // if meshes are not watertight

        double mmass;
        ChVector<> mcog;
        ChMatrix33<> minertia;
        mmesh->ComputeMassProperties(true, mmass, mcog, minertia);
        mmesh->Transform(body_pos, ChMatrix33<>(body_rot));  // rotate the mesh based on the orientation of body

        // filename = rover_dir / ("./steerRod_" + std::to_string(i + 1) + "_" + std::string(f_name) + ".obj");
        // std::vector<geometry::ChTriangleMeshConnected> meshes = {*mmesh};
        // geometry::ChTriangleMeshConnected::WriteWavefront(filename.string(), meshes);
        meshes.push_back(*mmesh);
    }

    // save the lower rod to obj/vtk files
    for (int i = 0; i < 4; i++) {
        std::shared_ptr<ChBodyAuxRef> body;
        if (i == 0) {
            body = rover.GetLowerArm(ViperWheelID::V_LF)->GetBody();
        }
        if (i == 1) {
            body = rover.GetLowerArm(ViperWheelID::V_RF)->GetBody();
        }
        if (i == 2) {
            body = rover.GetLowerArm(ViperWheelID::V_LB)->GetBody();
        }
        if (i == 3) {
            body = rover.GetLowerArm(ViperWheelID::V_RB)->GetBody();
        }
        ChFrame<> body_ref_frame = body->GetFrame_REF_to_abs();
        ChVector<> body_pos = body_ref_frame.GetPos();      // body->GetPos();
        ChQuaternion<> body_rot = body_ref_frame.GetRot();  // body->GetRot();

        auto mmesh = chrono_types::make_shared<ChTriangleMeshConnected>();
        std::string obj_path = "";
        if (i == 0 || i == 2) {
            obj_path = GetChronoDataFile("robot/viper/obj/viper_L_bt_sus.obj");
        }
        if (i == 1 || i == 3) {
            obj_path = GetChronoDataFile("robot/viper/obj/viper_R_bt_sus.obj");
        }
        double scale_ratio = 1.0;
        mmesh->LoadWavefrontMesh(obj_path, false, true);
        mmesh->Transform(ChVector<>(0, 0, 0), ChMatrix33<>(scale_ratio));  // scale to a different size
        mmesh->RepairDuplicateVertexes(1e-9);                              // if meshes are not watertight

        double mmass;
        ChVector<> mcog;
        ChMatrix33<> minertia;
        mmesh->ComputeMassProperties(true, mmass, mcog, minertia);
        mmesh->Transform(body_pos, ChMatrix33<>(body_rot));  // rotate the mesh based on the orientation of body

        // filename = rover_dir / ("./lowerRod_" + std::to_string(i + 1) + "_" + std::string(f_name) + ".obj");
        // std::vector<geometry::ChTriangleMeshConnected> meshes = {*mmesh};
        // geometry::ChTriangleMeshConnected::WriteWavefront(filename.string(), meshes);
        meshes.push_back(*mmesh);
    }

    // save the upper rod to obj/vtk files
    for (int i = 0; i < 4; i++) {
        std::shared_ptr<ChBodyAuxRef> body;
        if (i == 0) {
            body = rover.GetUpperArm(ViperWheelID::V_LF)->GetBody();
        }
        if (i == 1) {
            body = rover.GetUpperArm(ViperWheelID::V_RF)->GetBody();
        }
        if (i == 2) {
            body = rover.GetUpperArm(ViperWheelID::V_LB)->GetBody();
        }
        if (i == 3) {
            body = rover.GetUpperArm(ViperWheelID::V_RB)->GetBody();
        }
        ChFrame<> body_ref_frame = body->GetFrame_REF_to_abs();
        ChVector<> body_pos = body_ref_frame.GetPos();      // body->GetPos();
        ChQuaternion<> body_rot = body_ref_frame.GetRot();  // body->GetRot();

        auto mmesh = chrono_types::make_shared<ChTriangleMeshConnected>();
        std::string obj_path = "";
        if (i == 0 || i == 2) {
            obj_path = GetChronoDataFile("robot/viper/obj/viper_L_up_sus.obj");
        }
        if (i == 1 || i == 3) {
            obj_path = GetChronoDataFile("robot/viper/obj/viper_R_up_sus.obj");
        }

        double scale_ratio = 1.0;
        mmesh->LoadWavefrontMesh(obj_path, false, true);
        mmesh->Transform(ChVector<>(0, 0, 0), ChMatrix33<>(scale_ratio));  // scale to a different size
        mmesh->RepairDuplicateVertexes(1e-9);                              // if meshes are not watertight

        double mmass;
        ChVector<> mcog;
        ChMatrix33<> minertia;
        mmesh->ComputeMassProperties(true, mmass, mcog, minertia);
        mmesh->Transform(body_pos, ChMatrix33<>(body_rot));  // rotate the mesh based on the orientation of body

        // filename = rover_dir / ("./upperRod_" + std::to_string(i + 1) + "_" + std::string(f_name) + ".obj");
        // std::vector<geometry::ChTriangleMeshConnected> meshes = {*mmesh};
        // geometry::ChTriangleMeshConnected::WriteWavefront(filename.string(), meshes);
        meshes.push_back(*mmesh);
    }

    geometry::ChTriangleMeshConnected::WriteWavefront(filename.string(), meshes);
}

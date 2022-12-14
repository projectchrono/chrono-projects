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
// Demo to show Rover operating on incline of GRC-1 simulant, with DEM-Engine
// providing the DEM simulation support
//
// =============================================================================

#include "chrono_models/robot/viper/Viper.h"

#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#include "chrono_vehicle/ChVehicleModelData.h"

#include <DEM/API.h>
// #include <core/utils/DEMEPaths.hpp>
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

const double math_PI = 3.1415927;

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
    SetChronoDataPath(CHRONO_DATA_DIR);
    SetDEMEDataPath(DEME_DATA_DIR);
    std::cout << "DEME dir is " << GetDEMEDataPath() << std::endl;

    // `World'
    float G_mag = 9.81;
    float step_size = 1e-6;

    // Define the wheel geometry
    float wheel_rad = 0.25;
    float wheel_width = 0.2;
    float wheel_mass = 11.;
    float total_pressure = 22. * 9.81;
    float added_pressure = (total_pressure - wheel_mass * G_mag);
    float wheel_IYY = wheel_mass * wheel_rad * wheel_rad / 2;
    float wheel_IXX = (wheel_mass / 12) * (3 * wheel_rad * wheel_rad + wheel_width * wheel_width);
    float3 wheel_MOI = make_float3(wheel_IXX, wheel_IYY, wheel_IXX);

    // std::string wheel_obj_path = GetChronoDataFile("robot/viper/obj/viper_wheel.obj");
    std::string wheel_obj_path = "./Moon_rover_wheel.obj";

    // Global parameter for moving patch size:
    double wheel_range = 0.5;
    ////double body_range = 1.2;

    // Create a Chrono::Engine physical system
    float Slope_deg = 5;
    double G_ang = Slope_deg * math_PI / 180.;
    ChSystemSMC sys;
    ChVector<double> G = ChVector<double>(-G_mag * std::sin(G_ang), 0, -G_mag * std::cos(G_ang));
    sys.Set_G_acc(G);

    const int nW = 4;  // 4 wheels

    // Create the rover
    // auto driver = chrono_types::make_shared<ViperDCMotorControl>();
    float w_r = 0.8;
    auto driver = chrono_types::make_shared<ViperSpeedDriver>(0, w_r);
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

    viper.Initialize(ChFrame<>(ChVector<>(-0.9, -0.0, -0.13), QUNIT));

    // Get wheels and bodies to set up SCM patches
    std::vector<std::shared_ptr<ChBodyAuxRef>> Wheels;
    std::vector<ChVector<>> wheel_pos;
    Wheels.push_back(viper.GetWheel(ViperWheelID::V_LF)->GetBody());
    Wheels.push_back(viper.GetWheel(ViperWheelID::V_RF)->GetBody());
    Wheels.push_back(viper.GetWheel(ViperWheelID::V_LB)->GetBody());
    Wheels.push_back(viper.GetWheel(ViperWheelID::V_RB)->GetBody());

    auto Body_1 = viper.GetChassis()->GetBody();
    std::cout << "Rover mass: " << viper.GetRoverMass() << std::endl;
    std::cout << "Wheel mass: " << viper.GetWheelMass() << std::endl;

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
    float m_cm_cov = 1;
    // Define materials
    auto mat_type_terrain =
        DEMSim.LoadMaterial({{"E", 1e9 * kg_g_conv / m_cm_cov}, {"nu", 0.3}, {"CoR", 0.5}, {"mu", 0.9}, {"Crr", 0.0}});
    auto mat_type_wheel =
        DEMSim.LoadMaterial({{"E", 1e9 * kg_g_conv / m_cm_cov}, {"nu", 0.3}, {"CoR", 0.5}, {"mu", 0.9}, {"Crr", 0.0}});

    // Define the simulation world
    double world_x_size = 4.0 * m_cm_cov;
    double world_y_size = 2.0 * m_cm_cov;
    DEMSim.InstructBoxDomainDimension(world_x_size, world_y_size, world_y_size);
    // DEMSim.InstructBoxDomainNumVoxel(22, 21, 21, (world_y_size) / std::pow(2, 16) / std::pow(2, 21));
    float bottom = -0.5 * m_cm_cov;
    DEMSim.AddBCPlane(make_float3(0, 0, bottom), make_float3(0, 0, 1), mat_type_terrain);
    DEMSim.AddBCPlane(make_float3(0, world_y_size / 2, 0), make_float3(0, -1, 0), mat_type_terrain);
    DEMSim.AddBCPlane(make_float3(0, -world_y_size / 2, 0), make_float3(0, 1, 0), mat_type_terrain);
    // X-dir bounding planes
    DEMSim.AddBCPlane(make_float3(-world_x_size / 2., 0, 0), make_float3(1, 0, 0), mat_type_terrain);
    DEMSim.AddBCPlane(make_float3(world_x_size / 2., 0, 0), make_float3(-1, 0, 0), mat_type_terrain);

    // Then the ground particle template
    DEMClumpTemplate shape_template;
    shape_template.ReadComponentFromFile(GetDEMEDataFile("clumps/triangular_flat.csv"));
    // Calculate its mass and MOI
    float mass = 2.6e3 * 5.5886717 * kg_g_conv;  // in kg or g
    float3 MOI = make_float3(1.8327927, 2.1580013, 0.77010059) * (double)2.6e3 * kg_g_conv / m_cm_cov / m_cm_cov;
    std::for_each(shape_template.radii.begin(), shape_template.radii.end(), [m_cm_cov](float& r) { r *= m_cm_cov; });
    std::for_each(shape_template.relPos.begin(), shape_template.relPos.end(), [m_cm_cov](float3& r) { r *= m_cm_cov; });
    // Scale the template we just created
    std::vector<std::shared_ptr<DEMClumpTemplate>> ground_particle_templates;
    std::vector<double> scales = {0.00063, 0.00033, 0.00022, 0.00015, 0.00009};
    std::for_each(scales.begin(), scales.end(), [](double& r) { r *= 20.; });
    for (double scaling : scales) {
        auto this_template = shape_template;
        this_template.mass = (double)mass * scaling * scaling * scaling;
        this_template.MOI.x = (double)MOI.x * (double)(scaling * scaling * scaling * scaling * scaling);
        this_template.MOI.y = (double)MOI.y * (double)(scaling * scaling * scaling * scaling * scaling);
        this_template.MOI.z = (double)MOI.z * (double)(scaling * scaling * scaling * scaling * scaling);
        std::cout << "Mass: " << this_template.mass << std::endl;
        std::cout << "MOIX: " << this_template.MOI.x << std::endl;
        std::cout << "MOIY: " << this_template.MOI.y << std::endl;
        std::cout << "MOIZ: " << this_template.MOI.z << std::endl;
        std::cout << "=====================" << std::endl;
        std::for_each(this_template.radii.begin(), this_template.radii.end(), [scaling](float& r) { r *= scaling; });
        std::for_each(this_template.relPos.begin(), this_template.relPos.end(), [scaling](float3& r) { r *= scaling; });
        this_template.materials = std::vector<std::shared_ptr<DEMMaterial>>(this_template.nComp, mat_type_terrain);
        ground_particle_templates.push_back(DEMSim.LoadClumpType(this_template));
    }

    // Now we load clump locations from a checkpointed file
    {
        std::cout << "Making terrain..." << std::endl;
        auto clump_xyz = DEMSim.ReadClumpXyzFromCsv("./GRC_20e6.csv");
        auto clump_quaternion = DEMSim.ReadClumpQuatFromCsv("./GRC_20e6.csv");
        std::vector<float3> in_xyz;
        std::vector<float4> in_quat;
        std::vector<std::shared_ptr<DEMClumpTemplate>> in_types;
        unsigned int t_num = 0;
        for (int i = 0; i < scales.size(); i++) {
            char t_name[20];
            sprintf(t_name, "%04d", t_num);

            auto this_type_xyz = clump_xyz[std::string(t_name)];
            auto this_type_quat = clump_quaternion[std::string(t_name)];

            size_t n_clump_this_type = this_type_xyz.size();
            std::cout << "Loading clump " << std::string(t_name) << " which has particle num: " << n_clump_this_type
                      << std::endl;
            // Prepare clump type identification vector for loading into the system (don't forget type 0 in
            // ground_particle_templates is the template for rover wheel)
            std::vector<std::shared_ptr<DEMClumpTemplate>> this_type(n_clump_this_type,
                                                                     ground_particle_templates.at(t_num));

            // Add them to the big long vector
            std::for_each(this_type_xyz.begin(), this_type_xyz.end(), [m_cm_cov](float3& r) { r *= m_cm_cov; });
            in_xyz.insert(in_xyz.end(), this_type_xyz.begin(), this_type_xyz.end());
            in_quat.insert(in_quat.end(), this_type_quat.begin(), this_type_quat.end());
            in_types.insert(in_types.end(), this_type.begin(), this_type.end());
            std::cout << "Added clump type " << t_num << std::endl;
            // Our template names are 0000, 0001 etc.
            t_num++;
        }
        // Finally, load the info into this batch
        DEMClumpBatch base_batch(in_xyz.size());
        base_batch.SetTypes(in_types);
        base_batch.SetPos(in_xyz);
        base_batch.SetOriQ(in_quat);

        DEMSim.AddClumps(base_batch);
    }

    // Right now, this demo loads particles in one batch, so no need to have the following batch-by-batch addition
    /*
    std::vector<float> x_shift_dist = {-1.5, -0.5, 0.5, 1.5};
    std::vector<float> y_shift_dist = {-0.5, 0.5};
    for (float x_shift : x_shift_dist) {
        for (float y_shift : y_shift_dist) {
            auto part1_clump_xyz = DEMSim.ReadClumpXyzFromCsv("./GRC.csv");
            auto part1_clump_quaternion = DEMSim.ReadClumpQuatFromCsv("./GRC.csv");
            std::vector<float3> in_xyz;
            std::vector<float4> in_quat;
            std::vector<std::shared_ptr<DEMClumpTemplate>> in_types;
            unsigned int t_num = 0;
            for (int i = 0; i < scales.size(); i++) {
                // Our template names are 0001, 0002 etc.
                t_num++;
                char t_name[20];
                sprintf(t_name, "%04d", t_num);

                auto this_type_xyz = part1_clump_xyz[std::string(t_name)];
                auto this_type_quat = part1_clump_quaternion[std::string(t_name)];

                size_t n_clump_this_type = this_type_xyz.size();
                std::cout << "Loading clump " << std::string(t_name) << " which has particle num: " << n_clump_this_type
                          << std::endl;
                // Prepare clump type identification vector for loading into the system (don't forget type 0 in
                // ground_particle_templates is the template for rover wheel)
                std::vector<std::shared_ptr<DEMClumpTemplate>> this_type(n_clump_this_type,
                                                                         ground_particle_templates.at(t_num - 1));

                // Add them to the big long vector
                std::for_each(this_type_xyz.begin(), this_type_xyz.end(), [m_cm_cov](float3& r) { r *= m_cm_cov; });
                in_xyz.insert(in_xyz.end(), this_type_xyz.begin(), this_type_xyz.end());
                in_quat.insert(in_quat.end(), this_type_quat.begin(), this_type_quat.end());
                in_types.insert(in_types.end(), this_type.begin(), this_type.end());
                std::cout << "Added clump type " << t_num << std::endl;
            }

            // Now, apply x_ and y-shifts
            std::for_each(in_xyz.begin(), in_xyz.end(), [x_shift, y_shift](float3& xyz) {
                xyz.x += x_shift;
                xyz.y += y_shift;
            });

            // Finally, load the info into this batch
            DEMClumpBatch base_batch(in_xyz.size());
            base_batch.SetTypes(in_types);
            base_batch.SetPos(in_xyz);
            base_batch.SetOriQ(in_quat);

            DEMSim.AddClumps(base_batch);
        }
    }
    */

    ///////////////////////
    // Add wheel in DEM
    ///////////////////////

    // Instantiate this wheel
    std::cout << "Making wheels..." << std::endl;
    DEMSim.SetFamilyPrescribedAngVel(100);
    DEMSim.SetFamilyPrescribedLinVel(100);
    std::vector<std::shared_ptr<DEMTracker>> trackers;
    // std::vector<std::shared_ptr<DEMClumpBatch>> DEM_Wheels;
    std::vector<std::shared_ptr<DEMMeshConnected>> DEM_Wheels;
    for (int i = 0; i < nW; i++) {
        // DEM_Wheels.push_back(DEMSim.AddClumps(wheel_template, make_float3(wheel_pos[i].x(), wheel_pos[i].y(),
        // wheel_pos[i].z())));
        // DEM_Wheels[i]->InformCentroidPrincipal(make_float3(0), make_float4(0.7071, 0, 0, 0.7071));
        DEM_Wheels.push_back(DEMSim.AddWavefrontMeshObject(wheel_obj_path, mat_type_wheel));
        // If this is one of the left wheels, mirror it. It's not a big difference, but we should do...
        // if (i == 0 || i == 2) {
        //     DEM_Wheels[i]->Mirror(make_float3(0,0,0), make_float3(0,1,0));
        // }

        DEM_Wheels[i]->SetFamily(100);
        DEM_Wheels[i]->SetMass(wheel_mass);
        DEM_Wheels[i]->SetMOI(wheel_MOI);
        trackers.push_back(DEMSim.Track(DEM_Wheels[i]));
    }
    DEMSim.DisableFamilyOutput(100);  // no need outputting wheels (if it's mesh, actually won't be outputted anyway)
    // DEMSim.DisableContactBetweenFamilies(100, 0);  // No wheel-ground contact while settling
    std::cout << "Total num of triangles in a wheel: " << DEM_Wheels[0]->GetNumTriangles() << std::endl;

    //////
    // Make ready for DEM simulation
    ///////
    // Some inspectors
    auto max_z_finder = DEMSim.CreateInspector("clump_max_z");
    auto min_z_finder = DEMSim.CreateInspector("clump_min_z");
    auto total_mass_finder = DEMSim.CreateInspector("clump_mass");
    auto partial_mass_finder = DEMSim.CreateInspector("clump_mass", "return (Z <= -0.41);");
    auto max_v_finder = DEMSim.CreateInspector("clump_max_absv");

    /*
    // Now add a plane to compress the `road'
    auto compressor = DEMSim.AddExternalObject();
    compressor->AddPlane(make_float3(0, 0, 0), make_float3(0, 0, -1), mat_type_terrain);
    compressor->SetFamily(90);
    DEMSim.SetFamilyFixed(90);
    DEMSim.DisableContactBetweenFamilies(90, 100);  // no contact between compressor and wheels
    auto compressor_tracker = DEMSim.Track(compressor);
    */

    float base_vel = 0.4;
    DEMSim.SetCoordSysOrigin(make_float3(world_x_size / 2., world_y_size / 2., world_y_size / 2.));
    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetGravitationalAcceleration(ChVec2Float(G));
    // If you want to use a large UpdateFreq then you have to expand spheres to ensure safety
    DEMSim.SetCDUpdateFreq(20);
    DEMSim.SetMaxVelocity(50.);
    DEMSim.SetInitBinSize(scales.at(2));
    DEMSim.SetIntegrator(TIME_INTEGRATOR::EXTENDED_TAYLOR);

    DEMSim.Initialize();
    for (const auto& tracker : trackers) {
        std::cout << "A tracker is tracking owner " << tracker->obj->ownerID << std::endl;
    }
    std::cout << "End initialization" << std::endl;

    ///////////////////////////////////////////
    // Compress the road first
    ///////////////////////////////////////////

    float time_end = 15.0;
    unsigned int fps = 10;
    unsigned int report_freq = 2000;
    unsigned int param_update_freq = 10000;
    unsigned int out_steps = (unsigned int)(1.0 / (fps * step_size));
    float frame_accu_thres = 1.0 / fps;
    unsigned int report_steps = (unsigned int)(1.0 / (report_freq * step_size));
    unsigned int param_update_steps = (unsigned int)(1.0 / (param_update_freq * step_size));

    // path out_dir = current_path();
    path out_dir = GetChronoOutputPath();
    create_directory(out_dir);
    out_dir += "/DEME";
    create_directory(out_dir);
    out_dir += "/Rover_on_GRC_incline_" + std::to_string(Slope_deg) + "deg";
    path rover_dir = out_dir / "./rover";
    create_directory(out_dir);
    create_directory(rover_dir);
    unsigned int currframe = 0;
    unsigned int curr_step = 0;

    /*
    step_size = 1e-6;
    DEMSim.SetInitTimeStep(step_size);
    DEMSim.UpdateSimParams();
    // Settle for a while...
    for (float t = 0; t < 1.5; t += step_size, curr_step++) {
        if (curr_step % out_steps == 0) {
            std::cout << "Frame: " << currframe << std::endl;
            DEMSim.ShowThreadCollaborationStats();
            char filename[200];
            sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe);
            DEMSim.WriteSphereFile(std::string(filename));
            SaveParaViewFiles(viper, rover_dir, currframe);
            currframe++;
        }
        DEMSim.DoDynamics(step_size);
    }

    // Start compressing
    DEMSim.DoDynamicsThenSync(0);
    step_size = 1e-6;
    DEMSim.SetInitTimeStep(step_size);
    DEMSim.UpdateSimParams();
    float matter_volume = void_ratio_finder->GetValue();
    std::cout << "Void ratio before simulation: " << (total_volume - matter_volume) / matter_volume << std::endl;

    double now_z = -0.38;
    compressor_tracker->SetPos(make_float3(0, 0, now_z));
    float compress_time = 0.3;
    double compressor_final_dist = (now_z > -0.42) ? now_z - (-0.42) : 0.0;
    std::cout << "Compressor is going to travel for " << compressor_final_dist << " meters" << std::endl;
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
    step_size = 1e-6;
    DEMSim.SetInitTimeStep(step_size);
    DEMSim.UpdateSimParams();

    DEMSim.EnableContactBetweenFamilies(100, 0);  // Re-enable wheel-ground contact

    matter_volume = void_ratio_finder->GetValue();
    std::cout << "Void ratio now: " << (total_volume - matter_volume) / matter_volume << std::endl;
    std::cout << "========================" << std::endl;
    */

    ///////////////////////////////////////////
    // Real simulation
    ///////////////////////////////////////////
    // Put the wheel in place, then let the wheel sink in initially
    float max_z = max_z_finder->GetValue();
    // wheel_tracker->SetPos(make_float3(init_x, 0, max_z + 0.03 + wheel_rad));

    float bulk_den_high = partial_mass_finder->GetValue() / ((-0.41 + 0.5) * world_x_size * world_y_size);
    float bulk_den_low = total_mass_finder->GetValue() / ((max_z + 0.5) * world_x_size * world_y_size);
    std::cout << "Bulk density high: " << bulk_den_high << std::endl;
    std::cout << "Bulk density low: " << bulk_den_low << std::endl;

    // Timers
    std::chrono::high_resolution_clock::time_point h_start, d_start;
    std::chrono::high_resolution_clock::time_point h_end, d_end;
    std::chrono::duration<double> h_total, d_total;
    h_total = std::chrono::duration<double>(0);
    d_total = std::chrono::duration<double>(0);

    std::vector<ChQuaternion<>> wheel_rot(4);
    std::vector<ChVector<>> wheel_vel(4);
    std::vector<ChVector<>> wheel_angVel(4);
    float max_v;
    int change_step = 0;
    float frame_accu = frame_accu_thres;

    // Find max z
    // float init_max_z =
    // 0.268923
    unsigned int chrono_update_freq = 10;
    for (float t = 0; t < time_end; t += step_size, curr_step++, frame_accu += step_size) {
        if (curr_step % chrono_update_freq == 0) {
            for (int i = 0; i < nW; i++) {
                wheel_pos[i] = Wheels[i]->GetFrame_REF_to_abs().GetPos();
                trackers[i]->SetPos(ChVec2Float(wheel_pos[i]));
                wheel_rot[i] = Wheels[i]->GetFrame_REF_to_abs().GetRot();
                trackers[i]->SetOriQ(ChQ2Float(wheel_rot[i]));
                wheel_vel[i] = Wheels[i]->GetFrame_REF_to_abs().GetPos_dt();
                trackers[i]->SetVel(ChVec2Float(wheel_vel[i]));
                wheel_angVel[i] = Wheels[i]->GetFrame_REF_to_abs().GetWvel_par();
                trackers[i]->SetAngVel(ChVec2Float(wheel_angVel[i]));
            }
        }
        if (curr_step % report_steps == 0) {
            // for (int i = 0; i < nW; i++) {
            //     std::cout << "Wheel " << i << " position: " << wheel_pos[i].x() << ", " << wheel_pos[i].y() << ", "
            //               << wheel_pos[i].z() << std::endl;

            //     std::cout << "Wheel " << i << " rotation: " << wheel_rot[i].e0() << ", " << wheel_rot[i].e1() << ", "
            //               << wheel_rot[i].e2() << ", " << wheel_rot[i].e3() << std::endl;

            //     std::cout << "Wheel " << i << " angVel: " << wheel_angVel[i].x() << ", " << wheel_angVel[i].y() << ",
            //     "
            //               << wheel_angVel[i].z() << std::endl;
            // }
            float3 body_pos = ChVec2Float(Body_1->GetFrame_REF_to_abs().GetPos());
            float3 body_vel = ChVec2Float(Body_1->GetFrame_REF_to_abs().GetPos_dt());
            float slip = 1.0 - body_vel.x / (w_r * wheel_rad);
            max_v = max_v_finder->GetValue();
            std::cout << "Current slope: " << Slope_deg << std::endl;
            std::cout << "Time is " << t << std::endl;
            std::cout << "X: " << body_pos.x << std::endl;
            std::cout << "V: " << body_vel.x << std::endl;
            std::cout << "Slip: " << slip << std::endl;
            std::cout << "Max vel in simulation is " << max_v << std::endl;
            std::cout << "========================" << std::endl;
        }

        if (frame_accu >= frame_accu_thres) {
            frame_accu = 0.;
            std::cout << "Frame: " << currframe << std::endl;
            DEMSim.ShowThreadCollaborationStats();
            char filename[200], meshfilename[200];
            sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe);
            // sprintf(meshfilename, "%s/DEMdemo_output_mesh_%04d.vtk", out_dir.c_str(), currframe);
            DEMSim.WriteSphereFile(std::string(filename));
            // DEMSim.WriteMeshFile(std::string(meshfilename));
            SaveParaViewFiles(viper, rover_dir, currframe);
            currframe++;
        }

        // Run DEM first
        d_start = std::chrono::high_resolution_clock::now();
        DEMSim.DoDynamics(step_size);
        d_end = std::chrono::high_resolution_clock::now();
        d_total += d_end - d_start;

        // Then feed force
        if (curr_step % chrono_update_freq == 0) {
            for (int i = 0; i < nW; i++) {
                float3 F = trackers[i]->ContactAcc();
                F *= wheel_mass;
                float3 tor = trackers[i]->ContactAngAccLocal();
                tor = wheel_MOI * tor;
                Wheels[i]->Empty_forces_accumulators();
                Wheels[i]->Accumulate_force(Float2ChVec(F), wheel_pos[i], false);
                Wheels[i]->Accumulate_torque(Float2ChVec(tor), true);  // torque in DEME is local
            }
            h_start = std::chrono::high_resolution_clock::now();
            sys.DoStepDynamics(step_size * chrono_update_freq);
            viper.Update();
            h_end = std::chrono::high_resolution_clock::now();
            h_total += h_end - h_start;
        }

        t += step_size;
        frame_accu += step_size;

        // if (curr_step % param_update_steps == 0 && t < 0.2) {
        // // if (t < 0.25) {
        //     DEMSim.DoDynamicsThenSync(0);
        //     max_v = max_v_finder->GetValue();
        //     float multiplier = max_v / base_vel;
        //     step_size = base_step_size / multiplier;
        //     DEMSim.SetInitTimeStep(step_size);
        //     // DEMSim.SetMaxVelocity(max_v * 1.2);
        //     DEMSim.UpdateSimParams();
        //     std::cout << "Max vel in simulation is " << max_v << std::endl;
        //     std::cout << "Step size in simulation is " << step_size << std::endl;
        // }

        // if (t > 1.5 && change_step == 0) {
        //     DEMSim.DoDynamicsThenSync(0);
        //     step_size = 1.5e-6;
        //     DEMSim.SetInitTimeStep(step_size);
        //     DEMSim.UpdateSimParams();
        //     change_step = 1;
        // }

        // else if (t > 3.0 && change_step == 1) {
        //     DEMSim.DoDynamicsThenSync(0);
        //     step_size = 2e-6;
        //     DEMSim.SetInitTimeStep(step_size);
        //     DEMSim.SetMaxVelocity(20.0);
        //     DEMSim.UpdateSimParams();
        //     change_step = 2;
        // }
        // else if (t > 2.0 && change_step == 2) {
        //     DEMSim.DoDynamicsThenSync(0);
        //     step_size = 3e-6;
        //     DEMSim.SetInitTimeStep(step_size);
        //     DEMSim.UpdateSimParams();
        //     change_step = 3;
        // }
        // else if (t > 3.0 && change_step == 3) {
        //     DEMSim.DoDynamicsThenSync(0);
        //     step_size = 5e-6;
        //     DEMSim.SetInitTimeStep(step_size);
        //     DEMSim.UpdateSimParams();
        //     change_step = 4;
        // }
    }
    std::cout << h_total.count() << " seconds spent on host" << std::endl;
    std::cout << d_total.count() << " seconds spent on device" << std::endl;

    DEMSim.ShowThreadCollaborationStats();
    DEMSim.ShowTimingStats();
    DEMSim.ShowAnomalies();

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
        std::string obj_path = "./Moon_rover_wheel.obj";
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

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
    float step_size = 2e-6;  // 2e-6; // 1e-6 for 15 deg and above, perhaps

    // Define the wheel geometry
    float wheel_rad = 0.25;
    float wheel_width = 0.25;
    float wheel_mass = 5.;
    float wheel_IYY = wheel_mass * wheel_rad * wheel_rad / 2;
    float wheel_IXX = (wheel_mass / 12) * (3 * wheel_rad * wheel_rad + wheel_width * wheel_width);
    float3 wheel_MOI = make_float3(wheel_IXX, wheel_IYY, wheel_IXX);

    std::string wheel_obj_path = GetDEMEDataFile("mesh/rover_wheels/viper_wheel_right.obj");
    // std::string wheel_obj_path = GetChronoDataFile("robot/viper/obj/viper_wheel.obj");
    // std::string wheel_obj_path = "./Moon_rover_wheel.obj";

    // Create a Chrono::Engine physical system
    float Slope_deg = 15;
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

    float x_offset = (Slope_deg > 19) ? 0.2 : 0.0;
    float z_offset = 0.085;
    viper.Initialize(ChFrame<>(ChVector<>(-0.9 + x_offset, -0.0, -0.135 + z_offset), QUNIT));

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

    // Family 1 is for fixed ground which does not participate the force calc.
    DEMSim.SetFamilyFixed(1);
    DEMSim.DisableContactBetweenFamilies(1, 1);
    DEMSim.DisableContactBetweenFamilies(1, 255);

    DEMSim.SetCollectAccRightAfterForceCalc(true);

    // E, nu, CoR, mu, Crr...
    float mu = 0.4;
    float mu_wheel = 0.8;
    float mu_wall = 1.;
    float CoR = 0.25;
    float E = 1e8;
    auto mat_type_wall = DEMSim.LoadMaterial({{"E", E}, {"nu", 0.3}, {"CoR", CoR}, {"mu", mu_wall}, {"Crr", 0.00}});
    auto mat_type_wheel = DEMSim.LoadMaterial({{"E", E}, {"nu", 0.3}, {"CoR", CoR}, {"mu", mu_wheel}, {"Crr", 0.00}});
    auto mat_type_terrain = DEMSim.LoadMaterial({{"E", E}, {"nu", 0.3}, {"CoR", CoR}, {"mu", mu}, {"Crr", 0.00}});
    DEMSim.SetMaterialPropertyPair("mu", mat_type_wheel, mat_type_terrain, mu_wheel);
    DEMSim.SetMaterialPropertyPair("mu", mat_type_wall, mat_type_terrain, mu_wall);

    // Define the simulation world
    double world_x_size = 4.0;
    double world_y_size = 2.0;
    DEMSim.InstructBoxDomainDimension(world_x_size, world_y_size, world_y_size);
    float bottom = -0.5;
    DEMSim.AddBCPlane(make_float3(0, 0, bottom), make_float3(0, 0, 1), mat_type_wall);
    DEMSim.AddBCPlane(make_float3(0, world_y_size / 2, 0), make_float3(0, -1, 0), mat_type_wall);
    DEMSim.AddBCPlane(make_float3(0, -world_y_size / 2, 0), make_float3(0, 1, 0), mat_type_wall);
    // X-dir bounding planes
    DEMSim.AddBCPlane(make_float3(-world_x_size / 2., 0, 0), make_float3(1, 0, 0), mat_type_wall);
    DEMSim.AddBCPlane(make_float3(world_x_size / 2., 0, 0), make_float3(-1, 0, 0), mat_type_wall);

    // Define the terrain particle templates
    // Calculate its mass and MOI
    float mass1 = 2.6e3 * 4.2520508;
    float3 MOI1 = make_float3(1.6850426, 1.6375114, 2.1187753) * 2.6e3;
    float mass2 = 2.6e3 * 2.1670011;
    float3 MOI2 = make_float3(0.57402126, 0.60616378, 0.92890173) * 2.6e3;
    // Scale the template we just created
    std::vector<double> scales = {0.0014, 0.00075833, 0.00044, 0.0003, 0.0002, 0.00018333, 0.00017};
    std::for_each(scales.begin(), scales.end(), [](double& r) { r *= 10.; });
    // Then load it to system
    std::shared_ptr<DEMClumpTemplate> my_template2 =
        DEMSim.LoadClumpType(mass2, MOI2, GetDEMEDataFile("clumps/triangular_flat_6comp.csv"), mat_type_terrain);
    std::shared_ptr<DEMClumpTemplate> my_template1 =
        DEMSim.LoadClumpType(mass1, MOI1, GetDEMEDataFile("clumps/triangular_flat.csv"), mat_type_terrain);
    std::vector<std::shared_ptr<DEMClumpTemplate>> ground_particle_templates = {my_template2,
                                                                                DEMSim.Duplicate(my_template2),
                                                                                my_template1,
                                                                                DEMSim.Duplicate(my_template1),
                                                                                DEMSim.Duplicate(my_template1),
                                                                                DEMSim.Duplicate(my_template1),
                                                                                DEMSim.Duplicate(my_template1)};
    // Now scale those templates
    for (int i = 0; i < scales.size(); i++) {
        std::shared_ptr<DEMClumpTemplate>& my_template = ground_particle_templates.at(i);
        // Note the mass and MOI are also scaled in the process, automatically. But if you are not happy with this, you
        // can always manually change mass and MOI afterwards.
        my_template->Scale(scales.at(i));
        // Give these templates names, 0000, 0001 etc.
        char t_name[20];
        sprintf(t_name, "%04d", i);
        my_template->AssignName(std::string(t_name));
    }

    // Now we load clump locations from a checkpointed file
    {
        std::cout << "Making terrain..." << std::endl;
        std::unordered_map<std::string, std::vector<float3>> clump_xyz;
        std::unordered_map<std::string, std::vector<float4>> clump_quaternion;
        try {
            clump_xyz = DEMSim.ReadClumpXyzFromCsv("./GRC_20e6.csv");
            clump_quaternion = DEMSim.ReadClumpQuatFromCsv("./GRC_20e6.csv");
        } catch (...) {
            std::cout << "You will need to finish the GRCPrep demos first to obtain the checkpoint file GRC_20e6.csv, "
                         "in order to run this demo. \nThat is a 4m by 2m GRC-1 terrain patch that is around 15cm "
                         "thick.\nIf you don't have access to it, you can go to the forum "
                         "(https://groups.google.com/g/projectchrono) to ask the authors for it."
                      << std::endl;
            return 1;
        }
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

        // Maybe we need to make the terrain a bit thicker...
        // I found when I created the initial terrain it wasn't thick enough for the rover simulation at steeper slopes.
        // I was too lazy to regenerate everything, so I just make and `partial' copy of the existing terrain, add it on
        // top of the existing terrain, to make it thicker. The following section of the code is doing that. It may have
        // been done in a more or less confusing way of coding, but it doesn't matter; if it is too arcane for you, just
        // get the idea and move on to reading the rest of the script.
        {
            float up_dist, remove_pos;
            up_dist = 0.1;
            remove_pos = 0.43;

            std::vector<float> x_shift_dist = {0};
            std::vector<float> y_shift_dist = {0};
            std::vector<float> z_shift_dist = {up_dist};
            // Add some patches of such graular bed
            for (float x_shift : x_shift_dist) {
                for (float y_shift : y_shift_dist) {
                    for (float z_shift : z_shift_dist) {
                        std::vector<float3> my_xyz = in_xyz;
                        std::vector<float4> my_quat = in_quat;
                        std::vector<std::shared_ptr<DEMClumpTemplate>> my_types = in_types;
                        std::vector<notStupidBool_t> elem_to_remove(in_xyz.size(), 0);
                        for (size_t i = 0; i < in_xyz.size(); i++) {
                            if (in_xyz.at(i).z < -remove_pos)
                                elem_to_remove.at(i) = 1;
                        }
                        my_xyz.erase(std::remove_if(my_xyz.begin(), my_xyz.end(),
                                                    [&elem_to_remove, &my_xyz](const float3& i) {
                                                        return elem_to_remove.at(&i - my_xyz.data());
                                                    }),
                                     my_xyz.end());
                        my_quat.erase(std::remove_if(my_quat.begin(), my_quat.end(),
                                                     [&elem_to_remove, &my_quat](const float4& i) {
                                                         return elem_to_remove.at(&i - my_quat.data());
                                                     }),
                                      my_quat.end());
                        my_types.erase(std::remove_if(my_types.begin(), my_types.end(),
                                                      [&elem_to_remove, &my_types](const auto& i) {
                                                          return elem_to_remove.at(&i - my_types.data());
                                                      }),
                                       my_types.end());
                        DEMClumpBatch another_batch(my_xyz.size());
                        std::for_each(my_xyz.begin(), my_xyz.end(), [x_shift, y_shift, z_shift](float3& xyz) {
                            xyz.x += x_shift;
                            xyz.y += y_shift;
                            xyz.z += z_shift;
                        });
                        another_batch.SetTypes(my_types);
                        another_batch.SetPos(my_xyz);
                        another_batch.SetOriQ(my_quat);
                        DEMSim.AddClumps(another_batch);
                    }
                }
            }
        }
    }

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

    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetGravitationalAcceleration(ChVec2Float(G));
    // If you want to use a large UpdateFreq then you have to expand spheres to ensure safety
    DEMSim.SetCDUpdateFreq(30);
    // DEMSim.SetInitBinSize(scales.at(1));
    DEMSim.SetInitBinNumTarget(5e6);

    DEMSim.SetExpandSafetyAdder(0.2);
    DEMSim.SetExpandSafetyMultiplier(1.);
    DEMSim.SetErrorOutVelocity(100.);

    DEMSim.Initialize();
    for (const auto& tracker : trackers) {
        std::cout << "A tracker is tracking owner " << tracker->obj->ownerID << std::endl;
    }
    std::cout << "Time step size is " << step_size << std::endl;
    std::cout << "End initialization" << std::endl;

    ///////////////////////////////////////////
    // Compress the road first
    ///////////////////////////////////////////

    float time_end = 15.0;
    unsigned int fps = 10;
    // unsigned int move_box_ps = 1;
    unsigned int report_freq = 2000;
    unsigned int param_update_freq = 10000;
    unsigned int out_steps = (unsigned int)(1.0 / (fps * step_size));
    float frame_accu_thres = 1.0 / fps;
    unsigned int report_steps = (unsigned int)(1.0 / (report_freq * step_size));
    unsigned int param_update_steps = (unsigned int)(1.0 / (param_update_freq * step_size));

    // path out_dir = current_path();
    path out_dir = GetChronoOutputPath();
    create_directory(out_dir);
    // out_dir += "/DEME_FullDomain";
    out_dir += "/DEME";
    create_directory(out_dir);
    out_dir += "/Rover_on_GRC_incline_" + std::to_string(Slope_deg);
    path rover_dir = out_dir / "./rover";
    create_directory(out_dir);
    create_directory(rover_dir);
    unsigned int currframe = 0;
    unsigned int curr_step = 0;

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
    std::vector<ChVector<>> wheel_vel(4);
    std::vector<ChVector<>> wheel_angVel(4);
    float max_v;
    int change_step = 0;
    float frame_accu = frame_accu_thres;

    // Put the wheels somewhere that won't affect simulation
    {
        trackers[0]->SetPos(make_float3(1, 0.5, 0.75));
        trackers[1]->SetPos(make_float3(1, -0.5, 0.75));
        trackers[2]->SetPos(make_float3(-1, 0.5, 0.75));
        trackers[3]->SetPos(make_float3(-1, -0.5, 0.75));
    }

    // Settle first, then put the wheel in place, then let the wheel sink in initially
    for (float t = 0; t < 0.2; t += frame_accu_thres) {
        std::cout << "Num contacts: " << DEMSim.GetNumContacts() << std::endl;
        DEMSim.ShowThreadCollaborationStats();
        // std::cout << "Frame: " << currframe << std::endl;
        // char filename[200];
        // sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe);
        // DEMSim.WriteSphereFile(std::string(filename));
        // SaveParaViewFiles(viper, rover_dir, currframe);
        // currframe++;
        DEMSim.ShowTimingStats();
        DEMSim.DoDynamicsThenSync(frame_accu_thres);
    }

    float max_z = max_z_finder->GetValue();
    // wheel_tracker->SetPos(make_float3(init_x, 0, max_z + 0.03 + wheel_rad));

    float bulk_den_high = partial_mass_finder->GetValue() / ((-0.41 + 0.5) * world_x_size * world_y_size);
    float bulk_den_low = total_mass_finder->GetValue() / ((max_z + 0.5) * world_x_size * world_y_size);
    std::cout << "Bulk density high: " << bulk_den_high << std::endl;
    std::cout << "Bulk density low: " << bulk_den_low << std::endl;
    std::cout << "Max Z initially: " << max_z << std::endl;

    std::cout << "Num contacts: " << DEMSim.GetNumContacts() << std::endl;

    unsigned int chrono_update_freq = 10;
    // Active box sizes
    float box_halfsize_x = 0.5;
    float box_halfsize_y = 0.25;

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

        if (frame_accu >= frame_accu_thres) {
            frame_accu = 0.;
            std::cout << "Frame: " << currframe << std::endl;
            std::cout << "Num contacts: " << DEMSim.GetNumContacts() << std::endl;
            std::cout << h_total.count() << " seconds spent on host" << std::endl;
            std::cout << d_total.count() << " seconds spent on device" << std::endl;
            DEMSim.ShowThreadCollaborationStats();
            char filename[200];
            sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe);
            DEMSim.WriteSphereFile(std::string(filename));
            SaveParaViewFiles(viper, rover_dir, currframe);
            currframe++;
            DEMSim.ShowTimingStats();

            // Move the active box
            // if (Slope_deg < 19. && t > 0.2) {
            if (t > 0.2) {
                DEMSim.DoDynamicsThenSync(0.);
                DEMSim.ChangeClumpFamily(1);
                size_t num_changed = 0;
                for (int i = 0; i < nW; i++) {
                    wheel_pos[i] = Wheels[i]->GetFrame_REF_to_abs().GetPos();
                    float3 pos = ChVec2Float(wheel_pos[i]);
                    std::pair<float, float> Xrange =
                        std::pair<float, float>(pos.x - box_halfsize_x, pos.x + box_halfsize_x);
                    std::pair<float, float> Yrange =
                        std::pair<float, float>(pos.y - box_halfsize_y, pos.y + box_halfsize_y);
                    num_changed += DEMSim.ChangeClumpFamily(0, Xrange, Yrange);
                }
                std::cout << num_changed << " particles changed family number." << std::endl;
            }
        }

        // Run DEM
        d_start = std::chrono::high_resolution_clock::now();
        DEMSim.DoDynamics(step_size);
        d_end = std::chrono::high_resolution_clock::now();
        d_total += d_end - d_start;

        // Feed force
        if (curr_step % chrono_update_freq == 0) {
            for (int i = 0; i < nW; i++) {
                float3 F = trackers[i]->ContactAcc();
                F *= wheel_mass;
                float3 tor = trackers[i]->ContactAngAccLocal();
                tor = wheel_MOI * tor;
                wheel_pos[i] = Wheels[i]->GetFrame_REF_to_abs().GetPos();
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

        if (curr_step % report_steps == 0) {
            // float wheel_accu_vel = 0., wheel_accu_x = 0.;
            // for (int i = 0; i < nW; i++) {
            //     wheel_accu_vel += Wheels[i]->GetFrame_REF_to_abs().GetPos_dt().x();
            //     wheel_accu_x += Wheels[i]->GetFrame_REF_to_abs().GetPos().x();
            // }
            float rover_vel = Body_1->GetFrame_REF_to_abs().GetPos_dt().x();
            // float wheel_vel = wheel_accu_vel / nW;
            float rover_pos = Body_1->GetFrame_REF_to_abs().GetPos().x();
            // float wheel_pos = wheel_accu_x / nW;
            float slip = 1.0 - rover_vel / (w_r * wheel_rad);
            max_v = max_v_finder->GetValue();
            std::cout << "Current slope: " << Slope_deg << std::endl;
            std::cout << "Time is " << t << std::endl;
            std::cout << "X: " << rover_pos << std::endl;
            // std::cout << "Wheel X: " << wheel_pos << std::endl;
            std::cout << "V: " << rover_vel << std::endl;
            // std::cout << "Wheel V: " << wheel_vel << std::endl;
            std::cout << "Slip: " << slip << std::endl;
            std::cout << "Max vel in simulation is " << max_v << std::endl;
            std::cout << "========================" << std::endl;
            if (rover_pos > 1.) {
                std::cout << "This is far enough, stopping the simulation..." << std::endl;
                std::cout << "========================" << std::endl;
                DEMSim.DoDynamicsThenSync(0.);
                break;
            }
        }

        // if (curr_step % param_update_steps == 0 && t < 0.2) {
        // // if (t < 0.25) {
        //     DEMSim.DoDynamicsThenSync(0);
        //     max_v = max_v_finder->GetValue();
        //     float multiplier = max_v / base_vel;
        //     step_size = base_step_size / multiplier;
        //     DEMSim.SetInitTimeStep(step_size);
        //     // DEMSim.SetMaxVelocity(max_v * 1.2);
        //     DEMSim.UpdateStepSize();
        //     std::cout << "Max vel in simulation is " << max_v << std::endl;
        //     std::cout << "Step size in simulation is " << step_size << std::endl;
        // }

        // if (t > 1.5 && change_step == 0) {
        //     DEMSim.DoDynamicsThenSync(0);
        //     step_size = 1.5e-6;
        //     DEMSim.SetInitTimeStep(step_size);
        //     DEMSim.UpdateStepSize();
        //     change_step = 1;
        // }

        // else if (t > 3.0 && change_step == 1) {
        //     DEMSim.DoDynamicsThenSync(0);
        //     step_size = 2e-6;
        //     DEMSim.SetInitTimeStep(step_size);
        //     DEMSim.SetMaxVelocity(20.0);
        //     DEMSim.UpdateStepSize();
        //     change_step = 2;
        // }
        // else if (t > 2.0 && change_step == 2) {
        //     DEMSim.DoDynamicsThenSync(0);
        //     step_size = 3e-6;
        //     DEMSim.SetInitTimeStep(step_size);
        //     DEMSim.UpdateStepSize();
        //     change_step = 3;
        // }
        // else if (t > 3.0 && change_step == 3) {
        //     DEMSim.DoDynamicsThenSync(0);
        //     step_size = 5e-6;
        //     DEMSim.SetInitTimeStep(step_size);
        //     DEMSim.UpdateStepSize();
        //     change_step = 4;
        // }
    }
    std::cout << "Finishing up..." << std::endl;
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
        // std::string obj_path = "./Moon_rover_wheel.obj";
        std::string obj_path = GetChronoDataFile("robot/viper/obj/viper_wheel.obj");
        double scale_ratio = 1.0;
        mmesh->LoadWavefrontMesh(obj_path, false, true);
        mmesh->Transform(ChVector<>(0, 0, 0), ChMatrix33<>(scale_ratio));  // scale to a different size
        mmesh->RepairDuplicateVertexes(1e-9);                              // if meshes are not watertight
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
        mmesh->Transform(body_pos, ChMatrix33<>(body_rot));  // rotate the mesh based on the orientation of body

        // filename = rover_dir / ("./upperRod_" + std::to_string(i + 1) + "_" + std::string(f_name) + ".obj");
        // std::vector<geometry::ChTriangleMeshConnected> meshes = {*mmesh};
        // geometry::ChTriangleMeshConnected::WriteWavefront(filename.string(), meshes);
        meshes.push_back(*mmesh);
    }

    geometry::ChTriangleMeshConnected::WriteWavefront(filename.string(), meshes);
}

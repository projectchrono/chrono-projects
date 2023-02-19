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
// Demo to show a cone-drop simulation on complex DEM terrain
//
// =============================================================================

#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#include <DEM/API.h>
//#include <core/ApiVersion.h>
//#include <core/utils/DEMEPaths.hpp>
#include <DEM/HostSideHelpers.hpp>
#include <DEM/utils/Samplers.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <map>
#include <random>

using namespace deme;

using namespace chrono;
using namespace chrono::geometry;

const double math_PI = 3.14159;

int main() {
    SetChronoDataPath(CHRONO_DATA_DIR);
    SetDEMEDataPath(DEME_DATA_DIR);

    DEMSolver DEMSim;
    DEMSim.SetVerbosity(INFO);
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetOutputContent(OUTPUT_CONTENT::ABSV);
    DEMSim.SetMeshOutputFormat(MESH_FORMAT::VTK);
    DEMSim.SetContactOutputContent(OWNER | FORCE | POINT);

    // E, nu, CoR, mu, Crr...
    auto mat_type_cone = DEMSim.LoadMaterial({{"E", 1e9}, {"nu", 0.3}, {"CoR", 0.3}, {"mu", 0.5}, {"Crr", 0.00}});
    auto mat_type_terrain = DEMSim.LoadMaterial({{"E", 1e9}, {"nu", 0.3}, {"CoR", 0.3}, {"mu", 0.5}, {"Crr", 0.00}});

    float step_size = 1e-6;
    double world_size = 2;
    double soil_bin_diameter = 0.584;
    double cone_surf_area = 323e-6;
    double cone_diameter = std::sqrt(cone_surf_area / math_PI) * 2;
    double starting_height = -0.1;
    DEMSim.InstructBoxDomainDimension(world_size, world_size, world_size);
    // No need to add simulation `world' boundaries, b/c we'll add a cylinderical container manually
    DEMSim.InstructBoxDomainBoundingBC("none", mat_type_terrain);
    // Now add a cylinderical boundary along with a bottom plane
    double bottom = -0.5;
    auto walls = DEMSim.AddExternalObject();
    walls->AddCylinder(make_float3(0), make_float3(0, 0, 1), soil_bin_diameter / 2., mat_type_terrain, 0);
    walls->AddPlane(make_float3(0, 0, bottom), make_float3(0, 0, 1), mat_type_terrain);

    // Define the GRC terrain particle templates
    DEMClumpTemplate shape_template;
    shape_template.ReadComponentFromFile(GetDEMEDataFile("clumps/triangular_flat.csv"));
    // Calculate its mass and MOI
    float terrain_density = 2.6e3;
    float mass = terrain_density * 5.5886717;  // in kg or g
    float3 MOI = make_float3(1.8327927, 2.1580013, 0.77010059) * (double)2.6e3;
    double clump_vol = 5.5886717;
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
        this_template.SetVolume(clump_vol * scaling * scaling * scaling);
        ground_particle_templates.push_back(DEMSim.LoadClumpType(this_template));
    }

    // Now we load clump locations from a checkpointed file
    {
        std::cout << "Making terrain..." << std::endl;
        auto clump_xyz = DEMSim.ReadClumpXyzFromCsv("./GRC_3e6.csv");
        auto clump_quaternion = DEMSim.ReadClumpQuatFromCsv("./GRC_3e6.csv");
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
        // Now, we don't need all particles loaded... we just need a cylinderical portion out of it, to fill the soil
        // bin Remove the particles that are outside a cylinderical region
        std::vector<notStupidBool_t> elem_to_remove(in_xyz.size(), 0);
        for (size_t i = 0; i < in_xyz.size(); i++) {
            if (std::pow(in_xyz.at(i).x, 2) + std::pow(in_xyz.at(i).y, 2) >= std::pow(soil_bin_diameter / 2. - 0.02, 2))
                elem_to_remove.at(i) = 1;
        }
        in_xyz.erase(std::remove_if(
                         in_xyz.begin(), in_xyz.end(),
                         [&elem_to_remove, &in_xyz](const float3& i) { return elem_to_remove.at(&i - in_xyz.data()); }),
                     in_xyz.end());
        in_quat.erase(std::remove_if(in_quat.begin(), in_quat.end(),
                                     [&elem_to_remove, &in_quat](const float4& i) {
                                         return elem_to_remove.at(&i - in_quat.data());
                                     }),
                      in_quat.end());
        in_types.erase(std::remove_if(in_types.begin(), in_types.end(),
                                      [&elem_to_remove, &in_types](const auto& i) {
                                          return elem_to_remove.at(&i - in_types.data());
                                      }),
                       in_types.end());
        DEMClumpBatch base_batch(in_xyz.size());
        base_batch.SetTypes(in_types);
        base_batch.SetPos(in_xyz);
        base_batch.SetOriQ(in_quat);

        DEMSim.AddClumps(base_batch);

        // This batch is about 10cm thick... let's add another 2 batches, so we have something like 30cm
        float shift_dist = 0.13;
        for (int i = 0; i < 2; i++) {
            std::for_each(in_xyz.begin(), in_xyz.end(), [shift_dist](float3& xyz) { xyz.z += shift_dist; });
            DEMClumpBatch another_batch = base_batch;
            another_batch.SetPos(in_xyz);
            DEMSim.AddClumps(another_batch);
        }
    }

    // Load in the cone used for this penetration test
    auto cone_tip = DEMSim.AddWavefrontMeshObject(GetDEMEDataFile("mesh/cone.obj"), mat_type_cone);
    auto cone_body = DEMSim.AddWavefrontMeshObject(GetDEMEDataFile("mesh/cyl_r1_h2.obj"), mat_type_cone);
    std::cout << "Total num of triangles: " << cone_tip->GetNumTriangles() + cone_body->GetNumTriangles() << std::endl;

    // The initial cone mesh has base radius 1, and height 1. Let's stretch it a bit so it has a 60deg tip, instead of
    // 90deg.
    float tip_height = std::sqrt(3.);
    cone_tip->Scale(make_float3(1, 1, tip_height));
    // Then set mass properties
    float cone_mass = 7.8e3 * tip_height / 3 * math_PI;
    cone_tip->SetMass(cone_mass);
    // You can checkout https://en.wikipedia.org/wiki/List_of_moments_of_inertia
    cone_tip->SetMOI(make_float3(cone_mass * (3. / 20. + 3. / 80. * tip_height * tip_height),
                                 cone_mass * (3. / 20. + 3. / 80. * tip_height * tip_height), 3 * cone_mass / 10));
    // This cone mesh has its tip at the origin. And, float4 quaternion pattern is (x, y, z, w).
    cone_tip->InformCentroidPrincipal(make_float3(0, 0, 3. / 4. * tip_height), make_float4(0, 0, 0, 1));
    // Note the scale method will scale mass and MOI automatically. But this only goes for the case you scale xyz all
    // together; otherwise, the MOI scaling will not be accurate and you should manually reset them.
    cone_tip->Scale(cone_diameter / 2);
    // Note that position of objects is always the location of their centroid
    cone_tip->SetInitPos(make_float3(0, 0, starting_height));
    cone_tip->SetFamily(1);
    // The tip location, used to measure penetration length
    double tip_z = -cone_diameter / 2 * 3 / 4 * tip_height + starting_height;

    // The define the body that is connected to the tip
    float body_mass = 7.8e3 * math_PI;
    cone_body->SetMass(body_mass);
    cone_body->SetMOI(make_float3(body_mass * 7 / 12, body_mass * 7 / 12, body_mass / 2));
    // This cyl mesh (h = 2m, r = 1m) has its center at the origin. So the following call actually has no effect...
    cone_body->InformCentroidPrincipal(make_float3(0, 0, 0), make_float4(0, 0, 0, 1));
    cone_body->Scale(make_float3(cone_diameter / 2, cone_diameter / 2, 0.5));
    // Its initial position should be right above the cone tip...
    cone_body->SetInitPos(make_float3(0, 0, 0.5 + (cone_diameter / 2 / 4 * tip_height) + starting_height));
    cone_body->SetFamily(1);

    // Track the cone_tip
    auto tip_tracker = DEMSim.Track(cone_tip);

    // In fact, because the cone's motion is completely pre-determined, we can just prescribe family 1
    DEMSim.SetFamilyPrescribedLinVel(1, "0", "0", "-0.025");

    // Some inspectors
    auto max_z_finder = DEMSim.CreateInspector("clump_max_z");
    auto void_ratio_finder =
        DEMSim.CreateInspector("clump_volume", "return (X * X + Y * Y <= 0.25 * 0.25) && (Z <= -0.3);");
    float total_volume = 0.2 * math_PI * (0.25 * 0.25);

    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -9.81));
    DEMSim.SetCDUpdateFreq(20);
    DEMSim.SetMaxVelocity(40.);
    DEMSim.SetInitBinSize(2 * scales.at(2));
    DEMSim.Initialize();

    std::filesystem::path out_dir = GetChronoOutputPath();
    std::filesystem::create_directory(out_dir);
    out_dir += "/DEME";
    std::filesystem::create_directory(out_dir);
    out_dir += "/Cone_Penetration";
    std::filesystem::path mesh_dir = out_dir / "./cone";
    std::filesystem::create_directory(out_dir);
    std::filesystem::create_directory(mesh_dir);
    unsigned int currframe = 0;
    // unsigned int curr_step = 0;

    float sim_end = 15.0;
    unsigned int fps = 20;
    float frame_time = 1.0 / fps;
    std::cout << "Output at " << fps << " FPS" << std::endl;

    int step_size_marker = 0;
    double tip_z_when_first_hit;
    bool hit_terrain = false;
    for (float t = 0; t < sim_end; t += frame_time) {
        // if (step_size_marker == 0 && t > 0.5) {
        //     step_size = 2e-6;
        //     DEMSim.SetInitTimeStep(step_size);
        //     DEMSim.UpdateStepSize();
        //     step_size_marker = 1;
        // }
        std::cout << "Frame: " << currframe << std::endl;
        char filename[200], meshfilename[200], cnt_filename[200];
        sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe);
        sprintf(meshfilename, "%s/DEMdemo_mesh_%04d.vtk", mesh_dir.c_str(), currframe);
        sprintf(cnt_filename, "%s/Contact_pairs_%04d.csv", out_dir.c_str(), currframe);
        DEMSim.WriteSphereFile(std::string(filename));
        DEMSim.WriteMeshFile(std::string(meshfilename));
        DEMSim.WriteContactFile(std::string(cnt_filename));
        currframe++;

        float matter_volume = void_ratio_finder->GetValue();
        float void_ratio = (total_volume - matter_volume) / matter_volume;
        float bulk_density = matter_volume / total_volume * terrain_density;
        // float terrain_max_z = max_z_finder->GetValue();
        float3 forces = tip_tracker->ContactAcc();
        // Note cone_mass is not the true mass, b/c we scaled the the cone tip!
        forces *= cone_tip->mass;
        float pressure = std::abs(forces.z) / cone_surf_area;
        if (pressure > 1e-8 && !hit_terrain) {
            hit_terrain = true;
            tip_z_when_first_hit = tip_z;
        }
        float penetration = (hit_terrain) ? tip_z_when_first_hit - tip_z : 0;
        std::cout << "Time: " << t << std::endl;
        std::cout << "Void ratio: " << void_ratio << std::endl;
        std::cout << "Bulk density: " << bulk_density << std::endl;
        std::cout << "Penetration: " << penetration << std::endl;
        std::cout << "Force on cone: " << forces.x << ", " << forces.y << ", " << forces.z << std::endl;
        std::cout << "Pressure: " << pressure << std::endl;

        DEMSim.DoDynamicsThenSync(frame_time);
        DEMSim.ShowThreadCollaborationStats();

        tip_z -= 0.025 * frame_time;
    }

    std::cout << "ConeDrop Cosim exiting..." << std::endl;
    return 0;
}
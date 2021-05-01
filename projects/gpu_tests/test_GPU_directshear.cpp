// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2019 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Nic Olsen
// =============================================================================
// Chrono::Granular simulation performing a standard direct shear test.
// Material is settled in a rectangular mesh box, then compressed by a top
// plate. The test then measures the shear stress caused by moving the top half
// of the box at a constant velocity.
// =============================================================================

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

#include "chrono/core/ChGlobal.h"
#include "chrono/physics/ChBody.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/utils/ChUtilsSamplers.h"

#include "chrono_gpu/physics/ChSystemGpu.h"
#include "chrono_gpu/utils/ChGpuJsonParser.h"
#include "chrono_gpu/ChGpuData.h"

#include "chrono_thirdparty/filesystem/path.h"

using namespace chrono;
using namespace chrono::gpu;

#define FAM_ENTRIES_POS 7
#define FAM_ENTRIES_VEL 6
#define FAM_ENTRIES_FORCE 6

// Normal stress values for four tests (Pa)
double normal_stresses[] = {3.1e3, 6.4e3, 12.5e3, 24.2e3};
double plate_mass;

double shear_velocity_original = 0.1;  // 1 mm/s
double shear_displacement = 1;         // X displacement at which the test ends
double shear_velocity_inflation = 10;  // Multiplier on shear_velocity speed for shorter simulation
double shear_velocity = shear_velocity_inflation * shear_velocity_original;

double box_xy = 6;
double box_r = box_xy / 2;

// TODO tune these values
double time_settle = 0.4;
double time_compress = 1;
double time_shear = shear_displacement / shear_velocity;

// Indices of each object
const size_t bottom_i = 0;
const size_t top_i = 1;
const size_t plate_i = 2;

double fill_top;

// expected number of args for param sweep
constexpr int num_args_full = 3;


void ShowUsage(std::string name) {
    std::cout << "usage: " + name + " <json_file> <normal_stress_index>" << std::endl;
    std::cout << "must have either 1 or " << num_args_full - 1 << " arguments" << std::endl;

}

void SetupGranSystem(ChSystemGpuMesh& gpu_sys, ChGpuSimulationParameters& params) {

    gpu_sys.SetKn_SPH2SPH(params.normalStiffS2S);
    gpu_sys.SetKn_SPH2WALL(params.normalStiffS2W);
    gpu_sys.SetKn_SPH2MESH(params.normalStiffS2M);

    gpu_sys.SetKt_SPH2SPH(params.tangentStiffS2S);
    gpu_sys.SetKt_SPH2WALL(params.tangentStiffS2W);
    gpu_sys.SetKt_SPH2MESH(params.tangentStiffS2M);

    gpu_sys.SetGn_SPH2SPH(params.normalDampS2S);
    gpu_sys.SetGn_SPH2WALL(params.normalDampS2W);
    gpu_sys.SetGn_SPH2MESH(params.normalDampS2M);

    gpu_sys.SetGt_SPH2SPH(params.tangentDampS2S);
    gpu_sys.SetGt_SPH2WALL(params.tangentDampS2W);
    gpu_sys.SetGt_SPH2MESH(params.tangentDampS2M);

    gpu_sys.SetCohesionRatio(params.cohesion_ratio);
    gpu_sys.SetAdhesionRatio_SPH2MESH(params.adhesion_ratio_s2m);
    gpu_sys.SetAdhesionRatio_SPH2WALL(params.adhesion_ratio_s2w);

    gpu_sys.SetGravitationalAcceleration(ChVector<float>(params.grav_X, params.grav_Y, params.grav_Z));

    gpu_sys.SetFrictionMode(chrono::gpu::CHGPU_FRICTION_MODE::SINGLE_STEP);
    float mu = 0.5;
    gpu_sys.SetStaticFrictionCoeff_SPH2SPH(mu);
    gpu_sys.SetStaticFrictionCoeff_SPH2WALL(mu);

    gpu_sys.SetParticleOutputMode(CHGPU_OUTPUT_MODE::CSV);

    gpu_sys.SetTimeIntegrator(CHGPU_TIME_INTEGRATOR::FORWARD_EULER);
    gpu_sys.SetFixedStepSize(params.step_size);
    gpu_sys.SetBDFixed(true);

    double epsilon = 0.02 * params.sphere_radius;
    double spacing = 2 * params.sphere_radius + epsilon;

    std::vector<ChVector<float>> body_points;

    // utils::HCPSampler<float> sampler(spacing);
    utils::PDSampler<float> sampler(spacing);
    double fill_bottom = -box_r + spacing;
    fill_top = params.box_Z / 2 - spacing;  // TODO tune to roughly make a cube of material (6cm tall)

    ChVector<> hdims(box_r - params.sphere_radius - epsilon, box_r - params.sphere_radius - epsilon, 0);

    for (double z = fill_bottom; z < fill_top; z += spacing) {
        ChVector<> center(0, 0, z);
        auto points = sampler.SampleBox(center, hdims);
        body_points.insert(body_points.end(), points.begin(), points.end());
    }

    std::cout << "Created " << body_points.size() << " spheres" << std::endl;

    gpu_sys.SetParticles(body_points);

    // Mesh values
    std::vector<string> mesh_filenames;
    // TODO dull the corners and fix nans
    mesh_filenames.push_back(std::string(gpu::GetDataFile("meshes/directshear/shear_bottom.obj")));
    mesh_filenames.push_back(std::string(gpu::GetDataFile("meshes/directshear/shear_top.obj")));
    mesh_filenames.push_back(std::string(gpu::GetDataFile("meshes/directshear/downward_square.obj")));

    ChMatrix33<float> scale(ChVector<float>(box_r, box_r, box_r));
    std::vector<ChMatrix33<float>> mesh_rotscales = {scale, scale, scale};
    std::vector<ChVector<float>> mesh_translations = {ChVector<float>(0, 0, 0), ChVector<float>(0, 0, 0),
                                                      ChVector<float>(0, 0, 0)};
    std::vector<float> mesh_masses = {1000, 1000, (float)plate_mass};

    gpu_sys.AddMeshes(mesh_filenames, mesh_translations, mesh_rotscales, mesh_masses);
}

void SetInitialMeshes(ChSystemGpuMesh& gpu_sys, const std::shared_ptr<ChBody> plate) {
    // initial positions and velocity
    ChVector<float> mesh_pos(0, 0, 0);
    ChQuaternion<float> mesh_rot(1, 0, 0, 0);
    ChVector<float> mesh_lin_vel(0, 0, 0);
    ChVector<float> mesh_ang_vel(0, 0 , 0);
    
    // Bottom bin
    gpu_sys.ApplyMeshMotion(bottom_i, mesh_pos, mesh_rot, mesh_lin_vel, mesh_ang_vel);

    // Top bin
    gpu_sys.ApplyMeshMotion(top_i, mesh_pos, mesh_rot, mesh_lin_vel, mesh_ang_vel);

    // Plate
    ChVector<float> plate_pos(0, 0, (float)plate->GetPos().z());
    gpu_sys.ApplyMeshMotion(plate_i, plate_pos, mesh_rot, mesh_lin_vel, mesh_ang_vel);
}

int main(int argc, char* argv[]) {
    gpu::SetDataPath(std::string(PROJECTS_DATA_DIR) + "gpu/");

    ChGpuSimulationParameters params;
    // Some of the default values might be overwritten by user via command line
    if (argc < 2 || argc > 2 && argc != num_args_full || ParseJSON(gpu::GetDataFile(argv[1]), params) == false) {
        ShowUsage(argv[0]);
        return 1;
    }

    int normal_stress_id;
    if (argc == num_args_full){
        normal_stress_id = std::atoi(argv[2]);
    }
    else{
        normal_stress_id = 0;
    }

    float iteration_step = params.step_size;  // TODO

    ChSystemGpuMesh gran_sys(params.sphere_radius, params.sphere_density,
                            make_float3(params.box_X, params.box_Y, params.box_Z));

    filesystem::create_directory(filesystem::path(params.output_dir));


    SetupGranSystem(gran_sys, params);
    gran_sys.Initialize();

    unsigned int numMeshes = gran_sys.GetNumMeshes();
    std::cout << numMeshes << " soup families" << std::endl;


    unsigned int currframe = 0;
    float out_fps = 100;
    float frame_step = 1.f / out_fps;  // Duration of a frame
    unsigned int out_steps = (unsigned int)(frame_step / iteration_step);
    std::cout << "out_steps " << out_steps << std::endl;

    double m_time = 0;
    unsigned int step = 0;

    ChSystemSMC ch_sys;
    const double gx = params.grav_X;
    const double gy = params.grav_Y;
    const double gz = params.grav_Z;
    double grav_mag = std::sqrt(gx * gx + gy * gy + gz * gz);

    ch_sys.Set_G_acc(ChVector<>(gx, gy, gz));

    auto plate = std::make_shared<ChBody>();
    plate->SetBodyFixed(true);
    plate->SetPos(ChVector<>(0, 0, params.box_Z));  // Initially out of the way
    plate_mass = normal_stresses[normal_stress_id] * box_xy * box_xy / grav_mag;
    plate->SetMass(plate_mass);
    ch_sys.AddBody(plate);

    SetInitialMeshes(gran_sys, plate);

    std::cout << "Running settling..." << std::endl;
    for (; m_time < time_settle; m_time += iteration_step, step++) {
        if (step % out_steps == 0) {
            std::cout << "Rendering frame " << currframe << std::endl;
            char filename[100];
            sprintf(filename, "%s/step%06u", params.output_dir.c_str(), currframe++);
            gran_sys.WriteParticleFile(std::string(filename));
            gran_sys.WriteMeshes(std::string(filename));
        }
        gran_sys.AdvanceSimulation(iteration_step);
    }

    // Add a weighted top plate
    double plate_z = gran_sys.GetMaxParticleZ() + 2 * params.sphere_radius;
    std::cout << "Adding plate at "
              << "(0, 0, " << plate_z << ")" << std::endl;
    plate->SetPos(ChVector<>(0, 0, plate_z));
    plate->SetBodyFixed(false);

    float* forces = new float[numMeshes * FAM_ENTRIES_FORCE];

    // Compress the material under the weight of the plate
    std::cout << "Running compression..." << std::endl;
    m_time = 0;
    ChVector<float> plate_pos(0, 0, 0);
    ChQuaternion<float> plate_quat(1, 0, 0, 0);
    ChVector<float> plate_lin_velo(0, 0, 0);
    ChVector<float> plate_ang_velo(0, 0, 0);
    ChVector<> plate_force;
    ChVector<> plate_torque;

    ChVector<float> top_pos(0, 0, 0);
    ChQuaternion<float> top_quat(1, 0, 0, 0);
    ChVector<float> top_lin_velo(0, 0, 0);
    ChVector<float> top_rot_velo(0, 0, 0);

    ChVector<float> bottom_pos(0, 0, 0);
    ChQuaternion<float> bottom_quat(1, 0, 0, 0);
    ChVector<float> bottom_lin_velo(0, 0, 0);
    ChVector<float> bottom_rot_velo(0, 0, 0);


    for (; m_time < time_compress; m_time += iteration_step, step++) {
        // Update Plate
        plate_pos.z() = plate->GetPos().z();
        plate_lin_velo.z() = plate->GetPos_dt().z();

        gran_sys.ApplyMeshMotion(plate_i, plate_pos, plate_quat, plate_lin_velo, plate_ang_velo);

        if (step % out_steps == 0) {
            std::cout << "Rendering frame " << currframe << std::endl;
            char filename[100];
            sprintf(filename, "%s/step%06u", params.output_dir.c_str(), currframe++);
            gran_sys.WriteParticleFile(std::string(filename));
            gran_sys.WriteMeshes(std::string(filename));
        }

        ch_sys.DoStepDynamics(iteration_step);
        gran_sys.AdvanceSimulation(iteration_step);

        gran_sys.CollectMeshContactForces(plate_i, plate_force, plate_torque);
        plate->Empty_forces_accumulators();
        // set force in x and y direction to zero
        plate_force.x() = 0;
        plate_force.y() = 0;
        plate->Accumulate_force(plate_force, plate->GetPos(), false);
    }

    std::cout << std::endl << "Running shear test..." << std::endl;
    // 5 Hz low pass filter
    // utils::ChButterworth_Lowpass fm_lowpass5(1, dt, 5.0);
    m_time = 0;
    for (; m_time < time_shear; step++, m_time += iteration_step) {
        double pos = m_time * shear_velocity;

        // Update Plate
        plate_pos.x() = pos;
        plate_pos.z() = plate->GetPos().z();
        plate_lin_velo.x() = shear_velocity;
        plate_lin_velo.z() = plate->GetPos_dt().z();
        gran_sys.ApplyMeshMotion(plate_i, plate_pos, plate_quat, plate_lin_velo, plate_ang_velo);
        
        // Update top bin
        top_pos.x() = pos;
        top_lin_velo.x() = shear_velocity;
        gran_sys.ApplyMeshMotion(top_i, top_pos, top_quat, top_lin_velo, top_rot_velo);        

        gran_sys.AdvanceSimulation(iteration_step);
        ch_sys.DoStepDynamics(iteration_step);

        gran_sys.CollectMeshContactForces(plate_i, plate_force, plate_torque);
        ChVector<> top_bin_force;
        ChVector<> top_bin_torque;
        gran_sys.CollectMeshContactForces(top_i, top_bin_force, top_bin_torque);
        double shear_force = plate_force.x() + top_bin_force.x();

        plate->Empty_forces_accumulators();
        plate->Accumulate_force(ChVector<>(0, 0, plate_force.z()), plate->GetPos(), false);

        // shear_force = fm_lowpass5.Filter(shear_force);

        // Output displacement and force
        if (step % out_steps == 0) {
            std::cout << "Rendering frame " << currframe << std::endl;
            char filename[100];
            sprintf(filename, "%s/step%06u", params.output_dir.c_str(), currframe++);
            gran_sys.WriteParticleFile(std::string(filename));
            gran_sys.WriteMeshes(std::string(filename));

            double shear_area = box_xy * (box_xy - m_time * shear_velocity * 2);
            double normal_stress = (plate_mass * grav_mag) / shear_area;
            double shear_stress = shear_force / shear_area;
            std::cout << std::setprecision(4) << "Time: " << m_time << std::endl;
            std::cout << std::setprecision(4) << "\tShear displacement: " << pos << std::endl;
            std::cout << std::setprecision(4) << "\tNormal stress: " << normal_stress << std::endl;
            std::cout << std::setprecision(4) << "\tShear stress: " << shear_stress << std::endl;
            std::cout << std::setprecision(4) << "\tShear stress / Normal stress: " << shear_stress / normal_stress
                      << std::endl;
        }
    }

    return 0;
}

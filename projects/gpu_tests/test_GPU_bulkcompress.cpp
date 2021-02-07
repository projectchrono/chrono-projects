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
// Chrono::Granular simulation of granular material being compressed by a mass
// allowed to oscillate on top of the material.
// =============================================================================

#include <cmath>
#include <iostream>
#include <string>

#include "chrono/physics/ChBody.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/utils/ChUtilsSamplers.h"

#include "chrono_gpu/ChGpuData.h"
#include "chrono_gpu/physics/ChSystemGpu.h"
#include "chrono_gpu/utils/ChGpuJsonParser.h"
#include "chrono_gpu/utils/ChGpuVisualization.h"

#include "chrono_thirdparty/filesystem/path.h"

using namespace chrono;
using namespace chrono::gpu;

/*
 * Lommen 2014:
 * Settle particles then drop box
 * 18,000 particles
 * Settled for 1 s
 * Poisson Ratio = 0.25
 * Coefficient of Restitution = 0.5
 * Static friction = 0.5
 * Rolling Friction = 0.04
 * Drop height (above material)= 0.2 m
 * Box mass = 50 kg
 * Periodic boundary condition in x and y
 * Time after release = 3 s
 * Varied Shear modulus
 *
 * Assumed:
 * Particle diameter = 4, 8, or 16 mm
 * Density 1500 kg/m3
 * Settled height ~0.25 m
 */

const double time_settle = 1;
double fill_top;

const double block_mass = 50000;  // 50kg

void ShowUsage(std::string name) {
    std::cout << "usage: " + name + " <json_file> <output_dir>" << std::endl;
}

void SetupGranSystem(ChSystemGpuMesh& gpu_sys, ChGpuSimulationParameters& params, std::string out_dir) {

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
    gpu_sys.SetFrictionMode(chrono::gpu::CHGPU_FRICTION_MODE::SINGLE_STEP);

    gpu_sys.SetGravitationalAcceleration(ChVector<float>(params.grav_X, params.grav_Y, params.grav_Z));

    gpu_sys.SetOutputMode(CHGPU_OUTPUT_MODE::CSV);
    filesystem::create_directory(filesystem::path(out_dir));

    gpu_sys.SetTimeIntegrator(CHGPU_TIME_INTEGRATOR::FORWARD_EULER);
    gpu_sys.SetFixedStepSize(params.step_size);
    gpu_sys.SetBDFixed(true);

    double epsilon = 0.2 * params.sphere_radius;
    double spacing = 2 * params.sphere_radius + epsilon;

    std::vector<ChVector<float>> body_points;

    // utils::HCPSampler<float> sampler(spacing);
    utils::PDSampler<float> sampler(2 * params.sphere_radius + epsilon);
    double fill_bottom = -params.box_Z / 2 + params.sphere_radius + epsilon;
    fill_top = params.box_Z / 2 - params.sphere_radius - epsilon;
    ChVector<> hdims(params.box_X / 2 - params.sphere_radius - epsilon,
                     params.box_Y / 2 - params.sphere_radius - epsilon, 0);
    for (double z = fill_bottom; z < fill_top; z += spacing) {
        ChVector<> center(0, 0, z);
        auto points = sampler.SampleBox(center, hdims);
        body_points.insert(body_points.end(), points.begin(), points.end());
    }

    std::cout << "Created " << body_points.size() << " spheres" << std::endl;

    gpu_sys.SetParticlePositions(body_points);

    // Mesh values
    std::vector<string> mesh_filenames;
    std::string mesh_filename = gpu::GetDataFile("test_GPU_bulkcompress/downward_square.obj");
    mesh_filenames.push_back(mesh_filename);

    std::vector<ChMatrix33<float>> mesh_rotscales;
    std::vector<float3> mesh_translations;

    ChMatrix33<float> scaling(ChVector<float>(params.box_X / 2, params.box_Y / 2, 1));
    mesh_rotscales.push_back(scaling);
    mesh_translations.push_back(make_float3(0, 0, 0));

    std::vector<float> mesh_masses;
    mesh_masses.push_back(block_mass);

    gpu_sys.LoadMeshes(mesh_filenames, mesh_rotscales, mesh_translations, mesh_masses);
}

int main(int argc, char* argv[]) {
    gpu::SetDataPath(std::string(PROJECTS_DATA_DIR) + "gpu/");
    ChGpuSimulationParameters params;

    if (argc != 3 || ParseJSON(argv[1], params) == false) {
        ShowUsage(argv[0]);
        return 1;
    }

    float iteration_step = params.step_size;
    ChSystemGpuMesh gpu_sys(params.sphere_radius, params.sphere_density,
                            make_float3(params.box_X, params.box_Y, params.box_Z));

    std::string out_dir(argv[2]);
    SetupGranSystem(gpu_sys, params, out_dir);

    ChSystemSMC ch_sys;
    ch_sys.Set_G_acc(ChVector<>(params.grav_X, params.grav_Y, params.grav_Z));
    auto block = std::make_shared<ChBody>();
    block->SetBodyFixed(true);
    block->SetPos(ChVector<>(0, 0, params.box_Z));
    block->SetMass(block_mass);
    ch_sys.AddBody(block);

    int numMeshes = gpu_sys.GetNumMeshes();
    std::cout << numMeshes << " meshes" << std::endl;

    unsigned int currframe = 0;
    double out_fps = 100;
    float frame_step = 1.f / out_fps;  // Duration of a frame
    unsigned int out_steps = frame_step / iteration_step;
    std::cout << "out_steps " << out_steps << std::endl;

    unsigned int step = 0;
    bool box_released = false;
    gpu_sys.EnableMeshCollision(true);
    gpu_sys.Initialize();

    for (float t = 0; t < params.time_end; t += iteration_step, step++) {
        if (t >= time_settle && box_released == false) {
            block->SetBodyFixed(false);
            double max_z = gpu_sys.GetMaxParticleZ();

            double max_velo = 100; // maximum velocity of the slab is 1m/s
            double drop_height = max_velo * max_velo / (2 * std::abs(params.grav_Z));    // 0.2m


            block->SetPos(ChVector<>(0, 0, max_z + params.sphere_radius + drop_height));

            box_released = true;
            std::cout << "Releasing box" << std::endl;
        }

        gpu_sys.ApplyMeshMotion(0, block->GetPos(), block->GetRot(), block->GetPos_dt(), block->GetWvel_par());
        
        ChVector<> box_force;
        ChVector<> box_torque;
        gpu_sys.CollectMeshContactForces(0, box_force, box_torque);

        block->Empty_forces_accumulators();
        block->Accumulate_force(box_force, block->GetPos(), false);
        block->Accumulate_torque(box_torque, false);

        if (step % out_steps == 0) {
            std::cout << "Rendering frame " << currframe << std::endl;
            char filename[100];
            sprintf(filename, "%s/step%06u", out_dir.c_str(), currframe++);
            gpu_sys.WriteFile(std::string(filename));
            gpu_sys.WriteMeshes(std::string(filename));
            if (box_released) {
                std::cout << block->GetPos().z() << std::endl;
            }
        }


        ch_sys.DoStepDynamics(iteration_step);
        gpu_sys.AdvanceSimulation(iteration_step);
    }

    return 0;
}
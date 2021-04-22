// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2020 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Nic Olsen
// =============================================================================
// Rolling ball down an incline
// =============================================================================

#include <cmath>
#include <iostream>
#include <string>

#include "GpuDemoUtils.h"
#include "chrono/utils/ChUtilsSamplers.h"

#include "chrono_gpu/physics/ChSystemGpu.h"
#include "chrono_gpu/utils/ChGpuJsonParser.h"
#include "chrono_gpu/ChGpuData.h"
#include "chrono_thirdparty/filesystem/path.h"

using namespace chrono;
using namespace chrono::gpu;

ChVector<float> sphere_pos(0, 0, 0);

enum RUN_MODE { NONE = 0, SCHWARTZ = 1 };

// expected number of args for param sweep
constexpr int num_args_full = 8;

void ShowUsage(std::string name) {
    std::cout << "usage: " + name +
                     " <json_file> <output_dir> <psi_L> <roll_mode: 0-none, 1-schwartz> <mu_roll> <angle> <v_init>"
              << std::endl;
    std::cout << "must have either 1 or " << num_args_full - 1 << " arguments" << std::endl;

}

int main(int argc, char* argv[]) {
    gpu::SetDataPath(std::string(PROJECTS_DATA_DIR) + "gpu/");

    // Some of the default values are overwritten by user via command line
    ChGpuSimulationParameters params;
    // Some of the default values might be overwritten by user via command line
    if (argc < 2 || argc > 2 && argc != num_args_full || ParseJSON(gpu::GetDataFile(argv[1]), params) == false) {
        ShowUsage(argv[0]);
        return 1;
    }

    RUN_MODE run_mode = RUN_MODE::SCHWARTZ;
    float mu_roll = 0.01;
    double theta = 5;
    float v_init_mag = 2;
    if (argc == num_args_full){
        params.output_dir = argv[2];
        params.psi_L = std::stoi(argv[3]);
        run_mode = (RUN_MODE)std::atoi(argv[4]);
        mu_roll = std::stof(argv[5]);
        theta = std::stod(argv[6]);
        v_init_mag = std::stof(argv[7]);

    }
    params.rolling_friction_coeffS2S = mu_roll;
    params.rolling_friction_coeffS2W = mu_roll;



    params.box_X = 60;
    params.box_Y = 60;
    params.box_Z = 60;

    // Setup simulation
    ChSystemGpu gran_sys(params.sphere_radius, params.sphere_density,
                        make_float3(params.box_X, params.box_Y, params.box_Z));
    gran_sys.DisableMinLength();

    switch (run_mode) {
        case RUN_MODE::NONE:
            gran_sys.SetRollingMode(CHGPU_ROLLING_MODE::NO_RESISTANCE);
            break;
        case RUN_MODE::SCHWARTZ:
            gran_sys.SetRollingMode(CHGPU_ROLLING_MODE::SCHWARTZ);
            gran_sys.SetRollingCoeff_SPH2SPH(params.rolling_friction_coeffS2S);
            gran_sys.SetRollingCoeff_SPH2WALL(params.rolling_friction_coeffS2W);
            break;
        default:
            std::cout << "Invalid run mode" << std::endl;
            return 1;
    }

    gran_sys.SetFrictionMode(CHGPU_FRICTION_MODE::MULTI_STEP);
    gran_sys.SetPsiFactors(params.psi_T, params.psi_L);

    // normal force model
    gran_sys.SetKn_SPH2SPH(params.normalStiffS2S);
    gran_sys.SetKn_SPH2WALL(params.normalStiffS2W);
    gran_sys.SetGn_SPH2SPH(params.normalDampS2S);
    gran_sys.SetGn_SPH2WALL(params.normalDampS2W);

    // tangential force model 
    gran_sys.SetKt_SPH2SPH(params.tangentStiffS2S);
    gran_sys.SetKt_SPH2WALL(params.tangentStiffS2W);
    gran_sys.SetGt_SPH2SPH(params.tangentDampS2S);
    gran_sys.SetGt_SPH2WALL(params.tangentDampS2W);

    gran_sys.SetStaticFrictionCoeff_SPH2SPH(params.static_friction_coeffS2S);
    gran_sys.SetStaticFrictionCoeff_SPH2WALL(params.static_friction_coeffS2W);

    gran_sys.SetCohesionRatio(0.0);
    gran_sys.SetAdhesionRatio_SPH2WALL(0.0);

    double x = -std::sqrt((2.0 * std::cos(theta) * std::cos(theta)) / (4.0 - 4.0 * std::cos(theta) * std::cos(theta)));
    // Vector parallel to the plane going up (ax, ax, az)
    double ax = x / std::sqrt(2.0 * x * x + 1.0);
    double az = 1.0 / std::sqrt(2.0 * x * x + 1.0);

    // Plane normal
    ChVector<float> n((float)(-az / (2.0 * ax)), (float)(-az / (2.0 * ax)), 1.f);
    n.Normalize();

    ChVector<float> v_init(ax, ax, az);
    v_init = (-v_init_mag / v_init.Length()) * v_init;

    ChVector<float> plane_pos(-params.sphere_radius * n.x(), -params.sphere_radius * n.y(), -params.sphere_radius * n.z());
    ChVector<float> plane_normal(n.x(), n.y(), n.z());
    bool track_forces = false;
    gran_sys.CreateBCPlane(plane_pos, plane_normal, track_forces);

    gran_sys.SetGravitationalAcceleration(ChVector<float>(params.grav_X, params.grav_Y, params.grav_Z));
    gran_sys.SetOutputMode(params.write_mode);
    gran_sys.SetOutputFlags(CHGPU_OUTPUT_FLAGS::VEL_COMPONENTS | CHGPU_OUTPUT_FLAGS::ANG_VEL_COMPONENTS);

    gran_sys.SetTimeIntegrator(CHGPU_TIME_INTEGRATOR::CENTERED_DIFFERENCE);

    std::vector<ChVector<float>> body_points;
    body_points.push_back(sphere_pos);
    std::vector<ChVector<float>> body_vels(1, v_init);

    gran_sys.SetParticlePositions(body_points, body_vels);

    gran_sys.SetFixedStepSize(params.step_size);

    filesystem::create_directory(filesystem::path(params.output_dir));

    gran_sys.SetBDFixed(true);

    gran_sys.SetVerbosity(params.verbose);
    gran_sys.Initialize();

    int fps = 1000;
    float frame_step = 1.f / fps;
    float curr_time = 0;
    int currframe = 0;

    // write an initial frame
    char filename[100];
    sprintf(filename, "%s/step%06d", params.output_dir.c_str(), currframe++);
    gran_sys.WriteParticleFile(std::string(filename));

    std::cout << "frame step is " << frame_step << std::endl;
    while (curr_time < params.time_end) {
        gran_sys.AdvanceSimulation(frame_step);
        curr_time += frame_step;
        printf("rendering frame %u\n", currframe);
        sprintf(filename, "%s/step%06d", params.output_dir.c_str(), currframe++);
        gran_sys.WriteParticleFile(std::string(filename));
    }

    return 0;
}

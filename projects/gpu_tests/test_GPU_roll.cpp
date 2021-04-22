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
// Authors: Nic Olsen, Ruochun Zhang
// =============================================================================
// Rolling ball on a horizontal plane varying rolling friction. Base case is
// a ball rolling with no slip and no rolling friction.
// =============================================================================

#include <cmath>
#include <iostream>
#include <string>

#include "GpuDemoUtils.h"
#include "chrono/utils/ChUtilsSamplers.h"
#include "chrono_gpu/ChGpuData.h"
#include "chrono_gpu/physics/ChSystemGpu.h"
#include "chrono_gpu/utils/ChGpuJsonParser.h"
#include "chrono_thirdparty/filesystem/path.h"

using namespace chrono;
using namespace chrono::gpu;

ChVector<float> sphere_pos(0, 0, 0);
ChVector<float> v_init(10, -10, 0);

enum RUN_MODE { NONE = 0, SCHWARTZ = 1 };

void ShowUsage(std::string name) {
    std::cout << "usage: " + name + " <json_file> <output_dir> <psi_L> <roll_mode: 0-none, 1-schwartz> <mu_roll>"
              << std::endl;
}

int main(int argc, char* argv[]) {
    gpu::SetDataPath(std::string(PROJECTS_DATA_DIR) + "gpu/");

    // Some of the default values are overwritten by user via command line
    ChGpuSimulationParameters params;
    if (argc != 6 || ParseJSON(gpu::GetDataFile(argv[1]), params) == false) {
        ShowUsage(argv[0]);
        return 1;
    }
    params.output_dir = argv[2];
    params.psi_L = std::stoi(argv[3]);
    RUN_MODE run_mode = (RUN_MODE)std::atoi(argv[4]);
    float mu_roll = std::stof(argv[5]);
    params.rolling_friction_coeffS2S = mu_roll;
    params.rolling_friction_coeffS2W = mu_roll;

    params.box_X = 15;
    params.box_Y = 15;
    params.box_Z = 15;

    // Angle gravity
    params.grav_X = -565.80;
    params.grav_Y = -565.80;
    params.grav_Z = -565.80;

    // Setup simulation
    ChSystemGpu gpu_sys(params.sphere_radius, params.sphere_density,
                                 make_float3(params.box_X, params.box_Y, params.box_Z));
    gpu_sys.DisableMinLength();
    switch (run_mode) {
        case RUN_MODE::NONE:
            gpu_sys.SetRollingMode(CHGPU_ROLLING_MODE::NO_RESISTANCE);
            break;
        case RUN_MODE::SCHWARTZ:
            gpu_sys.SetRollingMode(CHGPU_ROLLING_MODE::SCHWARTZ);
            gpu_sys.SetRollingCoeff_SPH2WALL(params.rolling_friction_coeffS2W);
            gpu_sys.SetRollingCoeff_SPH2SPH(params.rolling_friction_coeffS2S);
            break;
        default:
            std::cout << "Invalid run mode" << std::endl;
            return 1;
    }

    gpu_sys.SetFrictionMode(CHGPU_FRICTION_MODE::MULTI_STEP);

    gpu_sys.SetPsiFactors(params.psi_T, params.psi_L);

    gpu_sys.SetKn_SPH2SPH(params.normalStiffS2S);
    gpu_sys.SetKn_SPH2WALL(params.normalStiffS2W);
    gpu_sys.SetGn_SPH2SPH(params.normalDampS2S);
    gpu_sys.SetGn_SPH2WALL(params.normalDampS2W);

    gpu_sys.SetKt_SPH2SPH(params.tangentStiffS2S);
    gpu_sys.SetKt_SPH2WALL(params.tangentStiffS2W);

    gpu_sys.SetGt_SPH2SPH(params.tangentDampS2S);
    gpu_sys.SetGt_SPH2WALL(params.tangentDampS2W);

    gpu_sys.SetStaticFrictionCoeff_SPH2SPH(params.static_friction_coeffS2S);
    gpu_sys.SetStaticFrictionCoeff_SPH2WALL(params.static_friction_coeffS2W);

    gpu_sys.SetCohesionRatio(params.cohesion_ratio);
    gpu_sys.SetAdhesionRatio_SPH2WALL(params.adhesion_ratio_s2w);

    // Plane normal
    ChVector<float> n(1, 1, 1);
    n.Normalize();
    ChVector<> plane_pos(sphere_pos.x() - params.sphere_radius * n.x(), sphere_pos.y() - params.sphere_radius * n.y(),
                         sphere_pos.z() - params.sphere_radius * n.z());
    ChVector<> plane_normal(n.x(), n.y(), n.z());
    bool track_forces = false;
    gpu_sys.CreateBCPlane(plane_pos, plane_normal, track_forces);

    gpu_sys.SetGravitationalAcceleration(ChVector<>(params.grav_X, params.grav_Y, params.grav_Z));
    gpu_sys.SetOutputMode(params.write_mode);
    gpu_sys.SetOutputFlags(CHGPU_OUTPUT_FLAGS::VEL_COMPONENTS | CHGPU_OUTPUT_FLAGS::ANG_VEL_COMPONENTS);

    gpu_sys.SetTimeIntegrator(CHGPU_TIME_INTEGRATOR::CENTERED_DIFFERENCE);

    std::vector<ChVector<float>> body_points;
    body_points.push_back(sphere_pos);

    gpu_sys.SetParticlePositions(body_points, std::vector<ChVector<float>>(1, v_init));

    gpu_sys.SetFixedStepSize(params.step_size);

    filesystem::create_directory(filesystem::path(params.output_dir));

    gpu_sys.SetBDFixed(true);

    gpu_sys.SetVerbosity(params.verbose);
    gpu_sys.Initialize();

    int fps = 1000;
    float frame_step = 1.f / fps;
    float curr_time = 0;
    int currframe = 0;

    // write an initial frame
    char filename[100];
    sprintf(filename, "%s/step%06d", params.output_dir.c_str(), currframe++);
    gpu_sys.WriteParticleFile(std::string(filename));

    std::cout << "frame step is " << frame_step << std::endl;
    while (curr_time < params.time_end) {
        float real_dt = gpu_sys.AdvanceSimulation(frame_step);

        curr_time += frame_step;
        printf("rendering frame %u\n", currframe);
        sprintf(filename, "%s/step%06d", params.output_dir.c_str(), currframe++);
        gpu_sys.WriteParticleFile(std::string(filename));
    }

    return 0;
}

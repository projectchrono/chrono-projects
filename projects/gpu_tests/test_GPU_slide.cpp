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
// Authors: Nic Olsen, Ruochun
// =============================================================================
// Sliding / rolling ball on a horizontal plane
// =============================================================================

#include <cmath>
#include <iostream>
#include <string>

#include "GpuDemoUtils.h"
#include "chrono/utils/ChUtilsSamplers.h"
#include "chrono_gpu/physics/ChSystemGpu.h"
#include "chrono_gpu/utils/ChGpuJsonParser.h"
#include "chrono_thirdparty/filesystem/path.h"

#include "../utils.h"

using namespace chrono;
using namespace chrono::gpu;

enum RUN_MODE { FRICTIONLESS = 0, ZERO_FRICTION = 1, SMALL_FRICTION = 2, LARGE_FRICTION = 3 };
float mu_small = 1e-5f;
float mu_large = 0.5;

bool axis_aligned = false;
ChVector<float> sphere_pos(0, 0, 0);
ChVector<float> v_init(10, -10, 0);

int main(int argc, char* argv[]) {
    std::string inputJson = GetProjectsDataFile("gpu/Slide.json");
    RUN_MODE run_mode = RUN_MODE::FRICTIONLESS;
    if (argc == 2) {
        inputJson = std::string(argv[1]);
    } else if (argc == 3) {
        inputJson = std::string(argv[1]);
        run_mode = (RUN_MODE)std::atoi(argv[2]);
    } else if (argc > 1) {
        std::cout << "Usage:\n./test_GPU_slide <json_file> [<run_mode>]" << std::endl;
        std::cout << "  run_mode: 0-frictionless, 1-zero_friction, 2-small_friction, 3-large_friction" << std::endl;
        return 1;
    }

    ChGpuSimulationParameters params;
    if (!ParseJSON(inputJson, params)) {
        std ::cout << "ERROR: reading input file " << inputJson << std::endl;
        return 1;
    }

    params.box_X = 15;
    params.box_Y = 15;
    params.box_Z = 15;

    // Angle gravity
    params.grav_X = -565.80f;
    params.grav_Y = -565.80f;
    params.grav_Z = -565.80f;

    // Setup simulation
    ChSystemGpu gpu_sys(params.sphere_radius, params.sphere_density,
                        ChVector<float>(params.box_X, params.box_Y, params.box_Z));
    gpu_sys.DisableMinLength();
    switch (run_mode) {
        case RUN_MODE::FRICTIONLESS: {
            std::cout << "Frictionless" << std::endl;
            gpu_sys.SetFrictionMode(CHGPU_FRICTION_MODE::FRICTIONLESS);
            break;
        }
        case RUN_MODE::ZERO_FRICTION: {
            std::cout << "Zero Friction" << std::endl;
            gpu_sys.SetFrictionMode(CHGPU_FRICTION_MODE::MULTI_STEP);
            params.static_friction_coeffS2W = 0.f;
            break;
        }
        case RUN_MODE::SMALL_FRICTION: {
            std::cout << "Small Friction " << mu_small << std::endl;
            gpu_sys.SetFrictionMode(CHGPU_FRICTION_MODE::MULTI_STEP);
            params.static_friction_coeffS2W = mu_small;
            break;
        }
        case RUN_MODE::LARGE_FRICTION: {
            std::cout << "Large Friction " << mu_large << std::endl;
            gpu_sys.SetFrictionMode(CHGPU_FRICTION_MODE::MULTI_STEP);
            params.static_friction_coeffS2W = mu_large;
            break;
        }
    }

    if (axis_aligned) {
        std::cout << "Axis Aligned" << std::endl;
        params.grav_X = 0;
        params.grav_Y = 0;
        params.grav_Z = -980;

        sphere_pos.x() = 0;
        sphere_pos.y() = 0;
        sphere_pos.z() = -params.box_X / 2 + params.sphere_radius;

        v_init.x() = 1;
        v_init.y() = 0;
        v_init.z() = 0;
    }

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

    gpu_sys.SetRollingMode(CHGPU_ROLLING_MODE::NO_RESISTANCE);

    gpu_sys.SetCohesionRatio(params.cohesion_ratio);
    gpu_sys.SetAdhesionRatio_SPH2WALL(params.adhesion_ratio_s2w);

    // Plane normal
    ChVector<float> n(1, 1, 1);
    n.Normalize();

    ChVector<> plane_pos(sphere_pos.x() - params.sphere_radius * n.x(), sphere_pos.y() - params.sphere_radius * n.y(),
                         sphere_pos.z() - params.sphere_radius * n.z());
    ChVector<> plane_normal(n.x(), n.y(), n.z());
    bool track_forces = false;
    if (!axis_aligned) {
        gpu_sys.CreateBCPlane(plane_pos, plane_normal, track_forces);
    }
    gpu_sys.SetGravitationalAcceleration(ChVector<>(params.grav_X, params.grav_Y, params.grav_Z));
    gpu_sys.SetParticleOutputMode(params.write_mode);
    gpu_sys.SetParticleOutputFlags(CHGPU_OUTPUT_FLAGS::VEL_COMPONENTS | CHGPU_OUTPUT_FLAGS::ANG_VEL_COMPONENTS);

    gpu_sys.SetTimeIntegrator(CHGPU_TIME_INTEGRATOR::CENTERED_DIFFERENCE);

    std::vector<ChVector<float>> body_points;
    body_points.push_back(sphere_pos);

    gpu_sys.SetParticles(body_points, std::vector<ChVector<float>>(1, v_init));

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
        curr_time += frame_step;
        printf("rendering frame %u\n", currframe);
        sprintf(filename, "%s/step%06d", params.output_dir.c_str(), currframe++);
        gpu_sys.WriteParticleFile(std::string(filename));
    }

    return 0;
}

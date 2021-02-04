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
// Authors: Nic Olsen,  Ruochun Zhang
// =============================================================================
// Normal contact sphere-vs-wall and sphere-vs-fixed-sphere
// =============================================================================

#include <cmath>
#include <iostream>
#include <string>

#include "GpuDemoUtils.h"
#include "chrono/utils/ChUtilsSamplers.h"
#include "chrono_gpu/physics/ChSystemGpu.h"
#include "chrono_gpu/utils/ChGpuJsonParser.h"
#include "chrono_thirdparty/filesystem/path.h"

using namespace chrono;
using namespace chrono::gpu;

enum RUN_MODE { SPHERE_WALL = 0, SPHERE_SPHERE = 1, SPHERE_WALL_ANGLE = 2 };

void ShowUsage(std::string name) {
    std::cout << "usage: " + name + " <json_file> <output_dir> <psi_L> <run_mode> <gamma_n>" << std::endl;
}

int main(int argc, char* argv[]) {
    ChGpuSimulationParameters params;

    // Some of the default values are overwritten by user via command line
    if (argc != 6 || ParseJSON(argv[1], params) == false) {
        ShowUsage(argv[0]);
        return 1;
    }
    params.output_dir = argv[2];
    params.psi_L = std::stoi(argv[3]);
    RUN_MODE run_mode = (RUN_MODE)std::stoi(argv[4]);

    params.box_X = 10;
    params.box_Y = 10;
    params.box_Z = (run_mode == SPHERE_SPHERE) ? 8 * params.sphere_radius : 4 * params.sphere_radius;

    float gamma_n = std::stof(argv[5]);
    params.normalDampS2S = gamma_n;
    params.normalDampS2W = gamma_n;
    std::cout << "Gamma " << gamma_n << std::endl;

    // Setup simulation
    ChSystemGpu gpu_sys(params.sphere_radius, params.sphere_density,
                                 make_float3(params.box_X, params.box_Y, params.box_Z));
    gpu_sys.DisableMinLength();

    gpu_sys.SetPsiFactors(params.psi_T, params.psi_L);

    gpu_sys.SetKn_SPH2SPH(params.normalStiffS2S);
    gpu_sys.SetKn_SPH2WALL(params.normalStiffS2W);
    gpu_sys.SetGn_SPH2SPH(params.normalDampS2S);
    gpu_sys.SetGn_SPH2WALL(params.normalDampS2W);

    gpu_sys.SetCohesionRatio(params.cohesion_ratio);
    gpu_sys.SetAdhesionRatio_SPH2WALL(params.adhesion_ratio_s2w);
    if (run_mode == SPHERE_WALL_ANGLE) {
        params.grav_X = -565.80;
        params.grav_Y = -565.80;
        params.grav_Z = -565.80;

        ChVector<> plane_pos((float)(-2 * params.sphere_radius / std::sqrt(3)),
                             (float)(-2 * params.sphere_radius / std::sqrt(3)),
                             (float)(-2 * params.sphere_radius / std::sqrt(3)));

        ChVector<> plane_normal(1, 1, 1);
        bool track_forces = false;
        gpu_sys.CreateBCPlane(plane_pos, plane_normal, track_forces);
    }

    gpu_sys.SetGravitationalAcceleration(ChVector<>(params.grav_X, params.grav_Y, params.grav_Z));
    gpu_sys.SetOutputMode(params.write_mode);
    gpu_sys.SetOutputFlags(
        CHGPU_OUTPUT_FLAGS::VEL_COMPONENTS | CHGPU_OUTPUT_FLAGS::FIXITY |
        CHGPU_OUTPUT_FLAGS::FORCE_COMPONENTS);  // NOTE: original test used custom FORCE_COMPONENTS output

    gpu_sys.SetFrictionMode(CHGPU_FRICTION_MODE::FRICTIONLESS);
    gpu_sys.SetTimeIntegrator(CHGPU_TIME_INTEGRATOR::CENTERED_DIFFERENCE);

    std::vector<ChVector<float>> body_points;
    body_points.push_back(ChVector<float>(0, 0, 0));

    if (run_mode == SPHERE_SPHERE) {
        std::vector<bool> body_points_fixed;
        body_points.push_back(ChVector<float>(0, 0, -3 * params.sphere_radius));
        body_points_fixed.push_back(false);
        body_points_fixed.push_back(true);
        gpu_sys.SetParticleFixed(body_points_fixed);
    }

    gpu_sys.SetParticlePositions(body_points);

    gpu_sys.SetFixedStepSize(params.step_size);

    filesystem::create_directory(filesystem::path(params.output_dir));

    gpu_sys.SetBDFixed(true);

    gpu_sys.SetVerbosity(params.verbose);
    gpu_sys.Initialize();

    int fps = 10000;
    float frame_step = 1.f / fps;
    float curr_time = 0;
    int currframe = 0;

    // write an initial frame
    char filename[100];
    sprintf(filename, "%s/step%06d", params.output_dir.c_str(), currframe++);
    gpu_sys.WriteFile(std::string(filename));

    std::cout << "frame step is " << frame_step << std::endl;
    while (curr_time < params.time_end) {
        gpu_sys.AdvanceSimulation(frame_step);
        curr_time += frame_step;
        printf("rendering frame %u\n", currframe);
        sprintf(filename, "%s/step%06d", params.output_dir.c_str(), currframe++);
        gpu_sys.WriteFile(std::string(filename));
    }

    return 0;
}

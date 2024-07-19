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
#include "chrono_thirdparty/filesystem/path.h"

#include "../utils.h"

using namespace chrono;
using namespace chrono::gpu;

ChVector3f sphere_pos(0, 0, 0);

enum ROLL_MODE { NONE = 0, SCHWARTZ = 1 };

int main(int argc, char* argv[]) {
    std::string inputJson = GetProjectsDataFile("gpu/Incline.json");
    ROLL_MODE roll_mode = ROLL_MODE::SCHWARTZ;
    double theta = 5;
    float v_init_mag = 2;

    if (argc == 2) {
        inputJson = std::string(argv[1]);
    } else if (argc == 5) {
        inputJson = std::string(argv[1]);
        roll_mode = (ROLL_MODE)std::atoi(argv[2]);
        theta = std::stod(argv[3]);
        v_init_mag = std::stof(argv[4]);
    } else if (argc > 1) {
        std::cout << "Usage:\n./test_GPU_incline <json_file> [<roll_mode> <angle> <vinit>]" << std::endl;
        std::cout << "  roll_mode: 0 - none, 1 - schwartz" << std::endl;
        return 1;
    }

    ChGpuSimulationParameters params;
    if (!ParseJSON(inputJson, params)) {
        std ::cout << "ERROR: reading input file " << inputJson << std::endl;
        return 1;
    }

    params.box_X = 60;
    params.box_Y = 60;
    params.box_Z = 60;

    // Setup simulation
    ChSystemGpu gran_sys(params.sphere_radius, params.sphere_density,
                         ChVector3f(params.box_X, params.box_Y, params.box_Z));
    gran_sys.DisableMinLength();

    switch (roll_mode) {
        case ROLL_MODE::NONE:
            gran_sys.SetRollingMode(CHGPU_ROLLING_MODE::NO_RESISTANCE);
            break;
        case ROLL_MODE::SCHWARTZ:
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
    ChVector3f n((float)(-az / (2.0 * ax)), (float)(-az / (2.0 * ax)), 1.f);
    n.Normalize();

    ChVector3f v_init(ax, ax, az);
    v_init = (-v_init_mag / v_init.Length()) * v_init;

    ChVector3f plane_pos(-params.sphere_radius * n.x(), -params.sphere_radius * n.y(),
                              -params.sphere_radius * n.z());
    ChVector3f plane_normal(n.x(), n.y(), n.z());
    bool track_forces = false;
    gran_sys.CreateBCPlane(plane_pos, plane_normal, track_forces);

    gran_sys.SetGravitationalAcceleration(ChVector3f(params.grav_X, params.grav_Y, params.grav_Z));
    gran_sys.SetParticleOutputMode(params.write_mode);
    gran_sys.SetParticleOutputFlags(CHGPU_OUTPUT_FLAGS::VEL_COMPONENTS | CHGPU_OUTPUT_FLAGS::ANG_VEL_COMPONENTS);

    gran_sys.SetTimeIntegrator(CHGPU_TIME_INTEGRATOR::CENTERED_DIFFERENCE);

    std::vector<ChVector3f> body_points;
    body_points.push_back(sphere_pos);
    std::vector<ChVector3f> body_vels(1, v_init);

    gran_sys.SetParticles(body_points, body_vels);

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

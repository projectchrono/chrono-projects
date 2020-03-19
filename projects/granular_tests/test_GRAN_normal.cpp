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
// Normal contact sphere-vs-wall and sphere-vs-fixed-sphere
// =============================================================================

#include <cmath>
#include <iostream>
#include <string>

#include "ChGranularDemoUtils.hpp"
#include "chrono/utils/ChUtilsSamplers.h"
#include "chrono_granular/api/ChApiGranularChrono.h"
#include "chrono_granular/physics/ChGranular.h"
#include "chrono_granular/utils/ChGranularJsonParser.h"
#include "chrono_thirdparty/filesystem/path.h"

using namespace chrono;
using namespace chrono::granular;

enum RUN_MODE { SPHERE_WALL = 0, SPHERE_SPHERE = 1, SPHERE_WALL_ANGLE = 2 };

void ShowUsage(std::string name) {
    std::cout << "usage: " + name + " <json_file> <output_dir> <psi_L> <run_mode> <gamma_n>" << std::endl;
}

int main(int argc, char* argv[]) {
    sim_param_holder params;

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
    ChSystemGranularSMC gran_sys(params.sphere_radius, params.sphere_density,
                                 make_float3(params.box_X, params.box_Y, params.box_Z));
    gran_sys.disableMinLength();
    ChGranularSMC_API apiSMC;
    apiSMC.setGranSystem(&gran_sys);

    gran_sys.setPsiFactors(params.psi_T, params.psi_L);

    gran_sys.set_K_n_SPH2SPH(params.normalStiffS2S);
    gran_sys.set_K_n_SPH2WALL(params.normalStiffS2W);
    gran_sys.set_Gamma_n_SPH2SPH(params.normalDampS2S);
    gran_sys.set_Gamma_n_SPH2WALL(params.normalDampS2W);

    gran_sys.set_Cohesion_ratio(params.cohesion_ratio);
    gran_sys.set_Adhesion_ratio_S2W(params.adhesion_ratio_s2w);
    if (run_mode == SPHERE_WALL_ANGLE) {
        params.grav_X = -565.80;
        params.grav_Y = -565.80;
        params.grav_Z = -565.80;

        float plane_pos[] = {(float)(-2 * params.sphere_radius / std::sqrt(3)),
                             (float)(-2 * params.sphere_radius / std::sqrt(3)),
                             (float)(-2 * params.sphere_radius / std::sqrt(3))};

        float plane_normal[] = {1, 1, 1};
        bool track_forces = false;
        gran_sys.Create_BC_Plane(plane_pos, plane_normal, track_forces);
    }

    gran_sys.set_gravitational_acceleration(params.grav_X, params.grav_Y, params.grav_Z);
    gran_sys.setOutputMode(params.write_mode);
    gran_sys.setOutputFlags(
        GRAN_OUTPUT_FLAGS::VEL_COMPONENTS | GRAN_OUTPUT_FLAGS::FIXITY |
        GRAN_OUTPUT_FLAGS::FORCE_COMPONENTS);  // NOTE: original test used custom FORCE_COMPONENTS output

    gran_sys.set_friction_mode(GRAN_FRICTION_MODE::FRICTIONLESS);
    gran_sys.set_timeIntegrator(GRAN_TIME_INTEGRATOR::CENTERED_DIFFERENCE);

    std::vector<ChVector<float>> body_points;
    body_points.push_back(ChVector<float>(0, 0, 0));

    if (run_mode == SPHERE_SPHERE) {
        std::vector<bool> body_points_fixed;
        body_points.push_back(ChVector<float>(0, 0, -3 * params.sphere_radius));
        body_points_fixed.push_back(false);
        body_points_fixed.push_back(true);
        gran_sys.setParticleFixed(body_points_fixed);
    }

    apiSMC.setElemsPositions(body_points);

    gran_sys.set_fixed_stepSize(params.step_size);

    filesystem::create_directory(filesystem::path(params.output_dir));

    gran_sys.set_BD_Fixed(true);

    gran_sys.setVerbose(params.verbose);
    gran_sys.initialize();

    int fps = 10000;
    float frame_step = 1.f / fps;
    float curr_time = 0;
    int currframe = 0;

    // write an initial frame
    char filename[100];
    sprintf(filename, "%s/step%06d", params.output_dir.c_str(), currframe++);
    gran_sys.writeFile(std::string(filename));

    std::cout << "frame step is " << frame_step << std::endl;
    while (curr_time < params.time_end) {
        gran_sys.advance_simulation(frame_step);
        curr_time += frame_step;
        printf("rendering frame %u\n", currframe);
        sprintf(filename, "%s/step%06d", params.output_dir.c_str(), currframe++);
        gran_sys.writeFile(std::string(filename));
    }

    return 0;
}
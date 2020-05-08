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
// Sliding / rolling ball on a horizontal plane
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

enum RUN_MODE { FRICTIONLESS = 0, ZERO_FRICTION = 1, SMALL_FRICTION = 2, LARGE_FRICTION = 3 };
float mu_small = 1e-5f;
float mu_large = 0.5;

bool axis_aligned = false;
ChVector<float> sphere_pos(0, 0, 0);
ChVector<float> v_init(10, -10, 0);

void ShowUsage(std::string name) {
    std::cout << "usage: " + name +
                     " <json_file> <output_dir> <psi_L> <run_mode: 0-frictionless, 1-zero_friction, 2-small_friction, "
                     "3-large_friction>"
              << std::endl;
}

int main(int argc, char* argv[]) {
    sim_param_holder params;

    // Some of the default values are overwritten by user via command line
    if (argc != 5 || ParseJSON(argv[1], params) == false) {
        ShowUsage(argv[0]);
        return 1;
    }
    params.output_dir = argv[2];
    params.psi_L = std::stoi(argv[3]);

    RUN_MODE run_mode = (RUN_MODE)std::stoi(argv[4]);

    params.box_X = 15;
    params.box_Y = 15;
    params.box_Z = 15;

    // Angle gravity
    params.grav_X = -565.80;
    params.grav_Y = -565.80;
    params.grav_Z = -565.80;

    // Setup simulation
    ChSystemGranularSMC gran_sys(params.sphere_radius, params.sphere_density,
                                 make_float3(params.box_X, params.box_Y, params.box_Z));
    gran_sys.disableMinLength();
    switch (run_mode) {
        case RUN_MODE::FRICTIONLESS: {
            std::cout << "Frictionless" << std::endl;
            gran_sys.set_friction_mode(GRAN_FRICTION_MODE::FRICTIONLESS);
            break;
        }
        case RUN_MODE::ZERO_FRICTION: {
            std::cout << "Zero Friction" << std::endl;
            gran_sys.set_friction_mode(GRAN_FRICTION_MODE::MULTI_STEP);
            params.static_friction_coeffS2W = 0.f;
            break;
        }
        case RUN_MODE::SMALL_FRICTION: {
            std::cout << "Small Friction " << mu_small << std::endl;
            gran_sys.set_friction_mode(GRAN_FRICTION_MODE::MULTI_STEP);
            params.static_friction_coeffS2W = mu_small;
            break;
        }
        case RUN_MODE::LARGE_FRICTION: {
            std::cout << "Large Friction " << mu_large << std::endl;
            gran_sys.set_friction_mode(GRAN_FRICTION_MODE::MULTI_STEP);
            params.static_friction_coeffS2W = mu_large;
            break;
        }
        default: {
            ShowUsage(argv[0]);
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

    ChGranularSMC_API apiSMC;
    apiSMC.setGranSystem(&gran_sys);

    gran_sys.setPsiFactors(params.psi_T, params.psi_L);

    gran_sys.set_K_n_SPH2SPH(params.normalStiffS2S);
    gran_sys.set_K_n_SPH2WALL(params.normalStiffS2W);
    gran_sys.set_Gamma_n_SPH2SPH(params.normalDampS2S);
    gran_sys.set_Gamma_n_SPH2WALL(params.normalDampS2W);

    gran_sys.set_K_t_SPH2SPH(params.tangentStiffS2S);
    gran_sys.set_K_t_SPH2WALL(params.tangentStiffS2W);

    gran_sys.set_Gamma_t_SPH2SPH(params.tangentDampS2S);
    gran_sys.set_Gamma_t_SPH2WALL(params.tangentDampS2W);

    gran_sys.set_static_friction_coeff_SPH2SPH(params.static_friction_coeffS2S);
    gran_sys.set_static_friction_coeff_SPH2WALL(params.static_friction_coeffS2W);

    gran_sys.set_rolling_mode(GRAN_ROLLING_MODE::NO_RESISTANCE);

    gran_sys.set_Cohesion_ratio(params.cohesion_ratio);
    gran_sys.set_Adhesion_ratio_S2W(params.adhesion_ratio_s2w);

    // Plane normal
    ChVector<float> n(1, 1, 1);
    n.Normalize();

    float plane_pos[] = {sphere_pos.x() - params.sphere_radius * n.x(), sphere_pos.y() - params.sphere_radius * n.y(),
                         sphere_pos.z() - params.sphere_radius * n.z()};
    float plane_normal[] = {n.x(), n.y(), n.z()};
    bool track_forces = false;
    if (!axis_aligned) {
        gran_sys.Create_BC_Plane(plane_pos, plane_normal, track_forces);
    }
    gran_sys.set_gravitational_acceleration(params.grav_X, params.grav_Y, params.grav_Z);
    gran_sys.setOutputMode(params.write_mode);
    gran_sys.setOutputFlags(GRAN_OUTPUT_FLAGS::VEL_COMPONENTS | GRAN_OUTPUT_FLAGS::ANG_VEL_COMPONENTS);

    gran_sys.set_timeIntegrator(GRAN_TIME_INTEGRATOR::CENTERED_DIFFERENCE);

    std::vector<ChVector<float>> body_points;
    body_points.push_back(sphere_pos);

    apiSMC.setElemsPositions(body_points, std::vector<ChVector<float>>(1, v_init));

    gran_sys.set_fixed_stepSize(params.step_size);

    filesystem::create_directory(filesystem::path(params.output_dir));

    gran_sys.set_BD_Fixed(true);

    gran_sys.setVerbose(params.verbose);
    gran_sys.initialize();

    int fps = 1000;
    float frame_step = 1.f / fps;
    float curr_time = 0;
    int currframe = 0;

    // write an initial frame
    char filename[100];
    sprintf(filename, "%s/step%06d", params.output_dir.c_str(), currframe++);
    gran_sys.writeFile(std::string(filename));

    std::cout << "frame step is " << frame_step << std::endl;
    while (curr_time < params.time_end) {
        float real_dt = gran_sys.advance_simulation(frame_step);

        curr_time += frame_step;
        printf("rendering frame %u\n", currframe);
        sprintf(filename, "%s/step%06d", params.output_dir.c_str(), currframe++);
        gran_sys.writeFile(std::string(filename));
    }

    return 0;
}
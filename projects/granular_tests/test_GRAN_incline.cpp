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

#include <iostream>
#include <string>
#include <cmath>
#include "chrono_thirdparty/filesystem/path.h"
#include "chrono_granular/api/ChApiGranularChrono.h"
#include "chrono_granular/physics/ChGranular.h"
#include "chrono_granular/utils/ChGranularJsonParser.h"
#include "ChGranularDemoUtils.hpp"
#include "chrono/utils/ChUtilsSamplers.h"

using namespace chrono;
using namespace chrono::granular;

void ShowUsage(std::string name) {
    std::cout << "usage: " + name + " <json_file> <output_dir> <psi_L> <incline angle>" << std::endl;
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

    double theta = std::stod(argv[4]);

    // TODO
    params.box_X = 10;
    params.box_Y = 10;
    params.box_Z = 10;

    // Setup simulation
    ChSystemGranularSMC gran_sys(params.sphere_radius, params.sphere_density,
                                 make_float3(params.box_X, params.box_Y, params.box_Z));

    ChGranularSMC_API apiSMC;
    apiSMC.setGranSystem(&gran_sys);

    gran_sys.setPsiFactors(params.psi_T, params.psi_L);

    gran_sys.set_K_n_SPH2SPH(params.normalStiffS2S);
    gran_sys.set_K_n_SPH2WALL(params.normalStiffS2W);
    gran_sys.set_Gamma_n_SPH2SPH(params.normalDampS2S);
    gran_sys.set_Gamma_n_SPH2WALL(params.normalDampS2W);

    gran_sys.set_K_t_SPH2SPH(params.tangentStiffS2S); // TODO play with this???
    gran_sys.set_K_t_SPH2WALL(params.tangentStiffS2W);
    
    gran_sys.set_Gamma_t_SPH2SPH(params.tangentDampS2S);
    gran_sys.set_Gamma_t_SPH2WALL(params.tangentDampS2W);
    
    gran_sys.set_static_friction_coeff_SPH2SPH(params.static_friction_coeffS2S);
    gran_sys.set_static_friction_coeff_SPH2WALL(params.static_friction_coeffS2W);

    gran_sys.set_Cohesion_ratio(params.cohesion_ratio);
    gran_sys.set_Adhesion_ratio_S2W(params.adhesion_ratio_s2w);

    double x = -std::sqrt((2 * std::cos(theta) * std::cos(theta)) / (4 - 4 * std::cos(theta) * std::cos(theta)));

    double ax = x / std::sqrt(2 * x * x + 1);
    double az = 1 / std::sqrt(2 * x * x + 1);

    ChVector<float> n((float)-az / (2 * ax), (float)-az / (2 * ax), 1.f);
    n.Normalize();
    
    float plane_pos[] = {-params.sphere_radius * n.x(), -params.sphere_radius * n.y(), -params.sphere_radius * n.z()};
    float plane_normal[] = {n.x(), n.y(), n.z()};
    bool track_forces = false;
    gran_sys.Create_BC_Plane(plane_pos, plane_normal, track_forces);

    gran_sys.set_gravitational_acceleration(params.grav_X, params.grav_Y, params.grav_Z);
    gran_sys.setOutputMode(params.write_mode);
    gran_sys.setOutputFlags(GRAN_OUTPUT_FLAGS::VEL_COMPONENTS | GRAN_OUTPUT_FLAGS::ANG_VEL_COMPONENTS);

    gran_sys.set_friction_mode(GRAN_FRICTION_MODE::MULTI_STEP);
    gran_sys.set_timeIntegrator(GRAN_TIME_INTEGRATOR::CENTERED_DIFFERENCE);

    std::vector<ChVector<float>> body_points;
    body_points.push_back(ChVector<float>(0, 0, 0));

    apiSMC.setElemsPositions(body_points);

    gran_sys.set_fixed_stepSize(params.step_size);

    filesystem::create_directory(filesystem::path(params.output_dir));

    gran_sys.set_BD_Fixed(true);

    gran_sys.setVerbose(params.verbose);
    gran_sys.initialize();

    // int fps = 10000;
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
        gran_sys.advance_simulation(frame_step);
        curr_time += frame_step;
        printf("rendering frame %u\n", currframe);
        sprintf(filename, "%s/step%06d", params.output_dir.c_str(), currframe++);
        gran_sys.writeFile(std::string(filename));
    }

    return 0;
}
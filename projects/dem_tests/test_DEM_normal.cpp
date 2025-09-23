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

#include "DemDemoUtils.h"
#include "chrono/utils/ChUtilsSamplers.h"
#include "chrono_dem/physics/ChSystemDem.h"
#include "chrono_dem/utils/ChDemJsonParser.h"
#include "chrono_thirdparty/filesystem/path.h"

#include "../utils.h"

using namespace chrono;
using namespace chrono::dem;

enum RUN_MODE { SPHERE_WALL = 0, SPHERE_SPHERE = 1, SPHERE_WALL_ANGLE = 2 };

int main(int argc, char* argv[]) {
    std::string inputJson = GetProjectsDataFile("dem/Normal.json");
    RUN_MODE run_mode = RUN_MODE::SPHERE_WALL;
    if (argc == 2) {
        inputJson = std::string(argv[1]);
    } else if (argc == 3) {
        inputJson = std::string(argv[1]);
        run_mode = (RUN_MODE)std::atoi(argv[2]);
    } else if (argc > 1) {
        std::cout << "Usage:\n./test_DEM_normal <json_file> [<run_mode>]" << std::endl;
        return 1;
    }

    ChDemSimulationParameters params;
    if (!ParseJSON(inputJson, params)) {
        std ::cout << "ERROR: reading input file " << inputJson << std::endl;
        return 1;
    }

    params.box_X = 10;
    params.box_Y = 10;
    params.box_Z = (run_mode == SPHERE_SPHERE) ? 8 * params.sphere_radius : 4 * params.sphere_radius;

    // Setup simulation
    ChSystemDem dem_sys(params.sphere_radius, params.sphere_density,
                        ChVector3f(params.box_X, params.box_Y, params.box_Z));
    dem_sys.DisableMinLength();

    dem_sys.SetPsiFactors(params.psi_T, params.psi_L);

    dem_sys.SetKn_SPH2SPH(params.normalStiffS2S);
    dem_sys.SetKn_SPH2WALL(params.normalStiffS2W);
    dem_sys.SetGn_SPH2SPH(params.normalDampS2S);
    dem_sys.SetGn_SPH2WALL(params.normalDampS2W);

    dem_sys.SetCohesionRatio(params.cohesion_ratio);
    dem_sys.SetAdhesionRatio_SPH2WALL(params.adhesion_ratio_s2w);
    if (run_mode == SPHERE_WALL_ANGLE) {
        params.grav_X = -565.80;
        params.grav_Y = -565.80;
        params.grav_Z = -565.80;

        ChVector3d plane_pos((float)(-2 * params.sphere_radius / std::sqrt(3)),
                             (float)(-2 * params.sphere_radius / std::sqrt(3)),
                             (float)(-2 * params.sphere_radius / std::sqrt(3)));

        ChVector3d plane_normal(1, 1, 1);
        bool track_forces = false;
        dem_sys.CreateBCPlane(plane_pos, plane_normal, track_forces);
    }

    dem_sys.SetGravitationalAcceleration(ChVector3d(params.grav_X, params.grav_Y, params.grav_Z));
    dem_sys.SetParticleOutputMode(params.write_mode);
    dem_sys.SetParticleOutputFlags(
        CHDEM_OUTPUT_FLAGS::VEL_COMPONENTS | CHDEM_OUTPUT_FLAGS::FIXITY |
        CHDEM_OUTPUT_FLAGS::FORCE_COMPONENTS);  // NOTE: original test used custom FORCE_COMPONENTS output

    dem_sys.SetFrictionMode(CHDEM_FRICTION_MODE::FRICTIONLESS);
    dem_sys.SetTimeIntegrator(CHDEM_TIME_INTEGRATOR::CENTERED_DIFFERENCE);

    std::vector<ChVector3f> body_points;
    body_points.push_back(ChVector3f(0, 0, 0));

    if (run_mode == SPHERE_SPHERE) {
        std::vector<bool> body_points_fixed;
        body_points.push_back(ChVector3f(0, 0, -3 * params.sphere_radius));
        body_points_fixed.push_back(false);
        body_points_fixed.push_back(true);
        dem_sys.SetParticleFixed(body_points_fixed);
    }

    dem_sys.SetParticles(body_points);

    dem_sys.SetFixedStepSize(params.step_size);

    filesystem::create_directory(filesystem::path(params.output_dir));

    dem_sys.SetBDFixed(true);

    dem_sys.SetVerbosity(params.verbose);
    dem_sys.Initialize();

    int fps = 10000;
    float frame_step = 1.f / fps;
    float curr_time = 0;
    int currframe = 0;

    // write an initial frame
    char filename[100];
    sprintf(filename, "%s/step%06d", params.output_dir.c_str(), currframe++);
    dem_sys.WriteParticleFile(std::string(filename));

    std::cout << "frame step is " << frame_step << std::endl;
    while (curr_time < params.time_end) {
        dem_sys.AdvanceSimulation(frame_step);
        curr_time += frame_step;
        printf("rendering frame %u\n", currframe);
        sprintf(filename, "%s/step%06d", params.output_dir.c_str(), currframe++);
        dem_sys.WriteParticleFile(std::string(filename));
    }

    return 0;
}

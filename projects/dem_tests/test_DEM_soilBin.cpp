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
// Authors: Conlain Kelly, Nic Olsen
// =============================================================================
// Simple Chrono::Dem settling experiment which allows for sweeping various
// simulation parameters to produce scaling analyses.
// =============================================================================

#include <iostream>
#include <string>
#include <cmath>
#include <chrono>

#include "chrono/core/ChGlobal.h"
#include "chrono/utils/ChUtilsSamplers.h"

#include "chrono_dem/physics/ChSystemDem.h"
#include "chrono_dem/utils/ChDemJsonParser.h"

#include "chrono_thirdparty/filesystem/path.h"

#include "../utils.h"

using namespace chrono;
using namespace chrono::dem;

int main(int argc, char* argv[]) {
    std::string inputJson = GetProjectsDataFile("dem/soilBin.json");
    if (argc == 2) {
        inputJson = std::string(argv[1]);
    } else if (argc > 2) {
        std::cout << "Usage:\n./demo_DEM_soilBin <json_file>" << std::endl;
        return 1;
    }

    ChDemSimulationParameters params;
    if (!ParseJSON(inputJson, params)) {
        std ::cout << "ERROR: reading input file " << inputJson << std::endl;
        return 1;
    }

    // Setup simulation
    ChSystemDem dem_sys(params.sphere_radius, params.sphere_density,
                        ChVector3f(params.box_X, params.box_Y, params.box_Z));

    dem_sys.SetPsiFactors(params.psi_T, params.psi_L);

    dem_sys.SetKn_SPH2SPH(params.normalStiffS2S);
    dem_sys.SetKn_SPH2WALL(params.normalStiffS2W);
    dem_sys.SetGn_SPH2SPH(params.normalDampS2S);
    dem_sys.SetGn_SPH2WALL(params.normalDampS2W);

    dem_sys.SetKt_SPH2SPH(params.tangentStiffS2S);
    dem_sys.SetKt_SPH2WALL(params.tangentStiffS2W);
    dem_sys.SetGt_SPH2SPH(params.tangentDampS2S);
    dem_sys.SetGt_SPH2WALL(params.tangentDampS2W);

    dem_sys.SetStaticFrictionCoeff_SPH2SPH(params.static_friction_coeffS2S);
    dem_sys.SetStaticFrictionCoeff_SPH2WALL(params.static_friction_coeffS2W);

    dem_sys.SetCohesionRatio(params.cohesion_ratio);
    dem_sys.SetAdhesionRatio_SPH2WALL(params.adhesion_ratio_s2w);

    dem_sys.SetGravitationalAcceleration(ChVector3f(params.grav_X, params.grav_Y, params.grav_Z));
    dem_sys.SetParticleOutputMode(params.write_mode);

    dem_sys.SetRollingMode(CHDEM_ROLLING_MODE::NO_RESISTANCE);

    std::vector<ChVector3f> body_points;

    {
        // fill box, layer by layer
        ChVector3d hdims(params.box_X / 2.f - 2 * params.sphere_radius, params.box_Y / 2.f - 2 * params.sphere_radius,
                         params.box_Z / 2.f - 2 * params.sphere_radius);
        ChVector3d center(0, 0, 0);

        utils::ChHCPSampler<float> sampler(2.2f * params.sphere_radius);

        body_points = sampler.SampleBox(center, hdims);
    }

    dem_sys.SetParticles(body_points);
    std::cout << "Added " << body_points.size() << std::endl;

    switch (params.run_mode) {
        case CHDEM_RUN_MODE::MULTI_STEP:
            dem_sys.SetFrictionMode(CHDEM_FRICTION_MODE::MULTI_STEP);
            dem_sys.SetKt_SPH2SPH(params.tangentStiffS2S);
            dem_sys.SetKt_SPH2WALL(params.tangentStiffS2W);
            dem_sys.SetGt_SPH2SPH(params.tangentDampS2S);
            dem_sys.SetGt_SPH2WALL(params.tangentDampS2W);
            break;
        case CHDEM_RUN_MODE::ONE_STEP:
            dem_sys.SetFrictionMode(CHDEM_FRICTION_MODE::SINGLE_STEP);
            dem_sys.SetKt_SPH2SPH(params.tangentStiffS2S);
            dem_sys.SetKt_SPH2WALL(params.tangentStiffS2W);
            dem_sys.SetGt_SPH2SPH(params.tangentDampS2S);
            dem_sys.SetGt_SPH2WALL(params.tangentDampS2W);
            break;
        case CHDEM_RUN_MODE::FRICTIONLESS:
            dem_sys.SetFrictionMode(CHDEM_FRICTION_MODE::FRICTIONLESS);
            break;
        default:
            std::cout << "Invalid run mode" << std::endl;
            return 1;
    }

    dem_sys.SetTimeIntegrator(CHDEM_TIME_INTEGRATOR::EXTENDED_TAYLOR);
    dem_sys.SetFixedStepSize(params.step_size);

    std::string out_dir;
    if (params.write_mode != CHDEM_OUTPUT_MODE::NONE) {
        out_dir = GetChronoOutputPath() + "DEM/";
        filesystem::create_directory(filesystem::path(out_dir));
        out_dir = out_dir + params.output_dir;
        filesystem::create_directory(filesystem::path(out_dir));
    }
    dem_sys.SetBDFixed(true);

    dem_sys.SetVerbosity(params.verbose);
    dem_sys.Initialize();

    int fps = 50;
    float frame_step = 1.f / fps;
    float curr_time = 0;
    int currframe = 0;
    unsigned int total_frames = (unsigned int)((float)params.time_end * fps);

    // write an initial frame
    char filename[100];
    sprintf(filename, "%s/step%06d.csv", out_dir.c_str(), currframe++);
    dem_sys.WriteParticleFile(std::string(filename));

    std::cout << "frame step is " << frame_step << std::endl;

    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    while (curr_time < params.time_end) {
        dem_sys.AdvanceSimulation(frame_step);
        curr_time += frame_step;
        printf("rendering frame %u of %u\n", currframe, total_frames + 1);
        sprintf(filename, "%s/step%06d.csv", out_dir.c_str(), currframe++);
        dem_sys.WriteParticleFile(std::string(filename));
    }
    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_sec = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    std::cout << time_sec.count() << " seconds" << std::endl;

    return 0;
}

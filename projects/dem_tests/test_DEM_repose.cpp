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
// A column of granular material forms a mound
// =============================================================================

#include <iostream>
#include <string>

#include "chrono/core/ChGlobal.h"
#include "chrono/utils/ChUtilsSamplers.h"
#include "chrono_dem/physics/ChSystemDem.h"
#include "chrono_dem/utils/ChDemJsonParser.h"
#include "chrono_dem/visualization/ChDemVisualizationGL.h"
#include "chrono_thirdparty/filesystem/path.h"

#include "../utils.h"

using namespace chrono;
using namespace chrono::dem;

// Enable/disable run-time visualization (if Chrono::OpenGL is available)
bool render = true;

int main(int argc, char* argv[]) {
    std::string inputJson = GetProjectsDataFile("dem/Repose.json");
    if (argc == 2) {
        inputJson = std::string(argv[1]);
    } else if (argc > 1) {
        std::cout << "Usage:\n./test_DEM_repose <json_file>" << std::endl;
        return 1;
    }

    ChDemSimulationParameters params;
    if (!ParseJSON(inputJson, params)) {
        std ::cout << "ERROR: reading input file " << inputJson << std::endl;
        return 1;
    }

    std::string out_dir = GetChronoOutputPath() + "DEM/";
    filesystem::create_directory(filesystem::path(out_dir));
    out_dir = out_dir + params.output_dir;
    filesystem::create_directory(filesystem::path(out_dir));

    // Setup simulation
    ChSystemDem dem_sys(params.sphere_radius, params.sphere_density,
                        ChVector3f(params.box_X, params.box_Y, params.box_Z));

    dem_sys.SetKn_SPH2SPH(params.normalStiffS2S);
    dem_sys.SetKn_SPH2WALL(params.normalStiffS2W);
    dem_sys.SetGn_SPH2SPH(params.normalDampS2S);
    dem_sys.SetGn_SPH2WALL(params.normalDampS2W);

    // dem_sys.SetFrictionMode(CHDEM_FRICTION_MODE::FRICTIONLESS);
    dem_sys.SetFrictionMode(CHDEM_FRICTION_MODE::MULTI_STEP);
    dem_sys.SetKt_SPH2SPH(params.tangentStiffS2S);
    dem_sys.SetKt_SPH2WALL(params.tangentStiffS2W);
    dem_sys.SetGt_SPH2SPH(params.tangentDampS2S);
    dem_sys.SetGt_SPH2WALL(params.tangentDampS2W);
    dem_sys.SetStaticFrictionCoeff_SPH2SPH(params.static_friction_coeffS2S);
    dem_sys.SetStaticFrictionCoeff_SPH2WALL(params.static_friction_coeffS2W);

    // dem_sys.SetRollingMode(CHDEM_ROLLING_MODE::NO_RESISTANCE);
    dem_sys.SetRollingMode(CHDEM_ROLLING_MODE::SCHWARTZ);
    dem_sys.SetRollingCoeff_SPH2SPH(params.rolling_friction_coeffS2S);
    dem_sys.SetRollingCoeff_SPH2WALL(params.rolling_friction_coeffS2W);

    dem_sys.SetCohesionRatio(params.cohesion_ratio);
    dem_sys.SetAdhesionRatio_SPH2WALL(params.adhesion_ratio_s2w);
    dem_sys.SetGravitationalAcceleration(ChVector3f(params.grav_X, params.grav_Y, params.grav_Z));
    dem_sys.SetParticleOutputMode(params.write_mode);

    dem_sys.SetBDFixed(true);

    // padding in sampler
    float fill_epsilon = 2.02f;
    // padding at top of fill
    ////float drop_height = 0.f;
    float spacing = fill_epsilon * params.sphere_radius;
    chrono::utils::ChPDSampler<float> sampler(spacing);

    // Fixed points on the bottom for roughness
    float bottom_z = -params.box_Z / 2.f + params.sphere_radius;
    ChVector3d bottom_center(0, 0, bottom_z);
    std::vector<ChVector3f> roughness_points = sampler.SampleBox(
        bottom_center,
        ChVector3f(params.box_X / 2.f - params.sphere_radius, params.box_Y / 2.f - params.sphere_radius, 0.f));

    // Create column of material
    std::vector<ChVector3f> material_points;

    float fill_bottom = bottom_z + spacing;
    float fill_width = 5.f;
    float fill_height = 2.f * fill_width;
    ////float fill_top = fill_bottom + fill_height;

    ChVector3f center(0.f, 0.f, fill_bottom + fill_height / 2.f);
    material_points = sampler.SampleCylinderZ(center, fill_width, fill_height / 2.f);

    std::vector<ChVector3f> body_points;
    std::vector<bool> body_points_fixed;
    body_points.insert(body_points.end(), roughness_points.begin(), roughness_points.end());
    body_points_fixed.insert(body_points_fixed.end(), roughness_points.size(), true);

    body_points.insert(body_points.end(), material_points.begin(), material_points.end());
    body_points_fixed.insert(body_points_fixed.end(), material_points.size(), false);

    dem_sys.SetParticles(body_points);
    dem_sys.SetParticleFixed(body_points_fixed);

    std::cout << "Added " << roughness_points.size() << " fixed points" << std::endl;
    std::cout << "Added " << material_points.size() << " material points" << std::endl;

    std::cout << "Actually added " << body_points.size() << std::endl;

    dem_sys.SetTimeIntegrator(CHDEM_TIME_INTEGRATOR::EXTENDED_TAYLOR);
    dem_sys.SetFixedStepSize(params.step_size);

    dem_sys.SetVerbosity(params.verbose);
    std::cout << "verbose: " << static_cast<int>(params.verbose) << std::endl;
    dem_sys.SetRecordingContactInfo(true);

    dem_sys.Initialize();

    ChDemVisualizationGL dem_vis(&dem_sys);
    if (render) {
        dem_vis.SetTitle("Chrono::Dem repose demo");
        dem_vis.UpdateCamera(ChVector3d(0, -30, -10), ChVector3d(0, 0, -20));
        dem_vis.SetCameraMoveScale(1.0f);
        dem_vis.Initialize();
    }

    int fps = 60;
    float frame_step = 1.f / fps;
    float curr_time = 0.f;
    int currframe = 0;
    unsigned int total_frames = (unsigned int)((float)params.time_end * fps);

    // write an initial frame
    char filename[100];
    sprintf(filename, "%s/step%06d.csv", out_dir.c_str(), currframe);
    dem_sys.WriteParticleFile(std::string(filename));

    char contactFilename[100];
    sprintf(contactFilename, "%s/contact%06d.csv", out_dir.c_str(), currframe);
    dem_sys.WriteContactInfoFile(std::string(contactFilename));

    currframe++;

    std::cout << "frame step is " << frame_step << std::endl;
    while (curr_time < params.time_end) {
        dem_sys.AdvanceSimulation(frame_step);

        if (render && dem_vis.Render())
            break;

        printf("Output frame %u of %u\n", currframe, total_frames);
        sprintf(filename, "%s/step%06d.csv", out_dir.c_str(), currframe);
        dem_sys.WriteParticleFile(std::string(filename));

        sprintf(contactFilename, "%s/contact%06d.csv", out_dir.c_str(), currframe);
        dem_sys.WriteContactInfoFile(std::string(contactFilename));

        curr_time += frame_step;
        currframe++;
    }

    return 0;
}

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
// Authors: Luning Fang
// =============================================================================
// Chrono::Dem simulation of granular material settled in cylinder first, then
// compressed from a plate on top modelled as a boundary condition
// =============================================================================

#include <iostream>
#include <string>
#include <cmath>

#include "chrono/core/ChGlobal.h"
#include "chrono/core/ChVector3.h"
#include "chrono/utils/ChUtilsSamplers.h"

#include "chrono_thirdparty/filesystem/path.h"

#include "chrono_dem/physics/ChSystemDem.h"
#include "chrono_dem/utils/ChDemJsonParser.h"

#include "../utils.h"

using namespace chrono;
using namespace chrono::dem;

// unit conversion from cgs to si
float F_CGS_TO_SI = 1e-5f;
float KE_CGS_TO_SI = 1e-7f;
float L_CGS_TO_SI = 1e-2f;

int main(int argc, char* argv[]) {
    std::string inputJson = GetProjectsDataFile("dem/compression.json");
    if (argc == 2) {
        inputJson = std::string(argv[1]);
    } else if (argc > 2) {
        std::cout << "Usage:\n./demo_DEM_compression <json_file>" << std::endl;
        return 1;
    }

    ChDemSimulationParameters params;
    if (!ParseJSON(inputJson, params)) {
        std ::cout << "ERROR: reading input file " << inputJson << std::endl;
        return 1;
    }

    // Setup simulation, big domain: 10 by 10 by 20
    ChSystemDem dem_sys(params.sphere_radius, params.sphere_density,
                        ChVector3f(params.box_X, params.box_Y, params.box_Z));

    // One thing we can do is to move the Big Box Domain by (X/2, Y/2, Z/2) using SetBDCenter, so the
    // coordinate range we are now working with is (0,0,0) to (X,Y,Z), instead of (-X/2,-Y/2,-Z/2) to (X/2, Y/2, Z/2).
    dem_sys.SetBDCenter(ChVector3f(params.box_X / 2, params.box_Y / 2, params.box_Z / 2));

    // creat cylinder boundary of Radius 5 at the center of the box domain
    ChVector3f cyl_center(params.box_X / 2, params.box_Y / 2, params.box_Z / 2);
    float cyl_rad = std::min(params.box_X, params.box_Y) / 2.0f;
    dem_sys.CreateBCCylinderZ(cyl_center, cyl_rad, false, true);

    // initialize sampler, set distance between center of spheres as 2.1r
    utils::ChHCPSampler<float> sampler(2.1f * params.sphere_radius);
    std::vector<ChVector3f> initialPos;

    // randomize by layer
    ChVector3f center(params.box_X / 2, params.box_Y / 2, params.sphere_radius);
    // fill up each layer
    while (center.z() + params.sphere_radius < params.box_Z) {
        auto points = sampler.SampleCylinderZ(center, cyl_rad - params.sphere_radius, 0);
        initialPos.insert(initialPos.end(), points.begin(), points.end());
        center.z() += 2.1f * params.sphere_radius;
    }

    size_t numSpheres = initialPos.size();

    // create initial velocity vector
    std::vector<ChVector3f> initialVelo;
    for (size_t i = 0; i < numSpheres; i++) {
        ChVector3f velo(-initialPos.at(i).x() / cyl_rad, -initialPos.at(i).x() / cyl_rad, 0.0f);
        initialVelo.push_back(velo);
    }

    // assign initial position and velocity to the granular system
    dem_sys.SetParticles(initialPos, initialVelo);

    dem_sys.SetPsiFactors(params.psi_T, params.psi_L);

    // normal force model
    dem_sys.SetKn_SPH2SPH(params.normalStiffS2S);
    dem_sys.SetKn_SPH2WALL(params.normalStiffS2W);
    dem_sys.SetGn_SPH2SPH(params.normalDampS2S);
    dem_sys.SetGn_SPH2WALL(params.normalDampS2W);

    // assign tangential force model and its parameters
    dem_sys.SetFrictionMode(CHDEM_FRICTION_MODE::MULTI_STEP);
    dem_sys.SetKt_SPH2SPH(params.tangentStiffS2S);
    dem_sys.SetKt_SPH2WALL(params.tangentStiffS2W);
    dem_sys.SetGt_SPH2SPH(params.tangentDampS2S);
    dem_sys.SetGt_SPH2WALL(params.tangentDampS2W);

    dem_sys.SetStaticFrictionCoeff_SPH2SPH(params.static_friction_coeffS2W);
    dem_sys.SetStaticFrictionCoeff_SPH2WALL(params.static_friction_coeffS2W);

    dem_sys.SetCohesionRatio(params.cohesion_ratio);
    dem_sys.SetAdhesionRatio_SPH2WALL(params.adhesion_ratio_s2w);

    dem_sys.SetGravitationalAcceleration(ChVector3f(params.grav_X, params.grav_Y, params.grav_Z));
    dem_sys.SetParticleOutputMode(params.write_mode);

    std::string out_dir = GetChronoOutputPath() + "DEM/";
    filesystem::create_directory(filesystem::path(out_dir));
    out_dir = out_dir + params.output_dir;
    filesystem::create_directory(filesystem::path(out_dir));

    // Set the position of the BD fixed
    dem_sys.SetBDFixed(true);
    dem_sys.SetTimeIntegrator(CHDEM_TIME_INTEGRATOR::FORWARD_EULER);
    dem_sys.SetFixedStepSize(params.step_size);

    dem_sys.SetVerbosity(params.verbose);

    // create top plane boundary condition with its position and normal
    ChVector3f topWallPos(params.box_X / 2, params.box_Y / 2, params.box_Z);
    ChVector3f topWallNrm(0.0f, 0.0f, -1.0f);
    size_t topWall = dem_sys.CreateBCPlane(topWallPos, topWallNrm, true);

    float topWall_vel;       // top plane moving velocity
    float topWall_offset;    // position offset of the plane when it first starts to move
    float topWall_moveTime;  // time when the plane first starts to move

    // user defined offset position function for top wall
    std::function<double3(float)> topWall_posFunc = [&topWall_offset, &topWall_vel, &topWall_moveTime](float t) {
        double3 pos = {0, 0, 0};
        pos.z = topWall_offset + topWall_vel * (t - topWall_moveTime);
        return pos;
    };

    dem_sys.SetParticleOutputFlags(ABSV);
    dem_sys.Initialize();

    // output frames per second
    int fps = 100;
    // assume we run for at least one frame
    float frame_step = 1.0f / fps;
    float curr_time = 0;
    int curr_frame = 0;
    unsigned int total_frames = (unsigned int)(((float)params.time_end - 0.5f) * fps) - 1;
    // initialize values that I want to keep track of
    ChVector3f plane_reaction_force;
    ChVector3f platePos;
    int nc;

    // let system run for 0.5 second so the particles can settle
    while (curr_time < 0.5) {
        dem_sys.AdvanceSimulation(frame_step);
        curr_time += frame_step;
        printf("time = %.4f\n", curr_time);
    }

    // top plate move downward with velocity 1cm/s
    topWall_vel = -1.0f;
    // i would like it to start from the top most sphere
    topWall_offset = (float)dem_sys.GetMaxParticleZ() + params.sphere_radius - topWallPos[2];
    topWall_moveTime = curr_time;

    // sphere settled now push the plate downward
    dem_sys.SetBCOffsetFunction(topWall, topWall_posFunc);

    // continue simulation until the end
    while (curr_time < params.time_end) {
        printf("rendering frame: %u of %u, curr_time: %.4f, ", curr_frame + 1, total_frames, curr_time);

        // write position
        char filename[100];
        sprintf(filename, "%s/step%06d.csv", out_dir.c_str(), curr_frame);
        dem_sys.WriteParticleFile(std::string(filename));
        dem_sys.AdvanceSimulation(frame_step);

        platePos = dem_sys.GetBCPlanePosition(topWall);
        std::cout << "top plate pos_z: " << platePos.z() << " cm";

        nc = dem_sys.GetNumContacts();
        std::cout << ", numContacts: " << nc;

        dem_sys.GetBCReactionForces(topWall, plane_reaction_force);
        std::cout << ", top plate force: " << plane_reaction_force.z() * F_CGS_TO_SI << " Newton";
        std::cout << "\n";

        curr_frame++;
        curr_time += frame_step;
    }
    return 0;
}

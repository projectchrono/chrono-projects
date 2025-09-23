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
// Authors: Conlain Kelly
// =============================================================================
// Chrono::Granular simulation of material flowing out of a cylindrical hopper.
// =============================================================================

#include <iostream>
#include <string>

#include "chrono/core/ChGlobal.h"
#include "chrono/utils/ChUtilsSamplers.h"
#include "chrono_dem/physics/ChSystemDem.h"
#include "chrono_dem/utils/ChDemJsonParser.h"
#include "chrono_thirdparty/filesystem/path.h"

#include "../utils.h"

using namespace chrono;
using namespace chrono::dem;

std::string cyl_filename = GetProjectsDataFile("dem/meshes/Gran_cylinder_transparent.obj");

// Take a ChBody and write its
void writeZCylinderMesh(std::ostringstream& outstream, ChVector3d pos, float rad, float height) {
    // Get basis vectors
    ChVector3d vx(1, 0, 0);
    ChVector3d vy(0, 1, 0);
    ChVector3d vz(0, 0, 1);

    ChVector3d scaling(rad, rad, height / 2);

    // Write the mesh name to find
    outstream << cyl_filename << ",";
    // Output in order
    outstream << pos.x() << ",";
    outstream << pos.y() << ",";
    outstream << pos.z() << ",";
    outstream << vx.x() << ",";
    outstream << vx.y() << ",";
    outstream << vx.z() << ",";
    outstream << vy.x() << ",";
    outstream << vy.y() << ",";
    outstream << vy.z() << ",";
    outstream << vz.x() << ",";
    outstream << vz.y() << ",";
    outstream << vz.z() << ",";
    outstream << scaling.x() << ",";
    outstream << scaling.y() << ",";
    outstream << scaling.z();
    outstream << "\n";
}

// Take a ChBody and write its
void writeZConeMesh(std::ostringstream& outstream, ChVector3d pos, std::string mesh_filename) {
    // Get basis vectors
    ChVector3d vx(1, 0, 0);
    ChVector3d vy(0, 1, 0);
    ChVector3d vz(0, 0, 1);

    ChVector3d scaling(1, 1, 1);

    // Write the mesh name to find
    outstream << mesh_filename << ",";
    // Output in order
    outstream << pos.x() << ",";
    outstream << pos.y() << ",";
    outstream << pos.z() << ",";
    outstream << vx.x() << ",";
    outstream << vx.y() << ",";
    outstream << vx.z() << ",";
    outstream << vy.x() << ",";
    outstream << vy.y() << ",";
    outstream << vy.z() << ",";
    outstream << vz.x() << ",";
    outstream << vz.y() << ",";
    outstream << vz.z() << ",";
    outstream << scaling.x() << ",";
    outstream << scaling.y() << ",";
    outstream << scaling.z();
    outstream << "\n";
}

int main(int argc, char* argv[]) {
    std::string inputJson = GetProjectsDataFile("dem/Coneflow.json");
    float aperture_diameter = 16.f;
    if (argc == 2) {
        inputJson = std::string(argv[1]);
    } else if (argc == 3) {
        inputJson = std::string(argv[1]);
        aperture_diameter = (float)std::atof(argv[2]);
    } else if (argc > 1) {
        std::cout << "Usage:\n./test_DEM_Coneflow <json_file> [<aperture_diameter>]" << std::endl;
        return 1;
    }

    ChDemSimulationParameters params;
    if (!ParseJSON(inputJson, params)) {
        std ::cout << "ERROR: reading input file " << inputJson << std::endl;
        return 1;
    }

    // Setup simulation
    ChSystemDem gran_sys(params.sphere_radius, params.sphere_density,
                         ChVector3f(params.box_X, params.box_Y, params.box_Z));
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

    gran_sys.SetCohesionRatio(params.cohesion_ratio);
    gran_sys.SetAdhesionRatio_SPH2WALL(params.adhesion_ratio_s2w);
    gran_sys.SetGravitationalAcceleration(ChVector3f(params.grav_X, params.grav_Y, params.grav_Z));
    gran_sys.SetParticleOutputMode(params.write_mode);

    gran_sys.SetBDFixed(true);

    // Fill box with bodies
    std::vector<ChVector3f> body_points;

    // padding in sampler
    float fill_epsilon = 2.02f;
    // padding at top of fill
    float fill_gap = 1.f;

    chrono::utils::ChPDSampler<float> sampler(fill_epsilon * params.sphere_radius);

    ChVector3f center_pt(0, 0, -2 - params.box_Z / 6);

    // width we want to fill to
    float fill_width = params.box_Z / 3.f;
    // height that makes this width above the cone
    float fill_height = fill_width;

    // fill to top
    float fill_top = params.box_Z / 2 - fill_gap;
    float fill_bottom = fill_top - fill_height;

    printf("width is %f, bot is %f, top is %f, height is %f\n", fill_width, fill_bottom, fill_top, fill_height);
    // fill box, layer by layer
    ChVector3d center(0, 0, fill_bottom);
    // shift up for bottom of box
    center.z() += fill_gap;

    while (center.z() < fill_top) {
        std::cout << "Create layer at " << center.z() << std::endl;
        auto points = sampler.SampleCylinderZ(center, fill_width, 0);
        body_points.insert(body_points.end(), points.begin(), points.end());
        center.z() += fill_epsilon * params.sphere_radius;
    }

    std::vector<ChVector3f> body_points_first;
    body_points_first.push_back(body_points[0]);

    gran_sys.SetParticles(body_points);

    float sphere_mass =
        (4.f / 3.f) * params.sphere_density * params.sphere_radius * params.sphere_radius * params.sphere_radius;

    printf("%d spheres with mass %f \n", (int)body_points.size(), body_points.size() * sphere_mass);

    // set time integrator
    gran_sys.SetTimeIntegrator(CHDEM_TIME_INTEGRATOR::CENTERED_DIFFERENCE);
    gran_sys.SetFixedStepSize(params.step_size);

    // set friction mode
    gran_sys.SetFrictionMode(CHDEM_FRICTION_MODE::MULTI_STEP);

    filesystem::create_directory(filesystem::path(params.output_dir));

    constexpr float cone_slope = 1.0;

    float cone_offset = aperture_diameter / 2.f;

    gran_sys.SetVerbosity(params.verbose);
    float hmax = params.box_Z;
    float hmin = center_pt.z() + cone_offset;
    // Finalize settings and initialize for runtime
    gran_sys.CreateBCConeZ(center_pt, cone_slope, hmax, hmin, false, false);

    ChVector3d cone_top_pos(0, 0, center_pt.z() + fill_width + 8);

    float cyl_rad = fill_width + 8;
    printf("top of cone is at %f, cone tip is %f, top width is %f, bottom width is hmin %f\n", cone_top_pos.z(),
           fill_width + 8, hmax, cone_offset);

    ChVector3f zvec(0, 0, 0);
    {
        std::string meshes_file = "coneflow_meshes.csv";

        std::ofstream meshfile{params.output_dir + "/" + meshes_file};
        std::ostringstream outstream;
        outstream << "mesh_name,dx,dy,dz,x1,x2,x3,y1,y2,y3,z1,z2,z3\n";
        writeZConeMesh(outstream, cone_top_pos, GetProjectsDataFile("dem/meshes/gran_zcone.obj"));
        writeZCylinderMesh(outstream, zvec, cyl_rad, params.box_Z);

        meshfile << outstream.str();
    }

    gran_sys.CreateBCCylinderZ(zvec, cyl_rad, false, false);

    ChVector3f plane_center(0, 0, center_pt.z() + 2 * cone_slope + cone_slope * cone_offset);
    ChVector3f plane_normal(0, 0, 1);

    printf("center is %f, %f, %f, plane center is is %f, %f, %f\n", center_pt.x(), center_pt.y(), center_pt.z(),
           plane_center.x(), plane_center.y(), plane_center.z());
    size_t cone_plane_bc_id = gran_sys.CreateBCPlane(plane_center, plane_normal, false);

    // put a plane at the bottom of the box to count forces
    ChVector3f box_bottom(0, 0, -params.box_Z / 2.f + 2.f);

    size_t bottom_plane_bc_id = gran_sys.CreateBCPlane(box_bottom, plane_normal, true);

    gran_sys.Initialize();

    // number of times to capture force data per second
    int captures_per_second = 200;
    // number of times to capture force before we capture a frame
    int captures_per_frame = 4;

    // assume we run for at least one frame
    float frame_step = 1.0f / captures_per_second;
    float curr_time = 0;
    int currcapture = 0;
    int currframe = 0;

    std::cout << "capture step is " << frame_step << std::endl;

    float t_remove_plane = .5;
    bool plane_active = false;

    ChVector3f reaction_forces(0, 0, 0);

    constexpr float F_CGS_TO_SI = 1e-5f;
    constexpr float M_CGS_TO_SI = 1e-3f;
    float total_system_mass = 4.0f / 3.0f * (float)CH_PI * params.sphere_density * params.sphere_radius * params.sphere_radius *
                              params.sphere_radius * body_points.size();
    printf("total system mass is %f kg \n", total_system_mass * M_CGS_TO_SI);
    char filename[100];
    // sprintf(filename, "%s/step%06d", params.output_dir.c_str(), currframe++);
    // gran_sys.writeFile(std::string(filename));

    // Run settling experiments
    while (curr_time < params.time_end) {
        if (!plane_active && curr_time > t_remove_plane) {
            gran_sys.DisableBCbyID(cone_plane_bc_id);
        }

        bool success = gran_sys.GetBCReactionForces(bottom_plane_bc_id, reaction_forces);
        if (!success) {
            printf("ERROR! Get contact forces for plane failed\n");
        } else {
            printf("curr time is %f, plane force is (%f, %f, %f) Newtons\n", curr_time,
                   F_CGS_TO_SI * reaction_forces.x(), F_CGS_TO_SI * reaction_forces.y(),
                   F_CGS_TO_SI * reaction_forces.z());
        }
        gran_sys.AdvanceSimulation(frame_step);
        curr_time += frame_step;

        // if this frame is a render frame
        if (currcapture % captures_per_frame == 0) {
            printf("rendering frame %u\n", currframe);
            sprintf(filename, "%s/step%06d", params.output_dir.c_str(), currframe++);
            gran_sys.WriteParticleFile(std::string(filename));
        }
        currcapture++;
    }

    return 0;
}

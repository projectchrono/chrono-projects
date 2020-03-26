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
// Authors: Nic Olsen, Conlain Kelly
// =============================================================================
// Granular material of varying material properties is poured onto a rough
// surface to form a mound.
// =============================================================================

#include <iostream>
#include <string>

#include "ChGranularDemoUtils.hpp"
#include "chrono/core/ChGlobal.h"
#include "chrono/utils/ChUtilsSamplers.h"
#include "chrono_granular/ChGranularData.h"
#include "chrono_granular/api/ChApiGranularChrono.h"
#include "chrono_granular/physics/ChGranular.h"
#include "chrono_granular/utils/ChGranularJsonParser.h"
#include "chrono_thirdparty/filesystem/path.h"

using namespace chrono;
using namespace chrono::granular;

int num_args = 5;
void ShowUsage(std::string name) {
    std::cout << "usage: " + name + " <json_file> <static_friction> <rolling_friction> <output_dir>" << std::endl;
}

std::string cyl_filename = granular::GetDataFile("shared/Gran_cylinder_transparent.obj");

void writeZCylinderMesh(std::ostringstream& outstream, ChVector<> pos, float rad, float height) {
    // Get basis vectors
    ChVector<> vx(1, 0, 0);
    ChVector<> vy(0, 1, 0);
    ChVector<> vz(0, 0, 1);

    ChVector<> scaling(rad, rad, height / 2);

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

void writeZConeMesh(std::ostringstream& outstream, ChVector<> pos, std::string mesh_filename) {
    // Get basis vectors
    ChVector<> vx(1, 0, 0);
    ChVector<> vy(0, 1, 0);
    ChVector<> vz(0, 0, 1);

    ChVector<> scaling(1, 1, 1);

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
    granular::SetDataPath(std::string(PROJECTS_DATA_DIR) + "granular/");
    sim_param_holder params;

    // Some of the default values might be overwritten by user via command line
    if (argc != num_args || ParseJSON(argv[1], params) == false) {
        ShowUsage(argv[0]);
        return 1;
    }

    float aperture_diameter = 6.f;

    params.static_friction_coeffS2S = std::stof(argv[2]);
    params.static_friction_coeffS2W = std::stof(argv[2]);

    params.rolling_friction_coeffS2S = std::stof(argv[3]);
    params.rolling_friction_coeffS2W = std::stof(argv[3]);

    params.output_dir = std::string(argv[4]);

    // Setup simulation
    ChSystemGranularSMC gran_sys(params.sphere_radius, params.sphere_density,
                                 make_float3(params.box_X, params.box_Y, params.box_Z));
    gran_sys.set_K_n_SPH2SPH(params.normalStiffS2S);
    gran_sys.set_K_n_SPH2WALL(params.normalStiffS2W);
    gran_sys.set_Gamma_n_SPH2SPH(params.normalDampS2S);
    gran_sys.set_Gamma_n_SPH2WALL(params.normalDampS2W);

    gran_sys.set_friction_mode(GRAN_FRICTION_MODE::MULTI_STEP);
    gran_sys.set_K_t_SPH2SPH(params.tangentStiffS2S);
    gran_sys.set_K_t_SPH2WALL(params.tangentStiffS2W);
    gran_sys.set_Gamma_t_SPH2SPH(params.tangentDampS2S);
    gran_sys.set_Gamma_t_SPH2WALL(params.tangentDampS2W);

    gran_sys.set_rolling_mode(GRAN_ROLLING_MODE::SCHWARTZ);
    gran_sys.set_rolling_coeff_SPH2SPH(params.rolling_friction_coeffS2S);
    gran_sys.set_rolling_coeff_SPH2WALL(params.rolling_friction_coeffS2W);

    gran_sys.set_Cohesion_ratio(params.cohesion_ratio);
    gran_sys.set_Adhesion_ratio_S2W(params.adhesion_ratio_s2w);
    gran_sys.set_gravitational_acceleration(params.grav_X, params.grav_Y, params.grav_Z);
    gran_sys.setOutputMode(params.write_mode);

    gran_sys.set_BD_Fixed(true);

    // padding in sampler
    float fill_epsilon = 2.02f;
    // padding at top of fill
    float fill_gap = 1.f;
    chrono::utils::PDSampler<float> sampler(fill_epsilon * params.sphere_radius);

    // Fill box with bodies
    std::vector<ChVector<float>> material_points;

    // width we want to fill to
    float fill_width = params.box_Z / 3.f;
    // height that makes this width above the cone
    float fill_height = fill_width;

    // fill to top
    float fill_top = params.box_Z / 2 - fill_gap;
    float fill_bottom = fill_top - fill_height;

    // fill box, layer by layer
    ChVector<> center(0, 0, fill_bottom);
    // shift up for bottom of box
    center.z() += fill_gap;

    while (center.z() < fill_top) {
        std::cout << "Create layer at " << center.z() << std::endl;
        auto points = sampler.SampleCylinderZ(center, fill_width, 0);
        material_points.insert(material_points.end(), points.begin(), points.end());
        center.z() += fill_epsilon * params.sphere_radius;
    }

    // Fixed points on the bottom for roughness
    float bottom_z = -params.box_Z / 2.f + params.sphere_radius;
    ChVector<> bottom_center(0, 0, bottom_z);
    std::vector<ChVector<float>> roughness_points =
        sampler.SampleCylinderZ(bottom_center, params.box_Z / 2.f - params.sphere_radius, 0);

    ChGranularSMC_API apiSMC;
    apiSMC.setGranSystem(&gran_sys);

    std::vector<ChVector<float>> body_points;
    std::vector<bool> body_points_fixed;
    body_points.insert(body_points.end(), roughness_points.begin(), roughness_points.end());
    body_points_fixed.insert(body_points_fixed.end(), roughness_points.size(), true);

    body_points.insert(body_points.end(), material_points.begin(), material_points.end());
    body_points_fixed.insert(body_points_fixed.end(), material_points.size(), false);

    apiSMC.setElemsPositions(body_points);
    gran_sys.setParticleFixed(body_points_fixed);

    std::cout << "Added " << roughness_points.size() << " fixed points" << std::endl;
    std::cout << "Added " << material_points.size() << " material points" << std::endl;

    gran_sys.set_timeIntegrator(GRAN_TIME_INTEGRATOR::EXTENDED_TAYLOR);
    gran_sys.set_fixed_stepSize(params.step_size);

    filesystem::create_directory(filesystem::path(params.output_dir));
    gran_sys.setVerbose(params.verbose);

    constexpr float cone_slope = 1.0;

    float cone_offset = aperture_diameter / 2.f;

    float center_pt[3] = {0.f, 0.f, -4.f - params.box_Z / 6.f};
    float hmax = params.box_Z;
    float hmin = center_pt[2] + cone_offset;

    gran_sys.Create_BC_Cone_Z(center_pt, cone_slope, hmax, hmin, false, false);

    ChVector<> cone_top_pos(0, 0, center_pt[2] + fill_width + 8);

    float cyl_rad = params.box_X / 2.f;

    float zvec[3] = {0, 0, 0};
    {
        std::string meshes_file = "coneflow_meshes.csv";

        std::ofstream meshfile{params.output_dir + "/" + meshes_file};
        std::ostringstream outstream;
        outstream << "mesh_name,dx,dy,dz,x1,x2,x3,y1,y2,y3,z1,z2,z3\n";
        writeZConeMesh(outstream, cone_top_pos, granular::GetDataFile("shared/gran_zcone.obj"));
        writeZCylinderMesh(outstream, ChVector<>(zvec[0], zvec[1], zvec[2]), cyl_rad, params.box_Z);

        meshfile << outstream.str();
    }

    gran_sys.Create_BC_Cyl_Z(zvec, cyl_rad, false, false);

    float plane_center[3] = {0, 0, center_pt[2] + 2 * cone_slope + cone_slope * cone_offset};
    float plane_normal[3] = {0, 0, 1};

    size_t cone_plane_bc_id = gran_sys.Create_BC_Plane(plane_center, plane_normal, false);
    gran_sys.initialize();

    float t_remove_plane = 1.f;  // Really let it settle first
    bool plane_active = false;

    int fps = 60;
    float frame_step = 1.f / fps;
    float curr_time = 0.f;
    int currframe = 0;

    // write an initial frame
    char filename[100];
    sprintf(filename, "%s/step%06d", params.output_dir.c_str(), currframe++);
    gran_sys.writeFile(std::string(filename));

    std::cout << "frame step is " << frame_step << std::endl;
    while (curr_time < params.time_end) {
        if (!plane_active && curr_time > t_remove_plane) {
            gran_sys.disable_BC_by_ID(cone_plane_bc_id);
        }

        gran_sys.advance_simulation(frame_step);
        curr_time += frame_step;
        printf("rendering frame %u\n", currframe);
        sprintf(filename, "%s/step%06d", params.output_dir.c_str(), currframe++);
        gran_sys.writeFile(std::string(filename));
    }

    return 0;
}
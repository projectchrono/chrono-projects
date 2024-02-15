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
// Chrono::Granular test program using SMC method for frictional contact for a
// Dam Break Simulation
// =============================================================================

#include <iostream>
#include <string>

#include "chrono/core/ChGlobal.h"
#include "chrono/utils/ChUtilsSamplers.h"
#include "chrono_gpu/physics/ChSystemGpu.h"
#include "chrono_gpu/utils/ChGpuJsonParser.h"
#include "chrono_thirdparty/filesystem/path.h"

#include "../utils.h"

using namespace chrono;
using namespace chrono::gpu;

enum run_mode { FRICTIONLESS_NOCYL = 0, FRICTIONLESS_WITHCYL = 1, MULTI_STEP_NOCYL = 2, MULTI_STEP_WITHCYL = 3 };

// whether or not to have a cylinder blocking the flow. Set by run_mode.
bool use_cylinder = false;

std::string box_filename = GetProjectsDataFile("gpu/meshes/BD_Box.obj");
std::string cyl_filename = GetProjectsDataFile("gpu/meshes/Gran_cylinder.obj");

ChGpuSimulationParameters params;

void writeBoxMesh(std::ostringstream& outstream) {
    ChVector3d pos(0, 0, 0);
    // Get basis vectors
    ChVector3d vx(1, 0, 0);
    ChVector3d vy(0, 1, 0);
    ChVector3d vz(0, 0, 1);

    ChVector3d scaling(params.box_X / 2, params.box_Y / 2, params.box_Z / 2);

    // Write the mesh name to find
    outstream << box_filename << ",";
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

int main(int argc, char* argv[]) {
    std::string inputJson = GetProjectsDataFile("gpu/DamBreak.json");
    int run_mode = 0;
    if (argc == 2) {
        inputJson = std::string(argv[1]);
    } else if (argc == 3) {
        inputJson = std::string(argv[1]);
        run_mode = std::atoi(argv[2]);
    } else if (argc > 1) {
        std::cout << "Usage:\n./test_GPU_DamBreak <json_file> [run_mode]" << std::endl;
        std::cout << "run_mode:  0-FRICTIONLESS_NOCYL, 1-FRICTIONLESS_WITHCYL, 2-MULTI_STEP_NOCYL, 3-MULTI_STEP_WITHCYL"
                  << std::endl;
        return 1;
    }

    ChGpuSimulationParameters params;
    if (!ParseJSON(inputJson, params)) {
        std ::cout << "ERROR: reading input file " << inputJson << std::endl;
        return 1;
    }

    std::cout << "Radius " << params.sphere_radius << std::endl;
    std::cout << "Density " << params.sphere_density << std::endl;
    std::cout << "box_Y " << params.box_Y << std::endl;
    std::cout << "output_dir " << params.output_dir << std::endl;

    // Setup simulation
    ChSystemGpu gran_sys(params.sphere_radius, params.sphere_density,
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

    gran_sys.SetCohesionRatio(0.0);
    gran_sys.SetAdhesionRatio_SPH2WALL(0.0);
    gran_sys.SetGravitationalAcceleration(ChVector3f(params.grav_X, params.grav_Y, params.grav_Z));
    gran_sys.SetParticleOutputMode(params.write_mode);

    gran_sys.SetTimeIntegrator(CHGPU_TIME_INTEGRATOR::CENTERED_DIFFERENCE);
    gran_sys.SetFixedStepSize(params.step_size);
    gran_sys.SetVerbosity(params.verbose);

    gran_sys.SetBDFixed(true);

    switch (run_mode) {
        case run_mode::MULTI_STEP_WITHCYL:
            gran_sys.SetFrictionMode(CHGPU_FRICTION_MODE::MULTI_STEP);
            use_cylinder = true;
            break;
        case run_mode::MULTI_STEP_NOCYL:
            gran_sys.SetFrictionMode(CHGPU_FRICTION_MODE::MULTI_STEP);
            use_cylinder = false;
            break;
        case run_mode::FRICTIONLESS_WITHCYL:
            gran_sys.SetFrictionMode(CHGPU_FRICTION_MODE::FRICTIONLESS);
            use_cylinder = true;
            break;
        case run_mode::FRICTIONLESS_NOCYL:
            gran_sys.SetFrictionMode(CHGPU_FRICTION_MODE::FRICTIONLESS);
            use_cylinder = false;
            break;
        default:
            std::cout << "Invalid run_mode" << std::endl;
            return 1;
    }

    // offset of radius from walls
    ChVector3f rad_offset = 1.02f * params.sphere_radius * ChVector3f(1, 1, 1);

    // (2 x 1 x 1) box (x,y,z)
    float sphere_diam = 2.f * params.sphere_radius;

    float max_z_fill = 2. * 200.;
    ChVector3f hdims = .5f * ChVector3f(2. * 100., params.box_Y, max_z_fill) - rad_offset;

    // start at bottom left corner
    ChVector3f center =
        ChVector3f(-params.box_X / 2., -params.box_Y / 2., -params.box_Z / 2.) + hdims + rad_offset;

    // Fill box with bodies
    std::vector<ChVector3f> body_points =
        utils::PDLayerSampler_BOX<float>(center, hdims, 2. * params.sphere_radius, 1.02);

    std::vector<ChVector3f> first_points;

    std::cout << "Adding " << body_points.size() << " particles" << std::endl;
    gran_sys.SetParticles(body_points);

    // just at end of material
    ChVector3f plane_center(center.x() + hdims.x() + sphere_diam, 0, 0);

    // face in -y, hold material in
    ChVector3f plane_normal(-1, 0, 0);

    printf("fill center is %f, %f, %f, plane center is %f, %f, %f\n", center.x(), center.y(), center.z(), plane_center.x(),
           plane_center.y(), plane_center.z());
    size_t plane_bc_id = gran_sys.CreateBCPlane(plane_center, plane_normal, true);

    ChVector3f cyl_center(params.box_X / 2.f - 200.f, 0, 0);

    float cyl_rad = 30;

    size_t cyl_bc_id;
    if (use_cylinder) {
        cyl_bc_id = gran_sys.CreateBCCylinderZ(cyl_center, cyl_rad, true, true);
    }

    filesystem::create_directory(filesystem::path(params.output_dir));

    // Finalize settings and initialize for runtime
    gran_sys.Initialize();

    int fps = 50;
    float frame_step = 1.0f / fps;
    float curr_time = 0;
    int currframe = 0;

    std::cout << "frame step is " << frame_step << std::endl;
    bool plane_active = true;
    ChVector3f reaction_forces(0, 0, 0);

    constexpr float F_CGS_TO_SI = 1e-5f;
    constexpr float M_CGS_TO_SI = 1e-3f;
    float total_system_mass = 4.0f / 3.0f * (float)CH_C_PI * params.sphere_density * params.sphere_radius * params.sphere_radius *
                              params.sphere_radius * body_points.size();
    printf("total system mass is %f kg \n", total_system_mass * M_CGS_TO_SI);

    std::string meshes_file = "dambreakmeshes.csv";

    // write mesh transforms for ospray renderer
    {
        std::ofstream meshfile{params.output_dir + "/" + meshes_file};
        std::ostringstream outstream;
        outstream << "mesh_name,dx,dy,dz,x1,x2,x3,y1,y2,y3,z1,z2,z3\n";
        writeBoxMesh(outstream);
        writeZCylinderMesh(outstream, cyl_center, cyl_rad, params.box_Z);

        meshfile << outstream.str();
    }

    // Run settling experiments
    while (curr_time < params.time_end) {
        if (plane_active && curr_time > 1) {
            printf("disabling plane!\n");
            plane_active = false;
            gran_sys.DisableBCbyID(plane_bc_id);
        }

        if (plane_active) {
            bool success = gran_sys.GetBCReactionForces(plane_bc_id, reaction_forces);
            if (!success) {
                printf("ERROR! Get contact forces for plane failed\n");
            } else {
                printf("curr time is %f, plane force is (%f, %f, %f) Newtons\n", curr_time,
                       F_CGS_TO_SI * reaction_forces.x(), F_CGS_TO_SI * reaction_forces.y(),
                       F_CGS_TO_SI * reaction_forces.z());
            }
        } else {
            if (use_cylinder) {
                bool success = gran_sys.GetBCReactionForces(cyl_bc_id, reaction_forces);
                if (!success) {
                    printf("ERROR! Get contact forces for cyl failed\n");
                } else {
                    printf("curr time is %f, cyl force is (%f, %f, %f) Newtons\n", curr_time,
                           F_CGS_TO_SI * reaction_forces.x(), F_CGS_TO_SI * reaction_forces.y(),
                           F_CGS_TO_SI * reaction_forces.z());
                }
            }
        }

        gran_sys.AdvanceSimulation(frame_step);
        curr_time += frame_step;
        printf("rendering frame %u\n", currframe);
        char filename[100];
        sprintf(filename, "%s/step%06d", params.output_dir.c_str(), currframe++);
        gran_sys.WriteParticleFile(std::string(filename));
    }

    return 0;
}

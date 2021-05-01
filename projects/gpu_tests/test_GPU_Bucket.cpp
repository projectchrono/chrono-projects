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
// Basic Chrono::Granular simulation of a settling granular material composed of
// spherical particles. Units are in CGS.
// =============================================================================

#include <iostream>
#include <string>
#ifdef _WINDOWS
#define _USE_MATH_DEFINES
#include <math.h>
#endif
#include "chrono/core/ChGlobal.h"
#include "chrono/utils/ChUtilsSamplers.h"

#include "chrono_gpu/physics/ChSystemGpu.h"
#include "chrono_gpu/utils/ChGpuJsonParser.h"
#include "chrono_gpu/ChGpuData.h"

#include "chrono_thirdparty/filesystem/path.h"

using namespace chrono;
using namespace chrono::gpu;

// expected number of args for param sweep
constexpr int num_args_full = 5;

int currcapture = 0;
int currframe = 0;

// -----------------------------------------------------------------------------
// Show command line usage
// -----------------------------------------------------------------------------
void ShowUsage(std::string name) {
    std::cout << "usage: " + name + " <json_file> [<run_mode> <box_Z> <output_dir>]" << std::endl;
    std::cout << "must have either 1 or " << num_args_full - 1 << " arguments" << std::endl;
}

ChGpuSimulationParameters params;

std::string box_filename = gpu::GetDataFile("meshes/BD_Box.obj");

// Take a ChBody and write its
void writeBoxMesh(std::ostringstream& outstream) {
    ChVector<> pos(0, 0, 0);
    // Get basis vectors
    ChVector<> vx(1, 0, 0);
    ChVector<> vy(0, 1, 0);
    ChVector<> vz(0, 0, 1);

    ChVector<> scaling(params.box_X / 2, params.box_Y / 2, params.box_Z / 2);

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

enum run_mode { FRICTIONLESS = 0, ONE_STEP = 1, MULTI_STEP = 2, FRICLESS_CHUNG = 3, FRICLESS_CD = 4 };

std::vector<size_t> bc_ids;
std::vector<std::string> bc_names;

constexpr float F_CGS_TO_SI = 1e-5;
constexpr float M_CGS_TO_SI = 1e-3;

void writeForcesFile(ChSystemGpu& gran_sys) {
    char forcefile[100];
    sprintf(forcefile, "%s/force%06d.csv", (params.output_dir + "/forces").c_str(), currcapture++);
    printf("force file is %s\n", forcefile);
    std::ofstream ofile(forcefile, std::ios::out);

    std::ostringstream outstrstream;

    outstrstream << "plane_id,fx,fy,fz\n";

    ChVector<float> reaction_forces;
    for (unsigned int i = 0; i < bc_ids.size(); i++) {
        bool success = gran_sys.GetBCReactionForces(bc_ids[i], reaction_forces);
        if (!success) {
            printf("ERROR! Get contact forces for plane %u failed\n", i);
        } else {
            outstrstream << i << "," << F_CGS_TO_SI * reaction_forces.x() << "," << F_CGS_TO_SI * reaction_forces.y()
                         << "," << F_CGS_TO_SI * reaction_forces.z() << "\n";
            printf("force on plane %u is (%f, %f, %f) Newtons\n", i, F_CGS_TO_SI * reaction_forces.x(),
                   F_CGS_TO_SI * reaction_forces.y(), F_CGS_TO_SI * reaction_forces.z());
        }
    }
    // delimiter to make for easier reading
    printf("--------------------------------\n");
    ofile << outstrstream.str();
}

int main(int argc, char* argv[]) {
    gpu::SetDataPath(std::string(PROJECTS_DATA_DIR) + "gpu/");

    // Some of the default values might be overwritten by user via command line
    if (argc < 2 || argc > 2 && argc != num_args_full || ParseJSON(gpu::GetDataFile(argv[1]), params) == false) {
        ShowUsage(argv[0]);
        return 1;
    }

    if (argc == num_args_full) {
        params.run_mode = std::atoi(argv[2]);
        params.box_Z = std::atof(argv[3]);
        params.output_dir = std::string(argv[4]);
        printf("new parameters: run_mode is %d, height is %f, %s\n", params.run_mode, params.box_Z,
               params.output_dir.c_str());
    }

    // Setup simulation
    ChSystemGpu gran_sys(params.sphere_radius, params.sphere_density,
                        make_float3(params.box_X, params.box_Y, params.box_Z));
    gran_sys.SetPsiFactors(params.psi_T, params.psi_L);

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
    gran_sys.SetGravitationalAcceleration(ChVector<float>(params.grav_X, params.grav_Y, params.grav_Z));
    gran_sys.SetParticleOutputMode(params.write_mode);

    std::vector<ChVector<float>> body_points;

    {
        //fill box, layer by layer
        ChVector<> hdims(params.box_X / 2.f - 3, params.box_Y / 2.f - 3, params.box_Z / 2.f - 3);
        ChVector<> center(0, 0, 0);

        // Fill box with bodies
        body_points = utils::PDLayerSampler_BOX<float>(center, hdims, 2. * params.sphere_radius, 1.02);
    }

    gran_sys.SetParticles(body_points);

    // face in upwards
    ChVector<float> bottom_plane_normal_X(1, 0, 0);
    ChVector<float> top_plane_normal_X(-1, 0, 0);
    ChVector<float> bottom_plane_normal_Y(0, 1, 0);
    ChVector<float> top_plane_normal_Y(0, -1, 0);
    ChVector<float> bottom_plane_normal_Z(0, 0, 1);
    ChVector<float> top_plane_normal_Z(0, 0, -1);

    // where to place confining planes
    float box_dist_X = params.box_X / 2.f - 2.f;
    float box_dist_Y = params.box_Y / 2.f - 2.f;
    float box_dist_Z = params.box_Z / 2.f - 2.f;

    // put a plane on the box walls
    ChVector<float> box_bottom_X(-box_dist_X, 0, 0);
    ChVector<float> box_top_X(box_dist_X, 0, 0);
    ChVector<float> box_bottom_Y(0, -box_dist_Y, 0);
    ChVector<float> box_top_Y(0, box_dist_Y, 0);
    ChVector<float> box_bottom_Z(0, 0, -box_dist_Z);
    ChVector<float> box_top_Z(0, 0, box_dist_Z);

    size_t bottom_plane_bc_id_X = gran_sys.CreateBCPlane(box_bottom_X, bottom_plane_normal_X, true);
    size_t top_plane_bc_id_X = gran_sys.CreateBCPlane(box_top_X, top_plane_normal_X, true);
    size_t bottom_plane_bc_id_Y = gran_sys.CreateBCPlane(box_bottom_Y, bottom_plane_normal_Y, true);
    size_t top_plane_bc_id_Y = gran_sys.CreateBCPlane(box_top_Y, top_plane_normal_Y, true);
    size_t bottom_plane_bc_id_Z = gran_sys.CreateBCPlane(box_bottom_Z, bottom_plane_normal_Z, true);
    size_t top_plane_bc_id_Z = gran_sys.CreateBCPlane(box_top_Z, top_plane_normal_Z, true);

    bc_ids.push_back(bottom_plane_bc_id_X);
    bc_ids.push_back(top_plane_bc_id_X);
    bc_ids.push_back(bottom_plane_bc_id_Y);
    bc_ids.push_back(top_plane_bc_id_Y);
    bc_ids.push_back(bottom_plane_bc_id_Z);
    bc_ids.push_back(top_plane_bc_id_Z);

    bc_names.push_back("bottom_plane_bc_X");
    bc_names.push_back("top_plane_bc_X");
    bc_names.push_back("bottom_plane_bc_Y");
    bc_names.push_back("top_plane_bc_Y");
    bc_names.push_back("bottom_plane_bc_Z");
    bc_names.push_back("top_plane_bc_Z");

    gran_sys.SetTimeIntegrator(CHGPU_TIME_INTEGRATOR::FORWARD_EULER);

    switch (params.run_mode) {
        case run_mode::MULTI_STEP:
            gran_sys.SetFrictionMode(CHGPU_FRICTION_MODE::MULTI_STEP);
            gran_sys.SetTimeIntegrator(CHGPU_TIME_INTEGRATOR::FORWARD_EULER);
            break;
        case run_mode::ONE_STEP:
            gran_sys.SetFrictionMode(CHGPU_FRICTION_MODE::SINGLE_STEP);
            gran_sys.SetTimeIntegrator(CHGPU_TIME_INTEGRATOR::FORWARD_EULER);
            break;
        case run_mode::FRICLESS_CHUNG:
            gran_sys.SetFrictionMode(CHGPU_FRICTION_MODE::FRICTIONLESS);
            gran_sys.SetTimeIntegrator(CHGPU_TIME_INTEGRATOR::CHUNG);
            break;

        case run_mode::FRICLESS_CD:
            gran_sys.SetFrictionMode(CHGPU_FRICTION_MODE::FRICTIONLESS);
            gran_sys.SetTimeIntegrator(CHGPU_TIME_INTEGRATOR::CENTERED_DIFFERENCE);
            break;

        case run_mode::FRICTIONLESS:
            gran_sys.SetFrictionMode(CHGPU_FRICTION_MODE::FRICTIONLESS);
            gran_sys.SetTimeIntegrator(CHGPU_TIME_INTEGRATOR::FORWARD_EULER);
            break;

        default:
            // fall through to frictionless as default
            gran_sys.SetFrictionMode(CHGPU_FRICTION_MODE::FRICTIONLESS);
    }

    gran_sys.SetFixedStepSize(params.step_size);

    filesystem::create_directory(filesystem::path(params.output_dir));
    filesystem::create_directory(filesystem::path(params.output_dir + "/forces"));

    gran_sys.SetBDFixed(true);

    gran_sys.SetVerbosity(params.verbose);
    gran_sys.Initialize();

    // number of times to capture force data per second
    int captures_per_second = 50;
    
    // number of times to capture force before we capture a frame
    int captures_per_frame = 2;

    // assume we run for at least one frame
    float frame_step = 1. / captures_per_second;
    float curr_time = 0;

    std::cout << "capture step is " << frame_step << std::endl;

    float total_system_mass = 4. / 3. * CH_C_PI * params.sphere_density * params.sphere_radius * params.sphere_radius *
                              params.sphere_radius * body_points.size();
    printf("total system mass is %f kg \n", total_system_mass * M_CGS_TO_SI);

    // write an initial frame
    char filename[100];
    printf("rendering frame %u\n", currframe);
    sprintf(filename, "%s/step%06d", params.output_dir.c_str(), currframe);
    gran_sys.WriteParticleFile(std::string(filename));

    // write mesh transforms for ospray renderer
    {
        std::ofstream meshfile{params.output_dir + "/BD_Box_mesh.csv"};
        std::ostringstream outstream;
        outstream << "mesh_name,dx,dy,dz,x1,x2,x3,y1,y2,y3,z1,z2,z3\n";
        writeBoxMesh(outstream);
        meshfile << outstream.str();
    }
    std::cout << "frame step is " << frame_step << std::endl;

    currframe++;

    // Run settling experiments
    while (curr_time < params.time_end) {
        printf("curr time is %f\n", curr_time);

        writeForcesFile(gran_sys);
        gran_sys.AdvanceSimulation(frame_step);
        curr_time += frame_step;

        // if this frame is a render frame
        if (currcapture % captures_per_frame == 0) {
            // printf("rendering frame %u\n", currframe);
            sprintf(filename, "%s/step%06d", params.output_dir.c_str(), currframe++);
            // gran_sys.WriteFile(std::string(filename));
        }



    }

    return 0;
}

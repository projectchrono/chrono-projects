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
// Authors: Nic Olsen
// =============================================================================
// Chrono::Granular simulation in which a cylinder is filled with granular
// material and then raised slightly, allowing material to flow out and around
// the cylinder for comparison with the analytical hydrostatic result.
// =============================================================================

#include <cmath>
#include <iostream>
#include <string>

#include "chrono/core/ChGlobal.h"
#include "chrono/utils/ChUtilsSamplers.h"
#include "chrono_gpu/physics/ChSystemGpu.h"
#include "chrono_gpu/utils/ChGpuJsonParser.h"
#include "chrono_gpu/utils/ChGpuVisualization.h"
#include "chrono_thirdparty/filesystem/path.h"

#include "../utils.h"

using namespace chrono;
using namespace chrono::gpu;

void writeMeshFrames(std::ostringstream& outstream,
                     const std::string obj_name,
                     const ChVector3d& pos,
                     const ChVector3d& mesh_scaling) {
    outstream << obj_name << ",";

    // Get basis vectors
    ChVector3d vx(1, 0, 0);
    ChVector3d vy(0, 1, 0);
    ChVector3d vz(0, 0, 1);

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
    outstream << mesh_scaling.x() << "," << mesh_scaling.y() << "," << mesh_scaling.z();
    outstream << "\n";
}

int main(int argc, char* argv[]) {
    std::string inputJson = GetProjectsDataFile("gpu/comvessels.json");
    if (argc == 2) {
        inputJson = std::string(argv[1]);
    } else if (argc > 1) {
        std::cout << "Usage:\n./test_GPU_commvessels <json_file>" << std::endl;
        return 1;
    }

    ChGpuSimulationParameters params;
    if (!ParseJSON(inputJson, params)) {
        std ::cout << "ERROR: reading input file " << inputJson << std::endl;
        return 1;
    }

    float iteration_step = params.step_size;
    const float Bx = params.box_X;
    const float By = params.box_Y;
    const float Bz = params.box_Z;

    // Overwrite parameters from the command line
    std::cout << "sphere_radius " << params.sphere_radius << std::endl;
    std::cout << "sphere_density " << params.sphere_density << std::endl;

    ChSystemGpuMesh gran_sys(params.sphere_radius, params.sphere_density, ChVector3f(Bx, By, Bz));

    gran_sys.SetVerbosity(params.verbose);

    gran_sys.SetKn_SPH2SPH(params.normalStiffS2S);
    gran_sys.SetKn_SPH2WALL(params.normalStiffS2W);
    gran_sys.SetKn_SPH2MESH(params.normalStiffS2M);

    gran_sys.SetKt_SPH2SPH(params.tangentStiffS2S);
    gran_sys.SetKt_SPH2WALL(params.tangentStiffS2W);
    gran_sys.SetKt_SPH2MESH(params.tangentStiffS2M);

    gran_sys.SetGn_SPH2SPH(params.normalDampS2S);
    gran_sys.SetGn_SPH2WALL(params.normalDampS2W);
    gran_sys.SetGn_SPH2MESH(params.normalDampS2M);

    gran_sys.SetGt_SPH2SPH(params.tangentDampS2S);
    gran_sys.SetGt_SPH2WALL(params.tangentDampS2W);
    gran_sys.SetGt_SPH2MESH(params.tangentDampS2M);

    gran_sys.SetCohesionRatio(params.cohesion_ratio);
    gran_sys.SetAdhesionRatio_SPH2MESH(params.adhesion_ratio_s2m);
    gran_sys.SetAdhesionRatio_SPH2WALL(params.adhesion_ratio_s2w);
    gran_sys.SetFrictionMode(chrono::gpu::CHGPU_FRICTION_MODE::FRICTIONLESS);

    gran_sys.SetStaticFrictionCoeff_SPH2SPH(params.static_friction_coeffS2S);
    gran_sys.SetStaticFrictionCoeff_SPH2WALL(params.static_friction_coeffS2W);
    gran_sys.SetStaticFrictionCoeff_SPH2MESH(params.static_friction_coeffS2M);

    gran_sys.SetParticleOutputMode(params.write_mode);
    filesystem::create_directory(filesystem::path(params.output_dir));
    gran_sys.SetTimeIntegrator(CHGPU_TIME_INTEGRATOR::CENTERED_DIFFERENCE);

    gran_sys.SetFixedStepSize(params.step_size);

    gran_sys.SetBDFixed(true);
    gran_sys.SetGravitationalAcceleration(ChVector3f(params.grav_X, params.grav_Y, params.grav_Z));

    // Fill the entire height
    const float fill_bottom = -Bz / 2.f;
    const float fill_height = Bz;

    // Cylinder mesh has interior radius 1 and total radius 1.1
    ChVector3f cyl_center(0, 0, 0);
    ChVector3f scaling(Bx / 4.f, Bx / 4.f, Bz);
    std::cout << "Cylinder radius: " << scaling.x() << std::endl;

    utils::ChPDSampler<float> sampler(2.05 * params.sphere_radius);
    std::vector<ChVector3f> body_points;

    const float fill_radius = scaling.x() - 2.f * params.sphere_radius;
    const float fill_top = fill_bottom + fill_height;
    std::cout << "Fill radius " << fill_radius << std::endl;
    std::cout << "Fill bottom " << fill_bottom << std::endl;
    std::cout << "Fill top " << fill_top << std::endl;

    ChVector3f center(0, 0, fill_bottom);
    center.z() += 2 * params.sphere_radius;
    while (center.z() < fill_top - 2 * params.sphere_radius) {
        auto points = sampler.SampleCylinderZ(center, fill_radius, 0);
        body_points.insert(body_points.end(), points.begin(), points.end());
        center.z() += 2.05 * params.sphere_radius;
    }

    auto n_spheres = body_points.size();
    std::cout << "Adding " << n_spheres << " particles" << std::endl;
    gran_sys.SetParticles(body_points);

    // Add mesh
    std::string mesh_filename(GetProjectsDataFile("gpu/meshes/cylinder_refined.obj"));
    gran_sys.AddMesh(mesh_filename, ChVector3f(0), ChMatrix33<float>(scaling), 10.0f);

    gran_sys.Initialize();
    const double time_settling = std::sqrt(-2.0 * (params.box_Z) / params.grav_Z);
    const double raising_vel = 1.0;

    const double raising_dist = 10 * 2.0 * 0.2;  // Hard-coded to be the same height as the 0.2 radius run
    const double time_raising = raising_dist / raising_vel;
    const double time_sitting = 10.0;  // TODO no idea how much is enough

    std::cout << "Time settling " << time_settling << std::endl;
    std::cout << "Raising velocity " << raising_vel << std::endl;
    std::cout << "Raising distance " << raising_dist << std::endl;
    std::cout << "Time raising " << time_raising << std::endl;
    std::cout << "Time sitting " << time_sitting << std::endl;

    double mesh_z = 0.0;
    double mesh_vz = raising_vel;

    // Set initial mesh locations for the settling phase
    ChVector3f mesh_pos(0, 0, 0);
    ChQuaternion<float> mesh_rot(1, 0, 0, 0);
    ChVector3f mesh_lin_vel(0, 0, 0);
    ChVector3f mesh_ang_vel(0, 0, 0);

    gran_sys.ApplyMeshMotion(0, mesh_pos, mesh_rot, mesh_lin_vel, mesh_ang_vel);

    unsigned int currframe = 0;
    float out_fps = 60;
    float frame_step = 1.f / out_fps;  // Duration of a frame
    unsigned int out_steps = (unsigned int)(frame_step / iteration_step);
    std::cout << "Writing at " << out_fps << " FPS" << std::endl;

    unsigned int step = 0;
    bool settled = false;
    bool raised = false;
    //    ChVector3d pos_mesh(0, 0, mesh_z);

    std::cout << "Settling..." << std::endl;
    for (float t = 0; t < time_settling + time_raising + time_sitting; t += iteration_step, step++) {
        if (t >= time_settling && t <= time_settling + time_raising) {
            // Raising phase
            if (!settled) {
                std::cout << "Raising..." << std::endl;
                settled = true;
            }

            mesh_z += iteration_step * raising_vel;
            mesh_pos.z() = mesh_z;
            mesh_lin_vel.z() = mesh_vz;
            gran_sys.ApplyMeshMotion(0, mesh_pos, mesh_rot, mesh_lin_vel, mesh_ang_vel);
        } else if (t > time_settling + time_raising) {
            if (!raised) {
                std::cout << "Raised." << std::endl;
                raised = true;
                mesh_lin_vel.z() = 0;
                gran_sys.ApplyMeshMotion(0, mesh_pos, mesh_rot, mesh_lin_vel, mesh_ang_vel);
            }
        }

        if (step % out_steps == 0) {
            std::cout << "Rendering frame " << currframe << std::endl;
            char filename[100];
            sprintf(filename, "%s/step%06u", params.output_dir.c_str(), currframe++);
            gran_sys.WriteParticleFile(std::string(filename));
            // gran_sys.write_meshes(std::string(filename));
            std::string mesh_output = std::string(filename) + "_meshframes.csv";

            std::ofstream meshfile(mesh_output);
            std::ostringstream outstream;
            outstream << "mesh_name,dx,dy,dz,x1,x2,x3,y1,y2,y3,z1,z2,z3,sx,sy,sz\n";

            writeMeshFrames(outstream, mesh_filename, mesh_pos, scaling);

            meshfile << outstream.str();
            meshfile.close();
        }

        gran_sys.AdvanceSimulation(iteration_step);
    }

    return 0;
}

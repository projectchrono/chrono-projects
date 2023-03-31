// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2021 projectchrono.org
// All right reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Ruochun Zhang
// =============================================================================
//
// A meshed ball hitting a granular bed under gravity.
//
// =============================================================================

#include "chrono/assets/ChSphereShape.h"
#include "chrono/core/ChGlobal.h"
#include "chrono/physics/ChBody.h"
#include "chrono/physics/ChForce.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/timestepper/ChTimestepper.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsSamplers.h"

#include <DEM/API.h>
#include <core/ApiVersion.h>
#include <core/utils/ThreadManager.h>
#include <DEM/HostSideHelpers.hpp>
#include <DEM/utils/Samplers.hpp>

#include <chrono>
#include <cstdio>
#include <filesystem>

using namespace deme;
using namespace chrono;
using namespace std::filesystem;

// ChVector to/from float3
inline float3 ChVec2Float(const ChVector<>& vec) {
    return make_float3(vec.x(), vec.y(), vec.z());
}
inline ChVector<> Float2ChVec(float3 f3) {
    return ChVector<>(f3.x, f3.y, f3.z);
}

inline float4 ChQ2Float(const ChQuaternion<>& Q) {
    float4 f4;
    f4.w = Q.e0();
    f4.x = Q.e1();
    f4.y = Q.e2();
    f4.z = Q.e3();
    return f4;
}

int main() {
    // World info
    float step_size = 1e-5;
    double world_size = 10;

    // Set Chromo system to manage a solid ball
    float ball_radius = 1.f;
    float ball_density = 7.8e3;
    float ball_mass = ball_density * 4 / 3 * 3.1416;
    // Create rigid ball_body simulation
    ChSystemSMC sys_ball;
    sys_ball.SetContactForceModel(ChSystemSMC::ContactForceModel::Hooke);
    sys_ball.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT);
    sys_ball.Set_G_acc(ChVector<>(0, 0, -9.81));

    double inertia = ball_mass * 2 / 5;
    ChVector<> ball_initial_pos(world_size / 2., world_size / 2., world_size / 3. * 2.);

    std::shared_ptr<ChBody> ball_body(sys_ball.NewBody());
    ball_body->SetMass(ball_mass);
    ball_body->SetInertiaXX(ChVector<>(inertia, inertia, inertia));
    ball_body->SetPos(ball_initial_pos);
    auto sph = chrono_types::make_shared<ChSphereShape>();
    sph->GetSphereGeometry().rad = ball_radius;
    ball_body->AddVisualShape(sph);
    sys_ball.AddBody(ball_body);

    // Then set up DEME system
    DEMSolver DEMSim;
    DEMSim.SetVerbosity(INFO);
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetOutputContent(OUTPUT_CONTENT::ABSV);
    DEMSim.SetMeshOutputFormat(MESH_FORMAT::VTK);

    // E, nu, CoR, mu, Crr...
    auto mat_type_ball = DEMSim.LoadMaterial({{"E", 1e10}, {"nu", 0.3}, {"CoR", 0.6}, {"mu", 0.3}, {"Crr", 0.01}});
    auto mat_type_terrain = DEMSim.LoadMaterial({{"E", 5e9}, {"nu", 0.3}, {"CoR", 0.8}, {"mu", 0.3}, {"Crr", 0.01}});
    // If you don't have this line, then CoR between mixer material and granular material will be 0.7 (average of the
    // two).
    DEMSim.SetMaterialPropertyPair("CoR", mat_type_ball, mat_type_terrain, 0.6);
    // Should do the same for mu and Crr, but since they are the same across 2 materials, it won't have an effect...

    DEMSim.InstructBoxDomainDimension({0, world_size}, {0, world_size}, {0, world_size});
    DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_terrain);

    auto projectile = DEMSim.AddWavefrontMeshObject((GET_DATA_PATH() / "mesh/sphere.obj").string(), mat_type_ball);
    std::cout << "Total num of triangles: " << projectile->GetNumTriangles() << std::endl;

    projectile->SetInitPos(ChVec2Float(ball_initial_pos));
    projectile->SetMass(ball_mass);
    projectile->SetMOI(make_float3(inertia, inertia, inertia));
    // Defaulted to family 255 which is fixed
    // Track the projectile
    auto proj_tracker = DEMSim.Track(projectile);

    float terrain_rad = 0.05;
    auto template_terrain = DEMSim.LoadSphereType(terrain_rad * terrain_rad * terrain_rad * 2.6e3 * 4 / 3 * 3.14,
                                                  terrain_rad, mat_type_terrain);
    float sample_halfheight = world_size / 8;
    float3 sample_center = make_float3(world_size / 2, world_size / 2, sample_halfheight + 0.05);
    float sample_halfwidth = world_size / 2 * 0.95;
    auto input_xyz = DEMBoxHCPSampler(sample_center, make_float3(sample_halfwidth, sample_halfwidth, sample_halfheight),
                                      2.01 * terrain_rad);
    DEMSim.AddClumps(template_terrain, input_xyz);
    std::cout << "Total num of particles: " << input_xyz.size() << std::endl;

    DEMSim.SetInitTimeStep(step_size);
    DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -9.81));
    // Max velocity info is generally just for the solver's reference and the user do not have to set it. The solver
    // wouldn't take into account a vel larger than this when doing async-ed contact detection: but this vel won't
    // happen anyway and if it does, something already went wrong.
    DEMSim.SetMaxVelocity(15.);
    // In general you don't have to worry about SetExpandSafetyAdder, unless if an entity has the property that a point
    // on it can move much faster than its CoM. In this demo, you are dealing with a meshed ball and you in fact don't
    // have this problem. In the Centrifuge demo though, this can be a problem since the centrifuge's CoM is not moving,
    // but its pointwise velocity can be high, so it needs to be accounted for using this method.
    DEMSim.SetExpandSafetyAdder(5.);
    DEMSim.SetInitBinSize(4 * terrain_rad);
    DEMSim.Initialize();

    path out_dir = GetChronoOutputPath();
    create_directory(out_dir);
    out_dir += "/DEME";
    create_directory(out_dir);
    out_dir += "/BallDrop_Cosim";
    create_directory(out_dir);

    float sim_time = 6.0;
    float settle_time = 1.0;
    unsigned int fps = 10;
    float frame_time = 1.0 / fps;

    std::cout << "Output at " << fps << " FPS" << std::endl;
    unsigned int currframe = 0;

    // We can let it settle first
    for (float t = 0; t < settle_time; t += frame_time) {
        std::cout << "Frame: " << currframe << std::endl;
        char filename[200], meshfilename[200];
        sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe);
        sprintf(meshfilename, "%s/DEMdemo_mesh_%04d.vtk", out_dir.c_str(), currframe);
        DEMSim.WriteSphereFile(std::string(filename));
        DEMSim.WriteMeshFile(std::string(meshfilename));
        currframe++;

        DEMSim.DoDynamicsThenSync(frame_time);
    }
    DEMSim.ShowThreadCollaborationStats();

    // Then drop the ball, literally.
    DEMSim.ChangeFamily(2, 1);
    unsigned int out_steps = (unsigned int)(1.0 / (fps * step_size));
    unsigned int curr_step = 0;
    for (float t = 0; t < sim_time; t += step_size, curr_step++) {
        if (curr_step % out_steps == 0) {
            std::cout << "Frame: " << currframe << std::endl;
            char filename[200], meshfilename[200];
            sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe);
            sprintf(meshfilename, "%s/DEMdemo_mesh_%04d.vtk", out_dir.c_str(), currframe);
            DEMSim.WriteSphereFile(std::string(filename));
            DEMSim.WriteMeshFile(std::string(meshfilename));
            currframe++;
        }
        proj_tracker->SetPos(ChVec2Float(ball_body->GetPos()));
        proj_tracker->SetOriQ(ChQ2Float(ball_body->GetRot()));
        proj_tracker->SetVel(ChVec2Float(ball_body->GetPos_dt()));
        proj_tracker->SetAngVel(ChVec2Float(ball_body->GetWvel_par()));

        {
            ChVector<> ball_force;
            ChVector<> ball_torque;
            float3 F = proj_tracker->ContactAcc();
            F *= ball_mass;
            float3 tor = proj_tracker->ContactAngAccLocal();
            tor = make_float3(inertia) * tor;
            ball_force = Float2ChVec(F);
            ball_torque = Float2ChVec(tor);

            ball_body->Empty_forces_accumulators();
            ball_body->Accumulate_force(ball_force, ball_body->GetPos(), false);
            ball_body->Accumulate_torque(ball_torque, true);  // torque in DEME is local
        }

        sys_ball.DoStepDynamics(step_size);
        DEMSim.DoStepDynamics();
    }

    DEMSim.ShowThreadCollaborationStats();
    DEMSim.ShowAnomalies();
    std::cout << "DEMdemo_BallDrop exiting..." << std::endl;
    return 0;
}

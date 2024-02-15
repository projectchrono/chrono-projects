// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All right reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Author: Radu Serban
// =============================================================================
//
// Chrono::Multicore test program for active bbox.
//
// The global reference frame has Z up.
// All units SI.
// =============================================================================

#include <cstdio>
#include <cmath>

#include "chrono/ChConfig.h"
#include "chrono/assets/ChVisualShapeBox.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono_multicore/physics/ChSystemMulticore.h"
#include "chrono_opengl/ChVisualSystemOpenGL.h"

using namespace chrono;

// =============================================================================
int main(int argc, char* argv[]) {
    bool aabb_active = true;
    ChContactMethod method = ChContactMethod::SMC;

    // Create system and set method-specific solver settings
    ChSystemMulticore* system;
    double time_step;
    switch (method) {
        case ChContactMethod::SMC: {
            ChSystemMulticoreSMC* sys = new ChSystemMulticoreSMC;
            sys->GetSettings()->solver.contact_force_model = ChSystemSMC::Hertz;
            sys->GetSettings()->solver.tangential_displ_mode = ChSystemSMC::TangentialDisplacementModel::OneStep;
            sys->GetSettings()->solver.adhesion_force_model = ChSystemSMC::AdhesionForceModel::Constant;
            sys->GetSettings()->solver.use_material_properties = true;
            system = sys;
            time_step = 1e-3;
            break;
        }
        case ChContactMethod::NSC: {
            ChSystemMulticoreNSC* sys = new ChSystemMulticoreNSC;
            sys->GetSettings()->solver.solver_type = SolverType::BB;
            sys->GetSettings()->solver.solver_mode = SolverMode::SLIDING;
            sys->GetSettings()->solver.max_iteration_normal = 0;
            sys->GetSettings()->solver.max_iteration_sliding = 200;
            sys->GetSettings()->solver.max_iteration_spinning = 0;
            sys->GetSettings()->solver.alpha = 0;
            sys->GetSettings()->solver.contact_recovery_speed = -1;
            system = sys;
            time_step = 1e-3;
            break;
        }
    }

    system->SetCollisionSystemType(ChCollisionSystem::Type::MULTICORE);

    // Set method-independent solver settings
    system->Set_G_acc(ChVector3d(0, 0, -9.8));
    system->GetSettings()->solver.use_full_inertia_tensor = false;
    system->GetSettings()->solver.tolerance = 1e-3;
    system->GetSettings()->collision.bins_per_axis = vec3(10, 10, 10);
    system->GetSettings()->collision.narrowphase_algorithm = ChNarrowphase::Algorithm::HYBRID;

    // Set number of threads.
    int threads = 4;
    int max_threads = omp_get_num_procs();
    if (threads > max_threads)
        threads = max_threads;
    omp_set_num_threads(threads);

    // Create ground body
    std::shared_ptr<chrono::ChContactMaterial> material_g;
    switch (method) {
        case ChContactMethod::SMC: {
            auto mat_g = chrono_types::make_shared<ChContactMaterialSMC>();
            mat_g->SetYoungModulus(1e7f);
            mat_g->SetFriction(0.7f);
            mat_g->SetRestitution(0.5f);
            mat_g->SetAdhesion(0.0f);
            material_g = mat_g;
            break;
        }
        case ChContactMethod::NSC: {
            auto mat_g = chrono_types::make_shared<ChContactMaterialNSC>();
            mat_g->SetFriction(0.7f);
            mat_g->SetRestitution(0.5f);
            mat_g->SetCohesion(0.0f);
            material_g = mat_g;
            break;
        }
    }

    auto ground = chrono_types::make_shared<ChBody>();
    ground->SetPos(ChVector3d(0, 0, 0));
    ground->SetBodyFixed(true);
    ground->SetCollide(true);
    utils::AddBoxGeometry(ground.get(), material_g, ChVector3d(1, 1, 0.1), ChVector3d(0, 0, -0.1));
    system->AddBody(ground);

    // Create ball body
    std::shared_ptr<chrono::ChContactMaterial> material_b;
    switch (method) {
        case ChContactMethod::SMC: {
            auto mat_b = chrono_types::make_shared<ChContactMaterialSMC>();
            mat_b->SetYoungModulus(1e7f);
            mat_b->SetFriction(0.7f);
            mat_b->SetRestitution(0.5f);
            mat_b->SetAdhesion(0.0f);
            material_b = mat_b;
            break;
        }
        case ChContactMethod::NSC: {
            auto mat_b = chrono_types::make_shared<ChContactMaterialNSC>();
            mat_b->SetFriction(0.7f);
            mat_b->SetRestitution(0.5f);
            mat_b->SetCohesion(0.0f);
            material_b = mat_b;
            break;
        }
    }

    auto ball = chrono_types::make_shared<ChBody>();
    ball->SetPos(ChVector3d(0, 0, 1));
    ball->SetPos_dt(ChVector3d(0, 0, 8));
    ball->SetCollide(true);
    utils::AddSphereGeometry(ball.get(), material_b, 0.2);
    system->AddBody(ball);

    // Enable deactivation of bodies that exit a specified bounding box
    system->GetSettings()->collision.use_aabb_active = aabb_active;
    real3 bmin(-1.1, -1.1, -0.5);
    real3 bmax(+1.1, +1.1, +3.0);
    system->GetSettings()->collision.aabb_min = bmin;
    system->GetSettings()->collision.aabb_max = bmax;

    real3 bsize = (bmax - bmin);
    real3 bpos = 0.5 * (bmax + bmin);
    auto bbox = chrono_types::make_shared<ChVisualShapeBox>(bsize.x, bsize.y, bsize.z);
    ground->AddVisualShape(bbox, ChFrame<>(ChVector3d(bpos.x, bpos.y, bpos.z)));

    std::cout << "Bmin:  " << bmin.x << " " << bmin.y << " " << bmin.z << std::endl;
    std::cout << "Bmax:  " << bmax.x << " " << bmax.y << " " << bmax.z << std::endl;
    std::cout << "Bsize: " << bsize.x << " " << bsize.y << " " << bsize.z << std::endl;
    std::cout << "Bpos:  " << bpos.x << " " << bpos.y << " " << bpos.z << std::endl;

    // Initialize OpenGL
    opengl::ChVisualSystemOpenGL vis;
    vis.AttachSystem(system);
    vis.SetWindowTitle("Test AABB");
    vis.SetWindowSize(1280, 720);
    vis.SetRenderMode(opengl::WIREFRAME);
    vis.Initialize();
    vis.AddCamera(ChVector3d(0, -10, 0), ChVector3d(0, 0, 0));
    vis.SetCameraVertical(CameraVerticalDir::Z);

    // Run simulation loop
    while (vis.Run()) {
            system->DoStepDynamics(time_step);
            vis.Render();
    }

    return 0;
}

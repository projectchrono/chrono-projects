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
// ChronoParallel test program for active bbox.
//
// The global reference frame has Z up.
// All units SI.
// =============================================================================

#include <cstdio>
#include <cmath>

#include "chrono/ChConfig.h"
#include "chrono/assets/ChBoxShape.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_opengl/ChOpenGLWindow.h"

using namespace chrono;
using namespace chrono::collision;

// =============================================================================
int main(int argc, char* argv[]) {
    bool aabb_active = true;
    ChContactMethod method = ChContactMethod::SMC;

    // Create system and set method-specific solver settings
    ChSystemParallel* system;
    double time_step;
    switch (method) {
        case ChContactMethod::SMC: {
            ChSystemParallelSMC* sys = new ChSystemParallelSMC;
            sys->GetSettings()->solver.contact_force_model = ChSystemSMC::Hertz;
            sys->GetSettings()->solver.tangential_displ_mode = ChSystemSMC::TangentialDisplacementModel::OneStep;
            sys->GetSettings()->solver.adhesion_force_model = ChSystemSMC::AdhesionForceModel::Constant;
            sys->GetSettings()->solver.use_material_properties = true;
            system = sys;
            time_step = 1e-3;
            break;
        }
        case ChContactMethod::NSC: {
            ChSystemParallelNSC* sys = new ChSystemParallelNSC;
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

    // Set method-independent solver settings
    system->Set_G_acc(ChVector<>(0, 0, -9.8));
    system->GetSettings()->solver.use_full_inertia_tensor = false;
    system->GetSettings()->solver.tolerance = 1e-3;
    system->GetSettings()->collision.bins_per_axis = vec3(10, 10, 10);
    system->GetSettings()->collision.narrowphase_algorithm = NarrowPhaseType::NARROWPHASE_HYBRID_MPR;

    // Set number of threads.
    int threads = 4;
    int max_threads = omp_get_num_procs();
    if (threads > max_threads)
        threads = max_threads;
    omp_set_num_threads(threads);

    // Create ground body
    std::shared_ptr<chrono::ChMaterialSurface> material_g;
    switch (method) {
        case ChContactMethod::SMC: {
            auto mat_g = chrono_types::make_shared<ChMaterialSurfaceSMC>();
            mat_g->SetYoungModulus(1e7f);
            mat_g->SetFriction(0.7f);
            mat_g->SetRestitution(0.5f);
            mat_g->SetAdhesion(0.0f);
            material_g = mat_g;
            break;
        }
        case ChContactMethod::NSC: {
            auto mat_g = chrono_types::make_shared<ChMaterialSurfaceNSC>();
            mat_g->SetFriction(0.7f);
            mat_g->SetRestitution(0.5f);
            mat_g->SetCohesion(0.0f);
            material_g = mat_g;
            break;
        }
    }

    auto ground = std::shared_ptr<ChBody>(system->NewBody());
    ground->SetPos(ChVector<>(0, 0, 0));
    ground->SetBodyFixed(true);
    ground->SetCollide(true);
    ground->GetCollisionModel()->ClearModel();
    utils::AddBoxGeometry(ground.get(), material_g, ChVector<>(1, 1, 0.1), ChVector<>(0, 0, -0.1));
    ground->GetCollisionModel()->BuildModel();
    system->AddBody(ground);

    // Create ball body
    std::shared_ptr<chrono::ChMaterialSurface> material_b;
    switch (method) {
        case ChContactMethod::SMC: {
            auto mat_b = chrono_types::make_shared<ChMaterialSurfaceSMC>();
            mat_b->SetYoungModulus(1e7f);
            mat_b->SetFriction(0.7f);
            mat_b->SetRestitution(0.5f);
            mat_b->SetAdhesion(0.0f);
            material_b = mat_b;
            break;
        }
        case ChContactMethod::NSC: {
            auto mat_b = chrono_types::make_shared<ChMaterialSurfaceNSC>();
            mat_b->SetFriction(0.7f);
            mat_b->SetRestitution(0.5f);
            mat_b->SetCohesion(0.0f);
            material_b = mat_b;
            break;
        }
    }

    auto ball = std::shared_ptr<ChBody>(system->NewBody());
    ball->SetPos(ChVector<>(0, 0, 1));
    ball->SetPos_dt(ChVector<>(0, 0, 8));
    ball->SetCollide(true);
    ball->GetCollisionModel()->ClearModel();
    utils::AddSphereGeometry(ball.get(), material_b, 0.2);
    ball->GetCollisionModel()->BuildModel();
    system->AddBody(ball);

    // Enable deactivation of bodies that exit a specified bounding box
    system->GetSettings()->collision.use_aabb_active = aabb_active;
    real3 bmin(-1.1, -1.1, -0.5);
    real3 bmax(+1.1, +1.1, +3.0);
    system->GetSettings()->collision.aabb_min = bmin;
    system->GetSettings()->collision.aabb_max = bmax;

    auto bbox = chrono_types::make_shared<ChBoxShape>();
    real3 bsize = 0.5 * (bmax - bmin);
    real3 bpos = 0.5 * (bmax + bmin);
    bbox->GetBoxGeometry().Size = ChVector<>(bsize.x, bsize.y, bsize.z);
    bbox->Pos = ChVector<>(bpos.x, bpos.y, bpos.z);
    ground->AddAsset(bbox);

    std::cout << "Bmin:  " << bmin.x << " " << bmin.y << " " << bmin.z << std::endl;
    std::cout << "Bmax:  " << bmax.x << " " << bmax.y << " " << bmax.z << std::endl;
    std::cout << "Bsize: " << bsize.x << " " << bsize.y << " " << bsize.z << std::endl;
    std::cout << "Bpos:  " << bpos.x << " " << bpos.y << " " << bpos.z << std::endl;

    // Initialize OpenGL
    opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
    gl_window.Initialize(1280, 720, "Test AABB", system);
    gl_window.SetCamera(ChVector<>(0, -10, 0), ChVector<>(0, 0, 0), ChVector<>(0, 0, 1));
    gl_window.SetRenderMode(opengl::WIREFRAME);

    // Run simulation loop
    while (true) {
        if (gl_window.Active()) {
            gl_window.DoStepDynamics(time_step);
            gl_window.Render();
        } else {
            break;
        }
    }

    return 0;
}

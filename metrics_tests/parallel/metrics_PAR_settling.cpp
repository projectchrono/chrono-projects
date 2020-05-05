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
// ChronoParallel test program for settling process of granular material.
//
// The global reference frame has Z up.
// All units SI (CGS, i.e., centimeter - gram - second)
//
// =============================================================================

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <valarray>
#include <cstdlib>
#include <vector>

#include "chrono/ChConfig.h"
#include "chrono/core/ChTimer.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsGenerators.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_parallel/solver/ChIterativeSolverParallel.h"

#include "chrono_thirdparty/filesystem/path.h"

#ifdef CHRONO_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif

#include "../BaseTest.h"

using namespace chrono;

// --------------------------------------------------------------------------

void TimingHeader() {
    printf("    TIME    |");
    printf("    STEP |");
    printf("   BROAD |");
    printf("  NARROW |");
    printf("  SOLVER |");
    printf("  UPDATE |");
    printf("# BODIES |");
    printf("# CONTACT|");
    printf(" # ITERS |");
    printf("   RESID |");
    printf("\n\n");
}

void TimingOutput(chrono::ChSystem* mSys) {
    double TIME = mSys->GetChTime();
    double STEP = mSys->GetTimerStep();
    double BROD = mSys->GetTimerCollisionBroad();
    double NARR = mSys->GetTimerCollisionNarrow();
    double SOLVER = mSys->GetTimerAdvance();
    double UPDT = mSys->GetTimerUpdate();
    double RESID = 0;
    int REQ_ITS = 0;
    int BODS = mSys->GetNbodies();
    int CNTC = mSys->GetNcontacts();
    if (chrono::ChSystemParallel* parallel_sys = dynamic_cast<chrono::ChSystemParallel*>(mSys)) {
        RESID = std::static_pointer_cast<chrono::ChIterativeSolverParallel>(mSys->GetSolver())->GetResidual();
        REQ_ITS = std::static_pointer_cast<chrono::ChIterativeSolverParallel>(mSys->GetSolver())->GetIterations();
        BODS = parallel_sys->GetNbodies();
        CNTC = parallel_sys->GetNcontacts();
    }

    printf("   %8.5f | %7.4f | %7.4f | %7.4f | %7.4f | %7.4f | %7d | %7d | %7d | %7.4f |\n", TIME, STEP, BROD, NARR,
           SOLVER, UPDT, BODS, CNTC, REQ_ITS, RESID);
}

// ====================================================================================

// Test class
class PARSettlingTest : public BaseTest {
  public:
    PARSettlingTest(const std::string& testName,
                    const std::string& testProjectName,
                    ChContactMethod method,
                    int num_threads)
        : BaseTest(testName, testProjectName), m_method(method), m_execTime(0), m_num_threads(num_threads) {}

    ~PARSettlingTest() {}

    // Override corresponding functions in BaseTest
    virtual bool execute() override;
    virtual double getExecutionTime() const override { return m_execTime; }

  private:
    ChContactMethod m_method;
    double m_execTime;
    int m_num_threads;
};

// ====================================================================================

bool PARSettlingTest::execute() {
    bool use_mat_properties = true;
    bool render = false;

    std::cout << "Test: " << getTestName() << std::endl;
    std::cout << "Requested number of threads: " << m_num_threads << std::endl;

    // ----------------
    // Model parameters
    // ----------------

    // Container dimensions
    double hdimX = 2.0;
    double hdimY = 0.25;
    double hdimZ = 0.5;
    double hthick = 0.25;

    // Granular material properties
    double radius_g = 0.05;
    int Id_g = 10000;
    double rho_g = 2500;
    double vol_g = (4.0 / 3) * CH_C_PI * radius_g * radius_g * radius_g;
    double mass_g = rho_g * vol_g;
    ChVector<> inertia_g = 0.4 * mass_g * radius_g * radius_g * ChVector<>(1, 1, 1);
    int num_layers = 10;

    // Terrain contact properties
    float friction_terrain = 0.9f;
    float restitution_terrain = 0.0f;
    float Y_terrain = 8e5f;
    float nu_terrain = 0.3f;
    float kn_terrain = 1.0e7f;
    float gn_terrain = 1.0e3f;
    float kt_terrain = 2.86e6f;
    float gt_terrain = 1.0e3f;

    // Estimates for number of bins for broad-phase
    int factor = 2;
    int binsX = (int)std::ceil(hdimX / radius_g) / factor;
    int binsY = (int)std::ceil(hdimY / radius_g) / factor;
    int binsZ = 1;
    std::cout << "Broad-phase bins: " << binsX << " x " << binsY << " x " << binsZ << std::endl;

    // --------------------------
    // Create the parallel system
    // --------------------------

    // Create system and set method-specific solver settings
    chrono::ChSystemParallel* system;
    double time_step;

    switch (m_method) {
        case ChContactMethod::SMC: {
            time_step = 1e-4;
            ChSystemParallelSMC* sys = new ChSystemParallelSMC;
            sys->GetSettings()->solver.contact_force_model = ChSystemSMC::Hooke;
            sys->GetSettings()->solver.tangential_displ_mode = ChSystemSMC::TangentialDisplacementModel::OneStep;
            sys->GetSettings()->solver.use_material_properties = use_mat_properties;
            system = sys;

            break;
        }
        case ChContactMethod::NSC: {
            time_step = 1e-3;
            ChSystemParallelNSC* sys = new ChSystemParallelNSC;
            sys->GetSettings()->solver.solver_mode = SolverMode::SLIDING;
            sys->GetSettings()->solver.max_iteration_normal = 0;
            sys->GetSettings()->solver.max_iteration_sliding = 200;
            sys->GetSettings()->solver.max_iteration_spinning = 0;
            sys->GetSettings()->solver.alpha = 0;
            sys->GetSettings()->solver.contact_recovery_speed = -1;
            sys->GetSettings()->collision.collision_envelope = 0.1 * radius_g;
            sys->ChangeSolverType(SolverType::APGD);
            system = sys;

            break;
        }
    }

    double g = 9.81;
    system->Set_G_acc(ChVector<>(0, 0, -g));
    system->GetSettings()->perform_thread_tuning = false;
    system->GetSettings()->solver.use_full_inertia_tensor = false;
    system->GetSettings()->solver.tolerance = 0.1;
    system->GetSettings()->solver.max_iteration_bilateral = 100;
    system->GetSettings()->collision.narrowphase_algorithm = NarrowPhaseType::NARROWPHASE_HYBRID_MPR;
    system->GetSettings()->collision.bins_per_axis = vec3(binsX, binsY, binsZ);

    // Set number of threads
    CHOMPfunctions::SetNumThreads(m_num_threads);

// Sanity check: print number of threads in a parallel region
#pragma omp parallel
#pragma omp master
    { std::cout << "Actual number of OpenMP threads: " << CHOMPfunctions::GetNumThreads() << std::endl; }

    // ---------------------
    // Create terrain bodies
    // ---------------------

    // Create contact material for terrain
    std::shared_ptr<ChMaterialSurface> material_terrain;

    switch (m_method) {
        case ChContactMethod::SMC: {
            auto mat_ter = chrono_types::make_shared<ChMaterialSurfaceSMC>();
            mat_ter->SetFriction(friction_terrain);
            mat_ter->SetRestitution(restitution_terrain);
            mat_ter->SetYoungModulus(Y_terrain);
            mat_ter->SetPoissonRatio(nu_terrain);
            mat_ter->SetAdhesion(100.0f);
            mat_ter->SetKn(kn_terrain);
            mat_ter->SetGn(gn_terrain);
            mat_ter->SetKt(kt_terrain);
            mat_ter->SetGt(gt_terrain);

            material_terrain = mat_ter;

            break;
        }
        case ChContactMethod::NSC: {
            auto mat_ter = chrono_types::make_shared<ChMaterialSurfaceNSC>();
            mat_ter->SetFriction(friction_terrain);
            mat_ter->SetRestitution(restitution_terrain);

            material_terrain = mat_ter;

            break;
        }
    }

    // Create container body
    auto container = std::shared_ptr<ChBody>(system->NewBody());
    system->AddBody(container);
    container->SetIdentifier(-1);
    container->SetMass(1);
    container->SetBodyFixed(true);
    container->SetCollide(true);

    container->GetCollisionModel()->ClearModel();
    // Bottom box
    utils::AddBoxGeometry(container.get(), material_terrain, ChVector<>(hdimX, hdimY, hthick), ChVector<>(0, 0, -hthick),
                          ChQuaternion<>(1, 0, 0, 0), true);
    // Front box
    utils::AddBoxGeometry(container.get(), material_terrain, ChVector<>(hthick, hdimY, hdimZ + hthick),
                          ChVector<>(hdimX + hthick, 0, hdimZ - hthick), ChQuaternion<>(1, 0, 0, 0), false);
    // Rear box
    utils::AddBoxGeometry(container.get(), material_terrain, ChVector<>(hthick, hdimY, hdimZ + hthick),
                          ChVector<>(-hdimX - hthick, 0, hdimZ - hthick), ChQuaternion<>(1, 0, 0, 0), false);
    // Left box
    utils::AddBoxGeometry(container.get(), material_terrain, ChVector<>(hdimX, hthick, hdimZ + hthick),
                          ChVector<>(0, hdimY + hthick, hdimZ - hthick), ChQuaternion<>(1, 0, 0, 0), false);
    // Right box
    utils::AddBoxGeometry(container.get(), material_terrain, ChVector<>(hdimX, hthick, hdimZ + hthick),
                          ChVector<>(0, -hdimY - hthick, hdimZ - hthick), ChQuaternion<>(1, 0, 0, 0), false);
    container->GetCollisionModel()->BuildModel();

    // ----------------
    // Create particles
    // ----------------

    // Create a particle generator and a mixture entirely made out of spheres
    utils::Generator gen(system);
    std::shared_ptr<utils::MixtureIngredient> m1 = gen.AddMixtureIngredient(utils::MixtureType::SPHERE, 1.0);
    m1->setDefaultMaterial(material_terrain);
    m1->setDefaultDensity(rho_g);
    m1->setDefaultSize(radius_g);

    // Set starting value for body identifiers
    gen.setBodyIdentifier(Id_g);

    // Create particles in layers until reaching the desired number of particles
    double r = 1.01 * radius_g;
    ChVector<> hdims(hdimX - r, hdimY - r, 0);
    ChVector<> center(0, 0, 2 * r);

    for (int il = 0; il < num_layers; il++) {
        gen.createObjectsBox(utils::SamplingType::POISSON_DISK, 2 * r, center, hdims);
        center.z() += 2 * r;
    }

    unsigned int num_particles = gen.getTotalNumBodies();
    std::cout << "Generated particles:  " << num_particles << std::endl;
    double total_weight = num_particles * (4 * CH_C_PI / 3) * r * r * r * rho_g * g;
    std::cout << "Total weigth:  " << total_weight << std::endl;

#ifdef CHRONO_OPENGL
    // -------------------------------
    // Create the visualization window
    // -------------------------------
    if (render) {
        opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
        gl_window.Initialize(1280, 720, "Settling test", system);
        gl_window.SetCamera(ChVector<>(0, -1, 0), ChVector<>(0, 0, 0), ChVector<>(0, 0, 1), 0.05f);
        gl_window.SetRenderMode(opengl::WIREFRAME);
    }
#endif

    // ---------------
    // Simulate system
    // ---------------

    double sim_time = 0;
    double broad_time = 0;
    double narrow_time = 0;
    double update_time = 0;
    double solve_time = 0;
    int num_steps = 0;

    ////TimingHeader();
    double time_end = 0.5;
    while (system->GetChTime() < time_end) {
        system->DoStepDynamics(time_step);

        sim_time += system->GetTimerStep();
        broad_time += system->GetTimerCollisionBroad();
        narrow_time += system->GetTimerCollisionNarrow();
        update_time += system->GetTimerUpdate();
        solve_time += system->GetTimerAdvance();
        num_steps++;

        ////TimingOutput(system);

#ifdef CHRONO_OPENGL
        if (render) {
            opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
            if (gl_window.Active()) {
                gl_window.Render();
            } else {
                return false;
            }
        }
#endif
    }

    system->CalculateContactForces();
    real3 cforce = system->GetBodyContactForce(container);
    int ncontacts = system->GetNcontacts();
    std::cout << "Number of contacts:         " << ncontacts << std::endl;
    std::cout << "Contact force on container: " << cforce.x << "  " << cforce.y << "  " << cforce.z << std::endl;
    std::cout << "Total simulation time: " << sim_time << std::endl;
    std::cout << "    Broad phase:       " << broad_time << std::endl;
    std::cout << "    Narrow phase:      " << narrow_time << std::endl;
    std::cout << "    Update phase:      " << update_time << std::endl;
    std::cout << "    Solve phase:       " << solve_time << std::endl;

    m_execTime = sim_time;
    addMetric("number_contacts", ncontacts);
    addMetric("vertical_force", cforce.z);
    addMetric("avg_sim_time_per_step (ms)", 1000 * sim_time / num_steps);
    addMetric("avg_broad_time_per_step (ms)", 1000 * broad_time / num_steps);
    addMetric("avg_narrow_time_per_step (ms)", 1000 * narrow_time / num_steps);
    addMetric("avg_update_time_per_step (ms)", 1000 * update_time / num_steps);
    addMetric("avg_solve_time_per_step (ms)", 1000 * solve_time / num_steps);

    return true;
}

int main(int argc, char** argv) {
    std::string out_dir = "../METRICS";
    if (!filesystem::create_directory(filesystem::path(out_dir))) {
        std::cout << "Error creating directory " << out_dir << std::endl;
        return 1;
    }

    bool passed = true;

    PARSettlingTest testDEM2("metrics_PAR_settling_DEM_2", "Chrono::Parallel", ChContactMethod::SMC, 2);
    PARSettlingTest testDEM4("metrics_PAR_settling_DEM_4", "Chrono::Parallel", ChContactMethod::SMC, 4);
    PARSettlingTest testDVI2("metrics_PAR_settling_DVI_2", "Chrono::Parallel", ChContactMethod::NSC, 2);
    PARSettlingTest testDVI4("metrics_PAR_settling_DVI_4", "Chrono::Parallel", ChContactMethod::NSC, 4);

    testDEM2.setOutDir(out_dir);
    testDEM2.setVerbose(true);
    passed &= testDEM2.run();
    testDEM2.print();

    testDEM4.setOutDir(out_dir);
    testDEM4.setVerbose(true);
    passed &= testDEM4.run();
    testDEM4.print();

    testDVI2.setOutDir(out_dir);
    testDVI2.setVerbose(true);
    passed &= testDVI2.run();
    testDVI2.print();

    testDVI4.setOutDir(out_dir);
    testDVI4.setVerbose(true);
    passed &= testDVI4.run();
    testDVI4.print();

    return 0;
}

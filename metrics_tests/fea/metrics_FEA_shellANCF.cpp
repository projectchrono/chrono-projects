//
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2013 Project Chrono
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file at the top level of the distribution
// and at http://projectchrono.org/license-chrono.txt.
//

#include <algorithm>
#include <iomanip>
#include <string>

#include "chrono/ChConfig.h"
#include "chrono/core/ChTimer.h"
#include "chrono/solver/ChIterativeSolverLS.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#include "chrono/fea/ChElementShellANCF.h"
#include "chrono/fea/ChLinkDirFrame.h"
#include "chrono/fea/ChLinkPointFrame.h"
#include "chrono/fea/ChMesh.h"

#include "chrono_thirdparty/filesystem/path.h"

#include "../BaseTest.h"

#ifdef CHRONO_MKL
#include "chrono_mkl/ChSolverMKL.h"
#endif

#ifdef CHRONO_MUMPS
#include "chrono_mumps/ChSolverMumps.h"
#endif

#ifdef CHRONO_OPENMP_ENABLED
#include <omp.h>
#endif

using namespace chrono;
using namespace chrono::fea;

using std::cout;
using std::endl;

enum class solver_type { MINRES, MKL, MUMPS };

// -----------------------------------------------------------------------------

int num_threads = 4;      // default number of threads
double step_size = 1e-3;  // integration step size
int num_steps = 10;       // number of integration steps

int numDiv_x = 30;  // mesh divisions in X direction
int numDiv_y = 30;  // mesh divisions in Y direction
int numDiv_z = 1;    // mesh divisions in Z direction

// -----------------------------------------------------------------------------

// Test class
class FEAShellTest : public BaseTest {
  public:
    FEAShellTest(const std::string& testName,
                 const std::string& testProjectName,
                 solver_type solver,
                 bool use_modifiedNewton,
                 bool verbose_solver)
        : BaseTest(testName, testProjectName),
          m_execTime(0),
          m_solver(solver),
          m_use_modifiedNewton(use_modifiedNewton),
          m_use_adaptiveStep(true),
          m_verbose_solver(verbose_solver) {}

    ~FEAShellTest() {}

    // Override corresponding functions in BaseTest
    virtual bool execute() override;
    virtual double getExecutionTime() const override { return m_execTime; }

  private:
    double m_execTime;
    solver_type m_solver;       // use MKL m_solver (if available)
    bool m_use_adaptiveStep;    // allow step size reduction
    bool m_use_modifiedNewton;  // use modified Newton method
    bool m_verbose_solver;      // verbose output from underlying solver
};

bool FEAShellTest::execute() {
    cout << endl;
    cout << "===================================================================" << endl;
    cout << "Solver:          ";
    switch (m_solver) {
        case solver_type::MINRES:
            cout << "MINRES";
            break;
        case solver_type::MKL:
            cout << "MKL";
            break;
        case solver_type::MUMPS:
            cout << "MUMPS";
            break;
        default:
            break;
    }
    cout << endl;
    cout << "Adaptive step:   " << (m_use_adaptiveStep ? "Yes" : "No") << endl;
    cout << "Modified Newton: " << (m_use_modifiedNewton ? "Yes" : "No") << endl;
    cout << endl;
    cout << "Mesh divisions:  " << numDiv_x << " x " << numDiv_y << endl;
    cout << endl;

    // Create the physical system
    ChSystemNSC my_system;

    my_system.Set_G_acc(ChVector<>(0, 0, -9.81));

    // Create a mesh, that is a container for groups of elements and their referenced nodes.
    auto my_mesh = chrono_types::make_shared<ChMesh>();

    // Geometry of the plate
    double plate_lenght_x = 1.0;
    double plate_lenght_y = 1.0;
    double plate_lenght_z = 0.04;  // small thickness
                                   // Specification of the mesh
    int N_x = numDiv_x + 1;
    int N_y = numDiv_y + 1;
    int N_z = numDiv_z + 1;
    // Number of elements in the z direction is considered as 1
    int TotalNumElements = numDiv_x * numDiv_y;
    //(1+1) is the number of nodes in the z direction
    int TotalNumNodes = (numDiv_x + 1) * (numDiv_y + 1);  // Or *(numDiv_z+1) for multilayer
                                                          // Element dimensions (uniform grid)
    double dx = plate_lenght_x / numDiv_x;
    double dy = plate_lenght_y / numDiv_y;
    double dz = plate_lenght_z / numDiv_z;

    // Create and add the nodes
    for (int i = 0; i < TotalNumNodes; i++) {
        // Parametric location and direction of nodal coordinates
        double loc_x = (i % (numDiv_x + 1)) * dx;
        double loc_y = (i / (numDiv_x + 1)) % (numDiv_y + 1) * dy;
        double loc_z = (i) / ((numDiv_x + 1) * (numDiv_y + 1)) * dz;

        double dir_x = 0;
        double dir_y = 0;
        double dir_z = 1;

        // Create the node
        auto node = chrono_types::make_shared<ChNodeFEAxyzD>(ChVector<>(loc_x, loc_y, loc_z), ChVector<>(dir_x, dir_y, dir_z));
        node->SetMass(0);
        // Fix all nodes along the axis X=0
        if (i % (numDiv_x + 1) == 0)
            node->SetFixed(true);

        // Add node to mesh
        my_mesh->AddNode(node);
    }

    // Create an isotropic material
    // Only one layer
    double rho = 500;
    double E = 2.1e7;
    double nu = 0.3;
    auto mat = chrono_types::make_shared<ChMaterialShellANCF>(rho, E, nu);

    // Create the elements
    for (int i = 0; i < TotalNumElements; i++) {
        // Definition of nodes forming an element
        int node0 = (i / (numDiv_x)) * (N_x) + i % numDiv_x;
        int node1 = (i / (numDiv_x)) * (N_x) + i % numDiv_x + 1;
        int node2 = (i / (numDiv_x)) * (N_x) + i % numDiv_x + 1 + N_x;
        int node3 = (i / (numDiv_x)) * (N_x) + i % numDiv_x + N_x;

        // Create the element and set its nodes.
        auto element = chrono_types::make_shared<ChElementShellANCF>();
        element->SetNodes(std::dynamic_pointer_cast<ChNodeFEAxyzD>(my_mesh->GetNode(node0)),
                          std::dynamic_pointer_cast<ChNodeFEAxyzD>(my_mesh->GetNode(node1)),
                          std::dynamic_pointer_cast<ChNodeFEAxyzD>(my_mesh->GetNode(node2)),
                          std::dynamic_pointer_cast<ChNodeFEAxyzD>(my_mesh->GetNode(node3)));

        // Element length is a fixed number in both direction. (uniform distribution of nodes in both directions)
        element->SetDimensions(dx, dy);
        // Single layer
        element->AddLayer(dz, 0 * CH_C_DEG_TO_RAD, mat);  // Thickness: dy;  Ply angle: 0.
                                                          // Set other element properties
        element->SetAlphaDamp(0.0);                       // Structural damping for this
        element->SetGravityOn(true);                      // element calculates its own gravitational load
                                                          // Add element to mesh
        my_mesh->AddElement(element);
    }

    // Switch off mesh class gravity (ANCF shell elements have a custom implementation)
    my_mesh->SetAutomaticGravity(false);

    // Remember to add the mesh to the system!
    my_system.Add(my_mesh);

#ifdef CHRONO_MKL
    std::shared_ptr<ChSolverMKL> mkl_solver;
#endif

#ifdef CHRONO_MUMPS
    std::shared_ptr<ChSolverMumps> mumps_solver;
#endif

    // Set up solver
    switch (m_solver) {
        case solver_type::MINRES: {
            auto solver = chrono_types::make_shared<ChSolverMINRES>();
            solver->SetMaxIterations(100);
            solver->EnableDiagonalPreconditioner(true);
            solver->SetVerbose(false);
            my_system.SetSolver(solver);
            my_system.SetSolverForceTolerance(1e-9);
        } break;
        case solver_type::MKL:
#ifdef CHRONO_MKL
            mkl_solver = chrono_types::make_shared<ChSolverMKL>();
            my_system.SetSolver(mkl_solver);
            mkl_solver->LockSparsityPattern(true);
            mkl_solver->SetVerbose(m_verbose_solver);
            mkl_solver->ForceSparsityPatternUpdate();
#endif
            break;
        case solver_type::MUMPS:
#ifdef CHRONO_MUMPS
            mumps_solver = chrono_types::make_shared<ChSolverMumps>();
            my_system.SetSolver(mumps_solver);
            mumps_solver->SetVerbose(m_verbose_solver);
#endif
            break;
        default:
            std::cout << "No m_solver set up" << std::endl;
            break;
    }

    // Set up integrator
    my_system.SetTimestepperType(ChTimestepper::Type::HHT);
    auto mystepper = std::static_pointer_cast<ChTimestepperHHT>(my_system.GetTimestepper());
    mystepper->SetAlpha(-0.2);
    mystepper->SetMaxiters(100);
    mystepper->SetAbsTolerances(1e-5);
    mystepper->SetMode(ChTimestepperHHT::POSITION);
    mystepper->SetStepControl(m_use_adaptiveStep);
    mystepper->SetModifiedNewton(m_use_modifiedNewton);
    mystepper->SetScaling(true);
    mystepper->SetVerbose(m_verbose_solver);

    // Initialize the output stream and set precision.
    utils::CSV_writer out("\t");
    out.stream().setf(std::ios::scientific | std::ios::showpos);
    out.stream().precision(6);

    // Get handle to tracked node.
    auto nodetip = std::dynamic_pointer_cast<ChNodeFEAxyzD>(my_mesh->GetNode(TotalNumNodes - 1));

    // Simulation loop
    double time_total = 0;
    double time_setup = 0;
    double time_solve = 0;
    double time_update = 0;
    double time_force = 0;
    double time_jacobian = 0;

    int num_iterations = 0;
    int num_setup_calls = 0;
    int num_solver_calls = 0;
    int num_force_calls = 0;
    int num_jacobian_calls = 0;

    for (int istep = 0; istep < num_steps; istep++) {
        if (m_verbose_solver) {
            cout << "-------------------------------------------------------------------" << endl;
            cout << "STEP: " << istep << endl;
        }

        my_mesh->ResetCounters();
        my_mesh->ResetTimers();

#ifdef CHRONO_MKL
        if (m_solver == solver_type::MKL)
            mkl_solver->ResetTimers();
#endif

#ifdef CHRONO_MUMPS
        if (m_solver == solver_type::MUMPS)
            mumps_solver->ResetTimers();
#endif

        my_system.DoStepDynamics(step_size);

        if (istep == 3 && m_solver == solver_type::MKL) {
#ifdef CHRONO_MKL
            mkl_solver->LockSparsityPattern(true);
#endif
        }

        time_total += my_system.GetTimerStep();
        time_setup += my_system.GetTimerSetup();
        time_solve += my_system.GetTimerSolver();
        time_update += my_system.GetTimerUpdate();

        time_force += my_mesh->GetTimeInternalForces();
        time_jacobian += my_mesh->GetTimeJacobianLoad();

        num_iterations += mystepper->GetNumIterations();
        num_setup_calls += mystepper->GetNumSetupCalls();
        num_solver_calls += mystepper->GetNumSolveCalls();

        num_force_calls += my_mesh->GetNumCallsInternalForces();
        num_jacobian_calls += my_mesh->GetNumCallsJacobianLoad();

        const ChVector<>& p = nodetip->GetPos();

        if (m_verbose_solver) {
            cout << endl;
            cout << "t = " << my_system.GetChTime() << "  ";
            cout << "node: [ " << p.x() << " " << p.y() << " " << p.z() << " ]  " << endl;
            cout << "step:  " << my_system.GetTimerStep() << endl;
            cout << "setup: " << my_system.GetTimerSetup();
            cout << endl << endl;
        }
    }

    double time_other = time_total - time_setup - time_solve - time_update - time_force - time_jacobian;

    cout << "-------------------------------------------------------------------" << endl;
    cout << "Total number of steps:        " << num_steps << endl;
    cout << "Total number of iterations:   " << num_iterations << endl;
    cout << "Total number of setup calls:  " << num_setup_calls << endl;
    cout << "Total number of m_solver calls: " << num_solver_calls << endl;
    cout << "Total number of internal force calls: " << num_force_calls << endl;
    cout << "Total number of Jacobian calls:       " << num_jacobian_calls << endl;
    cout << endl;
    cout << std::setprecision(3) << std::fixed;
    cout << "Total time: " << time_total << endl;
    cout << "  Setup:    " << time_setup << "\t (" << (time_setup / time_total) * 100 << "%)" << endl;
    cout << "  Solve:    " << time_solve << "\t (" << (time_solve / time_total) * 100 << "%)" << endl;
    cout << "  Forces:   " << time_force << "\t (" << (time_force / time_total) * 100 << "%)" << endl;
    cout << "  Jacobian: " << time_jacobian << "\t (" << (time_jacobian / time_total) * 100 << "%)" << endl;
    cout << "  Update:   " << time_update << "\t (" << (time_update / time_total) * 100 << "%)" << endl;
    cout << "  Other:    " << time_other << "\t (" << (time_other / time_total) * 100 << "%)" << endl;
    cout << endl;

    m_execTime = time_total;
    addMetric("number_iterations", num_iterations);
    addMetric("tip z displacement (mm)", 1000 * nodetip->GetPos().z());
    addMetric("time_setup", time_setup);
    addMetric("time_solve", time_solve);
    addMetric("time_jacobian", time_jacobian);
    addMetric("time_force", time_force);

    return true;
}

int main(int argc, char* argv[]) {
    std::string out_dir = "../METRICS";
    if (!filesystem::create_directory(filesystem::path(out_dir))) {
        std::cout << "Error creating directory " << out_dir << std::endl;
        return 1;
    }
#ifdef CHRONO_OPENMP_ENABLED
    // Set number of threads
    if (argc > 1)
        num_threads = std::stoi(argv[1]);
    num_threads = std::min(num_threads, CHOMPfunctions::GetNumProcs());
    CHOMPfunctions::SetNumThreads(num_threads);
    GetLog() << "Using " << num_threads << " thread(s)\n";
#else
    GetLog() << "No OpenMP\n";
#endif

    bool verbose_solver = false;
    bool verbose_test = false;

    FEAShellTest test_minres_full("metrics_FEA_shellANCF_MINRES_full", "Chrono::FEA", solver_type::MINRES, false, verbose_solver);
    test_minres_full.setOutDir(out_dir);
    test_minres_full.setVerbose(verbose_test);
    test_minres_full.run();
    test_minres_full.print();

    FEAShellTest test_minres_mod("metrics_FEA_shellANCF_MINRES_modified", "Chrono::FEA", solver_type::MINRES, true, verbose_solver);
    test_minres_mod.setOutDir(out_dir);
    test_minres_mod.setVerbose(verbose_test);
    test_minres_mod.run();
    test_minres_mod.print();

#ifdef CHRONO_MKL
    FEAShellTest test_mkl_full("metrics_FEA_shellANCF_MKL_full", "Chrono::FEA", solver_type::MKL, false, verbose_solver);
    test_mkl_full.setOutDir(out_dir);
    test_mkl_full.setVerbose(verbose_test);
    test_mkl_full.run();
    test_mkl_full.print();
     
    FEAShellTest test_mkl_mod("metrics_FEA_shellANCF_MKL_modified", "Chrono::FEA", solver_type::MKL, true, verbose_solver);
    test_mkl_mod.setOutDir(out_dir);
    test_mkl_mod.setVerbose(verbose_test);
    test_mkl_mod.run();
    test_mkl_mod.print();
#endif

    return 0;
}

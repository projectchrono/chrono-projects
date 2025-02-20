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
#include "chrono/utils/ChOpenMP.h"

#include "chrono/fea/ChElementShellANCF_3423.h"
#include "chrono/fea/ChLinkNodeSlopeFrame.h"
#include "chrono/fea/ChLinkNodeFrame.h"
#include "chrono/fea/ChMesh.h"

#include "chrono_thirdparty/filesystem/path.h"

#ifdef CHRONO_PARDISO_MKL
#include "chrono_pardisomkl/ChSolverPardisoMKL.h"
#endif

#undef CHRONO_MUMPS
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

// -----------------------------------------------------------------------------

int num_threads = 4;      // default number of threads
double step_size = 1e-3;  // integration step size
int num_steps = 20;       // number of integration steps
int skip_steps = 0;       // initial number of steps excluded from timing

int numDiv_x = 100;  // mesh divisions in X direction
int numDiv_y = 100;  // mesh divisions in Y direction
int numDiv_z = 1;   // mesh divisions in Z direction

std::string out_dir = "../TEST_SHELL_ANCF";  // name of output directory
bool output = true;                         // generate output file?
bool verbose = true;                         // verbose output?

// -----------------------------------------------------------------------------

void RunModel(int nthreads,              // number of OpenMP threads
              ChSolver::Type solver,     // linear solver type
              bool use_adaptiveStep,     // allow step size reduction
              bool use_modifiedNewton,   // use modified Newton method
              const std::string& suffix  // output filename suffix
) {
    cout << endl;
    cout << "===================================================================" << endl;
    cout << "Solver:          ";
    switch (solver) {
        case ChSolver::Type::MINRES:
            cout << "MINRES";
            break;
        case ChSolver::Type::PARDISO_MKL:
            cout << "PardisoMKL";
            break;
        case ChSolver::Type::MUMPS:
            cout << "MUMPS";
            break;
        default:
            break;
    }
    cout << endl;
    cout << "Adaptive step:   " << (use_adaptiveStep ? "Yes" : "No") << endl;
    cout << "Modified Newton: " << (use_modifiedNewton ? "Yes" : "No") << endl;
    cout << endl;
    cout << "Mesh divisions:  " << numDiv_x << " x " << numDiv_y << endl;
    cout << endl;

    // Create the physical system
    ChSystemNSC my_system;
    my_system.SetNumThreads(nthreads);
    my_system.SetGravitationalAcceleration(ChVector3d(0, 0, -9.81));

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
        auto node = chrono_types::make_shared<ChNodeFEAxyzD>(ChVector3d(loc_x, loc_y, loc_z), ChVector3d(dir_x, dir_y, dir_z));
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
        int node0 = (i / (numDiv_x)) * (N_x)+i % numDiv_x;
        int node1 = (i / (numDiv_x)) * (N_x)+i % numDiv_x + 1;
        int node2 = (i / (numDiv_x)) * (N_x)+i % numDiv_x + 1 + N_x;
        int node3 = (i / (numDiv_x)) * (N_x)+i % numDiv_x + N_x;

        // Create the element and set its nodes.
        auto element = chrono_types::make_shared<ChElementShellANCF_3423>();
        element->SetNodes(std::dynamic_pointer_cast<ChNodeFEAxyzD>(my_mesh->GetNode(node0)),
                          std::dynamic_pointer_cast<ChNodeFEAxyzD>(my_mesh->GetNode(node1)),
                          std::dynamic_pointer_cast<ChNodeFEAxyzD>(my_mesh->GetNode(node2)),
                          std::dynamic_pointer_cast<ChNodeFEAxyzD>(my_mesh->GetNode(node3)));

        // Element length is a fixed number in both direction. (uniform distribution of nodes in both directions)
        element->SetDimensions(dx, dy);
        element->AddLayer(dz, 0 * CH_DEG_TO_RAD, mat);  // Single layer; Thickness: dy;  Ply angle: 0.
        element->SetAlphaDamp(0.0);                       // Structural damping for this
        my_mesh->AddElement(element);
    }

    // Switch off mesh class gravity (ANCF shell elements have a custom implementation)
    my_mesh->SetAutomaticGravity(false);

    // Remember to add the mesh to the system!
    my_system.Add(my_mesh);

#ifdef CHRONO_PARDISO_MKL
    std::shared_ptr<ChSolverPardisoMKL> mkl_solver;
#endif

#ifdef CHRONO_MUMPS
    std::shared_ptr<ChSolverMumps> mumps_solver;
#endif


    // Set up solver
    switch (solver) {
        case ChSolver::Type::MINRES: {
            auto solver = chrono_types::make_shared<ChSolverMINRES>();
            solver->EnableDiagonalPreconditioner(true);
            solver->SetMaxIterations(100);
            solver->SetVerbose(true);
            solver->SetTolerance(1e-12);
            my_system.SetSolver(solver);
        } break;
        case ChSolver::Type::PARDISO_MKL:
#ifdef CHRONO_PARDISO_MKL
            mkl_solver = chrono_types::make_shared<ChSolverPardisoMKL>();
            my_system.SetSolver(mkl_solver);
            mkl_solver->LockSparsityPattern(true);
            mkl_solver->SetVerbose(verbose);
            mkl_solver->ForceSparsityPatternUpdate();
#endif
            break;
        case ChSolver::Type::MUMPS:
#ifdef CHRONO_MUMPS
            mumps_solver = chrono_types::make_shared<ChSolverMumps>();
            my_system.SetSolver(mumps_solver);
            mumps_solver->SetVerbose(verbose);
#endif
            break;
        default:
            std::cout << "No solver set up" << std::endl;
            break;
    }

    // Set up integrator
    my_system.SetTimestepperType(ChTimestepper::Type::HHT);
    auto mystepper = std::static_pointer_cast<ChTimestepperHHT>(my_system.GetTimestepper());
    mystepper->SetAlpha(-0.2);
    mystepper->SetMaxIters(100);
    mystepper->SetAbsTolerances(1e-3);
    mystepper->SetStepControl(use_adaptiveStep);
    mystepper->SetModifiedNewton(use_modifiedNewton);
    mystepper->SetVerbose(verbose);

    // Initialize the output stream and set precision.
    utils::ChWriterCSV out("\t");
    out.Stream().setf(std::ios::scientific | std::ios::showpos);
    out.Stream().precision(6);

    // Get handle to tracked node.
    auto nodetip = std::dynamic_pointer_cast<ChNodeFEAxyzD>(my_mesh->GetNode(TotalNumNodes - 1));

    // Simulation loop
    double time_total = 0;
    double time_setup = 0;
    double time_setup_assembly = 0;
    double time_setup_solvercall = 0;
    double time_solve = 0;
    double time_solve_assembly = 0;
    double time_solve_solvercall = 0;
    double time_update = 0;
    double time_force = 0;
    double time_jacobian = 0;
    double time_skipped = 0;

    int num_iterations = 0;
    int num_setup_calls = 0;
    int num_solver_calls = 0;
    int num_force_calls = 0;
    int num_jacobian_calls = 0;

    for (int istep = 0; istep < num_steps; istep++) {
        if (verbose) {
            cout << "-------------------------------------------------------------------" << endl;
            cout << "STEP: " << istep << endl;
        }

        my_mesh->ResetCounters();
        my_mesh->ResetTimers();

#ifdef CHRONO_PARDISO_MKL
        if (solver == ChSolver::Type::PARDISO_MKL)
            mkl_solver->ResetTimers();
#endif

#ifdef CHRONO_MUMPS
        if (solver == ChSolver::Type::MUMPS)
            mumps_solver->ResetTimers();
#endif

        my_system.DoStepDynamics(step_size);

        if (istep == 3 && solver == ChSolver::Type::PARDISO_MKL)
        {
#ifdef CHRONO_PARDISO_MKL
            mkl_solver->LockSparsityPattern(true);
#endif
        }

        if (istep == skip_steps) {
            if (verbose)
                cout << "Resetting counters at step = " << istep << endl;
            time_skipped = time_total;
            time_total = 0;
            time_setup = 0;
            time_setup_assembly = 0;
            time_setup_solvercall = 0;
            time_solve = 0;
            time_solve_assembly = 0;
            time_solve_solvercall = 0;
            time_update = 0;
            time_force = 0;
            time_jacobian = 0;
            num_iterations = 0;
            num_setup_calls = 0;
            num_solver_calls = 0;
            num_force_calls = 0;
            num_jacobian_calls = 0;
        }

        time_total += my_system.GetTimerStep();
        time_setup += my_system.GetTimerSetup();
        time_solve += my_system.GetTimerLSsolve();
        time_update += my_system.GetTimerUpdate();

        // TODO: if it is OK to move timer in ChSolver we can avoid this switch
#ifdef CHRONO_PARDISO_MKL
        if (solver == ChSolver::Type::PARDISO_MKL) {
            time_setup_assembly += mkl_solver->GetTimeSetup_Assembly();
            time_setup_solvercall += mkl_solver->GetTimeSetup_SolverCall();
            time_solve_assembly += mkl_solver->GetTimeSolve_Assembly();
            time_solve_solvercall += mkl_solver->GetTimeSolve_SolverCall();
        }
#endif
#ifdef CHRONO_MUMPS
        if (solver == ChSolver::Type::MUMPS) {
            time_setup_assembly += mumps_solver->GetTimeSetup_Assembly();
            time_setup_solvercall += mumps_solver->GetTimeSetup_SolverCall();
            time_solve_assembly += mumps_solver->GetTimeSolve_Assembly();
            time_solve_solvercall += mumps_solver->GetTimeSolve_SolverCall();
        }
#endif
        time_force += my_mesh->GetTimeInternalForces();
        time_jacobian += my_mesh->GetTimeJacobianLoad();

        num_iterations += mystepper->GetNumIterations();
        num_setup_calls += mystepper->GetNumSetupCalls();
        num_solver_calls += mystepper->GetNumSolveCalls();

        num_force_calls += my_mesh->GetNumCallsInternalForces();
        num_jacobian_calls += my_mesh->GetNumCallsJacobianLoad();

        const ChVector3d& p = nodetip->GetPos();

        if (verbose) {
            cout << endl;
            cout << "t = " << my_system.GetChTime() << "  ";
            cout << "node: [ " << p.x() << " " << p.y() << " " << p.z() << " ]  " << endl;
            cout << "step:  " << my_system.GetTimerStep() << endl;
            cout << "setup: " << my_system.GetTimerSetup();
#ifdef CHRONO_PARDISO_MKL
            if (solver == ChSolver::Type::PARDISO_MKL) {
                cout << "  [assembly: " << mkl_solver->GetTimeSetup_Assembly();
                cout << "  pardiso: " << mkl_solver->GetTimeSetup_SolverCall() <<"]";
            }
#endif
            cout << endl;
            cout << "solve: " << my_system.GetTimerLSsolve() << "  ";
#ifdef CHRONO_MUMPS
            if (solver == ChSolver::Type::MUMPS) {
                cout << "  [assembly: " << mumps_solver->GetTimeSolve_Assembly();
                cout << "  mumps: " << mumps_solver->GetTimeSolve_SolverCall() << "]";
            }
#endif
            cout << endl << endl;
        }

        if (output) {
            out << my_system.GetChTime() << my_system.GetTimerStep() << nodetip->GetPos() << endl;
        }
    }

    double time_other = time_total - time_setup - time_solve - time_update - time_force - time_jacobian;

    cout << "-------------------------------------------------------------------" << endl;
    cout << "Total number of steps:        " << num_steps - skip_steps << endl;
    cout << "Total number of iterations:   " << num_iterations << endl;
    cout << "Total number of setup calls:  " << num_setup_calls << endl;
    cout << "Total number of solver calls: " << num_solver_calls << endl;
    cout << "Total number of internal force calls: " << num_force_calls << endl;
    cout << "Total number of Jacobian calls:       " << num_jacobian_calls << endl;
    cout << endl;
    cout << std::setprecision(3) << std::fixed;
    cout << "Total time: " << time_total << endl;
    cout << "  Setup:    " << time_setup << "\t (" << (time_setup / time_total) * 100 << "%)" << endl;
    if (solver == ChSolver::Type::PARDISO_MKL || solver == ChSolver::Type::MUMPS) {
        cout << "    Assembly: " << time_setup_assembly << "\t (" << (time_setup_assembly / time_setup) * 100
            << "% setup)" << endl;
        cout << "    SolverCall:  " << time_setup_solvercall << "\t (" << (time_setup_solvercall / time_setup) * 100
            << "% setup)" << endl;
    }
    cout << "  Solve:    " << time_solve << "\t (" << (time_solve / time_total) * 100 << "%)" << endl;
    if (solver == ChSolver::Type::PARDISO_MKL || solver == ChSolver::Type::MUMPS) {
        cout << "    Assembly: " << time_solve_assembly << "\t (" << (time_solve_assembly / time_solve) * 100
            << "% solve)" << endl;
        cout << "    SolverCall:  " << time_solve_solvercall << "\t (" << (time_solve_solvercall / time_solve) * 100
            << "% solve)" << endl;
    }
    if (solver == ChSolver::Type::PARDISO_MKL || solver == ChSolver::Type::MUMPS) {
        cout << "  [TOT Assembly: " << time_setup_assembly+time_solve_assembly << "\t (" << ((time_setup_assembly + time_solve_assembly) / time_total) * 100
            << "% total)]" << endl;
        cout << "  [TOT SolverCall:  " << time_setup_solvercall + time_solve_solvercall << "\t (" << ((time_setup_solvercall + time_solve_solvercall) / time_total) * 100
            << "% total)]" << endl;
    }
    cout << "  Forces:   " << time_force << "\t (" << (time_force / time_total) * 100 << "%)" << endl;
    cout << "  Jacobian: " << time_jacobian << "\t (" << (time_jacobian / time_total) * 100 << "%)" << endl;
    cout << "  Update:   " << time_update << "\t (" << (time_update / time_total) * 100 << "%)" << endl;
    cout << "  Other:    " << time_other << "\t (" << (time_other / time_total) * 100 << "%)" << endl;
    cout << endl;
    cout << "Time for skipped steps (" << skip_steps << "): " << time_skipped << endl;

    if (output) {
        char name[100];
        std::sprintf(name, "%s/out_%s_%d.txt", out_dir.c_str(), suffix.c_str(), num_threads);
        cout << "Write output to: " << name << endl;
        out.WriteToFile(name);
    }
}

int main(int argc, char* argv[]) {
    // Set path to Chrono data directory
    SetChronoDataPath(CHRONO_DATA_DIR);

    // Create output directory (if it does not already exist).
    if (output) {
        if (!filesystem::create_directory(filesystem::path("../TEST_SHELL_ANCF"))) {
            std::cout << "Error creating directory ../TEST_SHELL_ANCF\n";
            return 1;
        }
    }

#ifdef CHRONO_OPENMP_ENABLED
    // Set number of threads
    if (argc > 1)
        num_threads = std::stoi(argv[1]);
    num_threads = std::min(num_threads, ChOMP::GetNumProcs());
    std::cout << "Using " << num_threads << " thread(s)\n";
#else
    std::cout << "No OpenMP\n";
#endif

    // Run simulations.
#ifdef CHRONO_PARDISO_MKL
    RunModel(num_threads, ChSolver::Type::PARDISO_MKL, true, false, "PardisoMKL_adaptive_full");
    RunModel(num_threads, ChSolver::Type::PARDISO_MKL, true, true, "PardisoMKL_adaptive_modified");
#endif

#ifdef CHRONO_MUMPS
    RunModel(num_threads, ChSolver::Type::MUMPS, true, false, "MUMPS_adaptive_full"); 
    RunModel(num_threads, ChSolver::Type::MUMPS, true, true, "MUMPS_adaptive_modified");
#endif

    RunModel(num_threads, ChSolver::Type::MINRES, true, false, "MINRES_adaptive_full");
    RunModel(num_threads, ChSolver::Type::MINRES, true, true, "MINRES_adaptive_modified");

    return 0;
}

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
// Authors: Radu Serban
// =============================================================================
//
// Tests for mesh-body constraints
//
// =============================================================================

#include <cstdio>
#include <cmath>

#include "chrono/physics/ChLoadContainer.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/solver/ChIterativeSolverLS.h"
#include "chrono/timestepper/ChTimestepper.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#include "chrono/fea/ChElementCableANCF.h"
#include "chrono/fea/ChLinkDirFrame.h"
#include "chrono/fea/ChLinkPointFrame.h"
#include "chrono/fea/ChLoadsBeam.h"
#include "chrono/fea/ChMesh.h"

////#undef CHRONO_PARDISO_MKL
#ifdef CHRONO_PARDISO_MKL
#include "chrono_pardisomkl/ChSolverPardisoMKL.h"
#endif

#undef CHRONO_POSTPROCESS
#ifdef CHRONO_POSTPROCESS
#include "chrono_postprocess/ChGnuPlot.h"
#endif

using namespace chrono;
using namespace chrono::fea;

// Test the reaction force and torque for an ANCF beam element with two nodes.
// The first node is at the origin and is constrained to ground using a ChLinkPointFrame
// constraint and, optionally, a ChLinkDirFrame constraint (as specified through the
// boolean argument 'constrain_dir').
// The argument 'dir' defines the initial beam configuration.
void test_beam(const std::string& name,  /// test name
               const ChVector3d& dir,    /// initial beam orientation
               double alpha,             /// damping coefficient
               bool constrain_dir,       /// if true, include ChLinkFirFrame constraints
               double duration           /// simulation length
               ) {
    std::cout << "=== Test " << name << " ===" << std::endl;
    std::cout << "Beam direction: " << dir.x() << " " << dir.y() << " " << dir.z() << std::endl;
    std::cout << "Constrain direction? " << constrain_dir << std::endl;

    // Create the system
    ChSystemNSC my_system;
    double g = 10;
    my_system.SetGravitationalAcceleration(ChVector3d(0, 0, -g));

    // Create the ground body
    auto ground = chrono_types::make_shared<ChBody>();
    ground->SetFixed(true);
    my_system.AddBody(ground);

    // Create the FEA mesh
    auto mesh = chrono_types::make_shared<ChMesh>();
    mesh->SetAutomaticGravity(true);
    my_system.Add(mesh);

    const double length = 1.0;
    double diam = 5e-3;
    double rho = 1000;
    auto msection_cable = chrono_types::make_shared<ChBeamSectionCable>();
    msection_cable->SetDiameter(diam);
    msection_cable->SetYoungModulus(1e7);
    msection_cable->SetInertia(CH_C_PI / 4.0 * pow(diam / 2, 4));
    msection_cable->SetDensity(rho);

    // Create the 2 nodes and 1 beam element
    auto node1 = chrono_types::make_shared<ChNodeFEAxyzD>(ChVector3d(0, 0, 0), dir);
    auto node2 = chrono_types::make_shared<ChNodeFEAxyzD>(length * dir, dir);
    mesh->AddNode(node1);
    mesh->AddNode(node2);

    auto beam_elem = chrono_types::make_shared<ChElementCableANCF>();
    beam_elem->SetNodes(node1, node2);
    beam_elem->SetSection(msection_cable);
    beam_elem->SetAlphaDamp(alpha);
    mesh->AddElement(beam_elem);

    // Create a hinge constraint
    auto point_cnstr = chrono_types::make_shared<ChLinkPointFrame>();
    point_cnstr->Initialize(node1, ground);
    my_system.Add(point_cnstr);

    // Create a direction constraint
    std::shared_ptr<ChLinkDirFrame> dir_cnstr;
    if (constrain_dir) {
        dir_cnstr = chrono_types::make_shared<ChLinkDirFrame>();
        dir_cnstr->Initialize(node1, ground);
        my_system.Add(dir_cnstr);
    }

#ifdef CHRONO_PARDISO_MKL
    // PardisoMKL solver + HHT
    std::cout << "Using HHT + PardisoMKL" << std::endl;

    auto mkl_solver = chrono_types::make_shared<ChSolverPardisoMKL>();
    my_system.SetSolver(mkl_solver);
    mkl_solver->LockSparsityPattern(true);
    mkl_solver->SetVerbose(false);

    my_system.SetTimestepperType(ChTimestepper::Type::HHT);
    auto mystepper = std::static_pointer_cast<ChTimestepperHHT>(my_system.GetTimestepper());
    mystepper->SetAlpha(-0.2);
    mystepper->SetMaxiters(100);
    mystepper->SetAbsTolerances(1e-3);
    mystepper->SetStepControl(true);
    mystepper->SetModifiedNewton(true);
    mystepper->SetVerbose(false);
#else
    // MINRES solver + Euler
    std::cout << "Using Euler + MINRES" << std::endl;
    auto solver = chrono_types::make_shared<ChSolverMINRES>();
    solver->EnableDiagonalPreconditioner(true);
    solver->EnableWarmStart(true);
    solver->SetMaxIterations(200);
    solver->SetTolerance(1e-12);
    solver->SetVerbose(false);
    my_system.SetSolver(solver);

    // my_system.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT_LINEARIZED);
    my_system.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT);
#endif

    // Simulation loop
    double step = 1e-4;

    utils::ChWriterCSV csv("  ");
    ChVector3d rforce(0);
    ChVector3d rtorque(0);

    while (my_system.GetChTime() < duration) {
        my_system.DoStepDynamics(step);

        ChVector3d pos = node2->GetPos();

        {
            auto frame = point_cnstr->GetFrameAbs();
            rforce = point_cnstr->GetReactionOnBody();
            rforce = frame.TransformDirectionLocalToParent(rforce);
        }

        if (constrain_dir) {
            auto frame = dir_cnstr->GetFrameAbs();
            rtorque = dir_cnstr->GetReactionOnBody();
            rtorque = frame.TransformDirectionLocalToParent(rtorque);
        }

        csv << my_system.GetChTime() << pos << rforce << rtorque << std::endl;
    }

    ChVector3d dir1 = node1->GetSlope1();
    ChVector3d pos2 = node2->GetPos();
    std::cout << "\nFinal configuration" << std::endl;
    std::cout << "  base direction: " << dir1.x() << " " << dir1.y() << " " << dir1.z() << std::endl;
    std::cout << "  tip position:   " << pos2.x() << " " << pos2.y() << " " << pos2.z() << std::endl;

    std::cout << "\nReaction force and torque on ground (final)" << std::endl;
    std::cout << "  force:  " << rforce.x() << " " << rforce.y() << " " << rforce.z() << std::endl;
    std::cout << "  torque: " << rtorque.x() << " " << rtorque.y() << " " << rtorque.z() << std::endl;

    // Final beam configuration
    int num_points = 51;
    ChVectorDynamic<> displ(beam_elem->GetNumCoordsPosLevel());
    beam_elem->GetStateBlock(displ);
    std::vector<ChVector3d> P(num_points);
    for (int i = 0; i < num_points; i++) {
        double eta = -1.0 + ((2.0 * i) / (num_points - 1));
        ChQuaternion<> R;
        beam_elem->EvaluateSectionFrame(eta, P[i], R);
    }

    // Estimate reaction force & torque for an equivalent rigid beam.
    double mass = CH_C_PI_4 * diam * diam * length * rho;
    rforce = ChVector3d(0, 0, -mass * g);
    rtorque = ChVector3d(0);
    if (constrain_dir) {
        double seg_mass = mass / num_points;
        ChVector3d seg_weight(0, 0, -seg_mass * g);
        for (int i = 0; i < num_points; i++) {
            rtorque += Vcross(P[i], seg_weight);
        }
    }

    std::cout << "\nRigid body approximation" << std::endl;
    std::cout << "  force:  " << rforce.x() << " " << rforce.y() << " " << rforce.z() << std::endl;
    std::cout << "  torque: " << rtorque.x() << " " << rtorque.y() << " " << rtorque.z() << std::endl;
    std::cout << std::endl << std::endl;

#ifdef CHRONO_POSTPROCESS
    std::string out_file = name + ".out";
    csv.WriteToFile(out_file);

    {
        std::string gpl_file = name + "_forces.gpl";
        std::string title = name + ": reaction forces on ground";
        postprocess::ChGnuPlot mplot(gpl_file.c_str());
        mplot.SetGrid();
        mplot.SetTitle(title.c_str());
        mplot.SetLabelX("t");
        mplot.Plot(out_file.c_str(), 1, 5, "F_x", " with lines");
        mplot.Plot(out_file.c_str(), 1, 6, "F_y", " with lines");
        mplot.Plot(out_file.c_str(), 1, 7, "F_z", " with lines");
    }

    {
        std::string gpl_file = name + "_torques.gpl";
        std::string title = name + ": reaction torques on ground";
        postprocess::ChGnuPlot mplot(gpl_file.c_str());
        mplot.SetGrid();
        mplot.SetTitle(title.c_str());
        mplot.SetLabelX("t");
        mplot.Plot(out_file.c_str(), 1, 8, "T_x", " with lines");
        mplot.Plot(out_file.c_str(), 1, 9, "T_y", " with lines");
        mplot.Plot(out_file.c_str(), 1, 10, "T_z", " with lines");
    }

    {
        utils::ChWriterCSV csv(" ");
        for (int i = 0; i < num_points; i++)
            csv << P[i] << std::endl;
        std::string out_file = name + "_beam.out";
        csv.WriteToFile(out_file);

        std::string gpl_file = name + "_beam.gpl";
        std::string title = name + ": final beam configuration";
        postprocess::ChGnuPlot mplot(gpl_file.c_str());
        mplot.SetGrid();
        mplot.SetCommand("set size ratio -1");
        mplot.SetTitle(title.c_str());
        mplot.SetLabelX("x");
        mplot.SetLabelY("z");
        mplot.Plot(out_file.c_str(), 1, 3, "", " with lines");
    }
#endif
}

int main(int argc, char* argv[]) {
    // Set path to Chrono data directory
    SetChronoDataPath(CHRONO_DATA_DIR);

    test_beam("hinge", ChVector3d(0, 0, -1), 0.005, false, 0.75);

    test_beam("cantilever1", ChVector3d(1, 0, 0), 0.1, true, 2.5);

    ChVector3d dir(1, 0, -1);
    dir.Normalize();
    test_beam("cantilever2", dir, 0.15, true, 2.5);

    return 0;
}

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
// Authors: Antonio Recuero, Milad Rakhsha, Conlain Kelly, Radu Serban
// =============================================================================
//
// EAS Brick Element
//
// This test checks the dynamics of a beam made up of 10 brick elements.
// This element is a regular 8-noded trilinear brick element with enhanced
// assumed strain that alleviates locking. More information on the validation of
// this element may be found in Chrono's documentation. This simulation excites
// the beam by applying the sudden action of a gravity field.
// =============================================================================

#include <cmath>
#include <algorithm>
#include <string>

#include "chrono/physics/ChSystemNSC.h"
#include "chrono/solver/ChIterativeSolverLS.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#include "chrono/fea/ChElementBar.h"
#include "chrono/fea/ChElementBrick.h"
#include "chrono/fea/ChElementSpring.h"
#include "chrono/fea/ChLinkDirFrame.h"
#include "chrono/fea/ChLinkPointFrame.h"
#include "chrono/fea/ChVisualizationFEAmesh.h"

#include "chrono_thirdparty/filesystem/path.h"

#include "../BaseTest.h"

#undef CHRONO_MKL

#ifdef CHRONO_MKL
#include "chrono_mkl/ChSolverMKL.h"
#endif

using namespace chrono;
using namespace fea;

bool use_mkl = true;            // Use the MKL solver (if available)
const double step_size = 1e-3;  // Step size
const int num_steps = 500;      // Number of time steps for test

// ====================================================================================

// Test class
class BrickIso_GravTest : public BaseTest {
  public:
    BrickIso_GravTest(const std::string& testName, const std::string& testProjectName)
        : BaseTest(testName, testProjectName), m_execTime(0) {}

    ~BrickIso_GravTest() {}

    // Override corresponding functions in BaseTest
    virtual bool execute() override;
    virtual double getExecutionTime() const override { return m_execTime; }

  private:
    double m_execTime;
};

bool BrickIso_GravTest::execute() {
    // Create the physical system
    ChSystemNSC my_system;
    my_system.Set_G_acc(ChVector<>(0, 0, -9.81));

    auto my_mesh = chrono_types::make_shared<ChMesh>();

    // Dimensions of the plate
    double plate_lenght_x = 1;
    double plate_lenght_y = 0.1;
    double plate_lenght_z = 0.01;  // small thickness
    // Specification of the mesh
    int numDiv_x = 10;  // 10 elements along the length of the beam
    int numDiv_y = 1;
    int numDiv_z = 1;
    int N_x = numDiv_x + 1;
    int N_y = numDiv_y + 1;
    int N_z = numDiv_z + 1;
    // Number of elements in the z direction is considered as 1
    int TotalNumElements = numDiv_x;
    int TotalNumNodes = (numDiv_x + 1) * 4;  // 4 nodes per brick face
    // For uniform mesh
    double dx = plate_lenght_x / numDiv_x;
    double dy = plate_lenght_y / numDiv_y;
    double dz = plate_lenght_z / numDiv_z;
    int MaxMNUM = 1;
    int MTYPE = 1;
    int MaxLayNum = 1;
    ChMatrixDynamic<double> COORDFlex(TotalNumNodes, 3);
    ChMatrixDynamic<double> VELCYFlex(TotalNumNodes, 3);
    ChMatrixDynamic<int> NumNodes(TotalNumElements, 8);
    ChMatrixDynamic<int> LayNum(TotalNumElements, 1);
    ChMatrixDynamic<int> NDR(TotalNumNodes, 3);
    ChMatrixDynamic<double> ElemLengthXY(TotalNumElements, 3);
    ChMatrixNM<double, 10, 12> MPROP;

    //!------------ Material Data-----------------

    for (int i = 0; i < MaxMNUM; i++) {
        MPROP(i, 0) = 500;      // Density [kg/m3]
        MPROP(i, 1) = 2.1E+07;  // E (Pa)
        MPROP(i, 2) = 0.3;      // nu
    }

    std::shared_ptr<ChContinuumElastic> mmaterial(new ChContinuumElastic);
    mmaterial->Set_RayleighDampingK(0.0);
    mmaterial->Set_RayleighDampingM(0.0);
    mmaterial->Set_density(MPROP(0, 0));
    mmaterial->Set_E(MPROP(0, 1));
    mmaterial->Set_G(MPROP(0, 1) / (2 + 2 * MPROP(0, 2)));
    mmaterial->Set_v(MPROP(0, 2));

    //!--------------- Element data--------------------

    for (int i = 0; i < TotalNumElements; i++) {
        // All the elements belong to the same layer, e.g layer number 1.
        LayNum(i, 0) = 1;

        NumNodes(i, 0) = i;
        NumNodes(i, 1) = i + 1;
        NumNodes(i, 2) = i + 1 + N_x;
        NumNodes(i, 3) = i + N_x;
        NumNodes(i, 4) = i + 2 * N_x;
        NumNodes(i, 5) = i + 2 * N_x + 1;
        NumNodes(i, 6) = i + 3 * N_x + 1;
        NumNodes(i, 7) = i + 3 * N_x;  // Numbering of nodes

        ElemLengthXY(i, 0) = dx;
        ElemLengthXY(i, 1) = dy;
        ElemLengthXY(i, 2) = dz;

        if (MaxLayNum < LayNum(i, 0)) {
            MaxLayNum = LayNum(i, 0);
        }
    }

    //!--------- Constraints and initial coordinates/velocities

    for (int i = 0; i < TotalNumNodes; i++) {
        // Constrain clamped nodes. Assigns (1) if constrained, (0) otherwise
        NDR(i, 0) = (i % N_x == 0) ? 1 : 0;
        NDR(i, 1) = (i % N_x == 0) ? 1 : 0;
        NDR(i, 2) = (i % N_x == 0) ? 1 : 0;

        // Coordinates
        COORDFlex(i, 0) = (i % (N_x)) * dx;
        if ((i / N_x + 1) % 2 == 0) {
            COORDFlex(i, 1) = dy;
        } else {
            COORDFlex(i, 1) = 0.0;
        };
        COORDFlex(i, 2) = (i) / ((N_x)*2) * dz;

        // Velocities
        VELCYFlex(i, 0) = 0;
        VELCYFlex(i, 1) = 0;
        VELCYFlex(i, 2) = 0;
    }

    // Adding the nodes to the mesh
    int i = 0;
    while (i < TotalNumNodes) {
        auto node = chrono_types::make_shared<ChNodeFEAxyz>(ChVector<>(COORDFlex(i, 0), COORDFlex(i, 1), COORDFlex(i, 2)));
        node->SetMass(0.0);
        my_mesh->AddNode(node);
        if (NDR(i, 0) == 1 && NDR(i, 1) == 1 && NDR(i, 2) == 1) {
            node->SetFixed(true);
        }
        i++;
    }

    // Create a node at the tip by dynamic casting
    auto nodetip = std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(TotalNumNodes - 1));

    int elemcount = 0;
    while (elemcount < TotalNumElements) {
        auto element = chrono_types::make_shared<ChElementBrick>();
        ChVectorN<double, 3> InertFlexVec;  // Read element length, used in ChElementBrick
        InertFlexVec.setZero();
        InertFlexVec(0) = ElemLengthXY(elemcount, 0);
        InertFlexVec(1) = ElemLengthXY(elemcount, 1);
        InertFlexVec(2) = ElemLengthXY(elemcount, 2);
        element->SetInertFlexVec(InertFlexVec);
        // Note we change the order of the nodes to comply with the arrangement of shape functions
        element->SetNodes(std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(NumNodes(elemcount, 0))),
                          std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(NumNodes(elemcount, 1))),
                          std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(NumNodes(elemcount, 2))),
                          std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(NumNodes(elemcount, 3))),
                          std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(NumNodes(elemcount, 4))),
                          std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(NumNodes(elemcount, 5))),
                          std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(NumNodes(elemcount, 6))),
                          std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(NumNodes(elemcount, 7))));
        element->SetMaterial(mmaterial);
        element->SetElemNum(elemcount);
        element->SetGravityOn(true);      // Turn gravity on/off from within the element
        element->SetMooneyRivlin(false);  // Turn on/off Mooney Rivlin (Linear Isotropic by default)
        // element->SetMRCoefficients(551584.0, 137896.0); // Set two coefficients for Mooney-Rivlin

        ChVectorN<double, 9> stock_alpha_EAS;
        stock_alpha_EAS.setZero();
        element->SetStockAlpha(stock_alpha_EAS(0), stock_alpha_EAS(1), stock_alpha_EAS(2),
                               stock_alpha_EAS(3), stock_alpha_EAS(4), stock_alpha_EAS(5),
                               stock_alpha_EAS(6), stock_alpha_EAS(7), stock_alpha_EAS(8));
        my_mesh->AddElement(element);
        elemcount++;
    }

    // Deactivate automatic gravity in mesh
    my_mesh->SetAutomaticGravity(false);
    // Remember to add the mesh to the system!
    my_system.Add(my_mesh);

#ifndef CHRONO_MKL
    use_mkl = false;
#endif

    // Setup solver
    if (use_mkl) {
#ifdef CHRONO_MKL
        auto mkl_solver = chrono_types::make_shared<ChSolverMKL>();
        mkl_solver->LockSparsityPattern(true);
        mkl_solver->SetVerbose(true);
        my_system.SetSolver(mkl_solver);
#endif
    } else {
        auto solver = chrono_types::make_shared<ChSolverMINRES>();
        solver->SetMaxIterations(200);
        solver->EnableDiagonalPreconditioner(true);
        solver->SetVerbose(false);
        my_system.SetSolver(solver);
        my_system.SetSolverForceTolerance(1e-9);
    }

    // Setup integrator
    my_system.SetTimestepperType(ChTimestepper::Type::HHT);
    auto mystepper = std::static_pointer_cast<ChTimestepperHHT>(my_system.GetTimestepper());
    mystepper->SetAlpha(-0.01);
    mystepper->SetMaxiters(10000);
    mystepper->SetAbsTolerances(1e-09);
    mystepper->SetMode(ChTimestepperHHT::POSITION);
    mystepper->SetScaling(true);

    ChTimer<> timer;
    int num_iterations = 0;

    // Simulation loop
    for (unsigned int it = 0; it < num_steps; it++) {
        timer.start();
        my_system.DoStepDynamics(step_size);
        timer.stop();

        num_iterations += mystepper->GetNumIterations();
        std::cout << "time = " << my_system.GetChTime() << "\t" << nodetip->GetPos().z() << std::endl;
    }

    // Report run time and total number of iterations
    std::cout << "sim time: " << timer.GetTimeSeconds() << " Num iterations: " << num_iterations << std::endl;

    m_execTime = timer.GetTimeSeconds();
    addMetric("tip_y_position (mm)", 1000 * nodetip->GetPos().z());
    addMetric("avg_num_iterations", (double)num_iterations / num_steps);
    addMetric("avg_time_per_step (ms)", 1000 * m_execTime / num_steps);

    return true;
}

// ====================================================================================

int main(int argc, char* argv[]) {
    std::string out_dir = "../METRICS";
    if (!filesystem::create_directory(filesystem::path(out_dir))) {
        std::cout << "Error creating directory " << out_dir << std::endl;
        return 1;
    }

    BrickIso_GravTest test("metrics_FEA_EASBrickIso_Grav", "Chrono::FEA");
    test.setOutDir(out_dir);
    test.setVerbose(true);
    bool passed = test.run();
    test.print();

    return 0;
}

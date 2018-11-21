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
// Author: Alessandro Tasora
// =============================================================================
//
// ChronoFEA demo program for basic FEA functionality
// =============================================================================

#include "chrono/physics/ChSystemNSC.h"
#include "chrono/solver/ChSolverMINRES.h"

#include "chrono/fea/ChElementSpring.h"
#include "chrono/fea/ChElementTetra_4.h"
#include "chrono/fea/ChElementTetra_10.h"
#include "chrono/fea/ChElementHexa_8.h"
#include "chrono/fea/ChElementHexa_20.h"
#include "chrono/fea/ChMesh.h"
#include "chrono/fea/ChLinkPointFrame.h"

using namespace chrono;
using namespace chrono::fea;

int main(int argc, char* argv[]) {
    GetLog() << "Example: FEA techology for finite elements \n";
    GetLog() << "-------------------------------------------------\n";
    GetLog() << "LINEAR hexahedral element FEA  \n\n";

    // The physical system: it contains all physical objects.
    ChSystemNSC my_system;

    // Create a mesh, that is a container for groups
    // of elements and their referenced nodes.
    auto my_mesh = std::make_shared<ChMesh>();

    // Create a material, that must be assigned to each element,
    // and set its parameters
    auto mmaterial = std::make_shared<ChContinuumElastic>();
    mmaterial->Set_E(207e6);
    mmaterial->Set_v(0.3);

    // Create some nodes. These are the classical point-like
    // nodes with x,y,z degrees of freedom, that can be used
    // for many types of FEM elements in space.
    // While creating them, also set X0 undeformed positions.
    double sx = 0.01;
    double sy = 0.10;
    double sz = 0.01;
    auto mnode1 = std::make_shared<ChNodeFEAxyz>(ChVector<>(0, 0, 0));
    auto mnode2 = std::make_shared<ChNodeFEAxyz>(ChVector<>(0, 0, sz));
    auto mnode3 = std::make_shared<ChNodeFEAxyz>(ChVector<>(sx, 0, sz));
    auto mnode4 = std::make_shared<ChNodeFEAxyz>(ChVector<>(sx, 0, 0));
    auto mnode5 = std::make_shared<ChNodeFEAxyz>(ChVector<>(0, sy, 0));
    auto mnode6 = std::make_shared<ChNodeFEAxyz>(ChVector<>(0, sy, sz));
    auto mnode7 = std::make_shared<ChNodeFEAxyz>(ChVector<>(sx, sy, sz));
    auto mnode8 = std::make_shared<ChNodeFEAxyz>(ChVector<>(sx, sy, 0));

    // For example, set applied forces to nodes:
    mnode5->SetForce(ChVector<>(0, -1000, 0));
    mnode6->SetForce(ChVector<>(0, -1000, 0));
    mnode7->SetForce(ChVector<>(0, -1000, 0));
    mnode8->SetForce(ChVector<>(0, -1000, 0));

    // Remember to add nodes and elements to the mesh!
    my_mesh->AddNode(mnode1);
    my_mesh->AddNode(mnode2);
    my_mesh->AddNode(mnode3);
    my_mesh->AddNode(mnode4);
    my_mesh->AddNode(mnode5);
    my_mesh->AddNode(mnode6);
    my_mesh->AddNode(mnode7);
    my_mesh->AddNode(mnode8);

    // Create the tetrahedron element, and assign
    // it nodes and material
    auto melement1 = std::make_shared<ChElementHexa_8>();
    melement1->SetNodes(mnode1, mnode2, mnode3, mnode4, mnode5, mnode6, mnode7, mnode8);
    melement1->SetMaterial(mmaterial);

    // Remember to add elements to the mesh!
    my_mesh->AddElement(melement1);

    // Remember to add the mesh to the system!
    my_system.Add(my_mesh);

    // Create also a truss
    auto truss = std::make_shared<ChBody>();
    my_system.Add(truss);
    truss->SetBodyFixed(true);

    // Create a constraint between a node and the truss
    auto constraint1 = std::make_shared<ChLinkPointFrame>();
    constraint1->Initialize(my_mesh,  // node container
                            0,        // index of node in node container
                            truss);   // body to be connected to

    auto constraint2 = std::make_shared<ChLinkPointFrame>();
    constraint2->Initialize(my_mesh,  // node container
                            1,        // index of node in node container
                            truss);   // body to be connected to

    auto constraint3 = std::make_shared<ChLinkPointFrame>();
    constraint3->Initialize(my_mesh,  // node container
                            2,        // index of node in node container
                            truss);   // body to be connected to

    auto constraint4 = std::make_shared<ChLinkPointFrame>();
    constraint4->Initialize(my_mesh,  // node container
                            3,        // index of node in node container
                            truss);   // body to be connected to

    my_system.Add(constraint1);
    my_system.Add(constraint2);
    my_system.Add(constraint3);
    my_system.Add(constraint4);

    // Set no gravity
    // my_system.Set_G_acc(VNULL);

    // Mark system completion
    my_system.SetupInitial();

    // Perform a linear static analysis
    my_system.SetSolverType(ChSolver::Type::MINRES);
    auto msolver = std::static_pointer_cast<ChSolverMINRES>(my_system.GetSolver());
    msolver->SetDiagonalPreconditioning(true);
    msolver->SetVerbose(true);
    my_system.SetMaxItersSolverSpeed(100);
    my_system.SetTolForce(1e-12);

    my_system.DoStaticLinear();

    // Output result
    // GetLog()<<melement1.GetStiffnessMatrix()<<"\n";
    // GetLog()<<melement1.GetMatrB()<<"\n";
    GetLog() << mnode1->GetPos() << "\n";
    GetLog() << mnode2->GetPos() << "\n";
    GetLog() << mnode3->GetPos() << "\n";
    GetLog() << mnode4->GetPos() << "\n";
    GetLog() << "node5 displ: " << mnode5->GetPos() - mnode5->GetX0() << "\n";
    GetLog() << "node6 displ: " << mnode6->GetPos() - mnode6->GetX0() << "\n";
    GetLog() << "node7 displ: " << mnode7->GetPos() - mnode7->GetX0() << "\n";
    GetLog() << "node8 displ: " << mnode8->GetPos() - mnode8->GetX0() << "\n";

    return 0;
}

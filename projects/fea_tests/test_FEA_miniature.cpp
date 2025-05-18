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
// Authors: Alessandro Tasora
// =============================================================================
//
// FEA for 3D beams and constrains using the Matlab engine
//
// =============================================================================

#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChLinkMate.h"
#include "chrono/physics/ChLinkLock.h"
#include "chrono/physics/ChLinkMate.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/solver/ChIterativeSolverLS.h"

#include "chrono/fea/ChElementBeamEuler.h"
#include "chrono/fea/ChBuilderBeam.h"
#include "chrono/fea/ChMesh.h"

#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"

#include "chrono_matlab/ChMatlabEngine.h"
#include "chrono_matlab/ChSolverMatlab.h"

using namespace chrono;
using namespace chrono::fea;
using namespace chrono::irrlicht;

int main(int argc, char* argv[]) {
    // Set path to Chrono data directory
    SetChronoDataPath(CHRONO_DATA_DIR);

    // Create a Chrono::Engine physical system
    ChSystemNSC sys;

    double scales = 100;


    double thickZ = scales * 0.00015;
    double hbarW = scales * 0.00070;
    double hbarL1 = scales * 0.00381;
    double hbarL2 = scales * 0.00387;
    double hbarL3 = scales * 0.00381;
    double vbarL = scales * 0.01137;
    double vbarW = scales * 0.00006;
    double Rpinion = scales * 0.00040;
    double OffPin = scales * 0.00050;
    double Rbalance = scales * 0.00500;
    double Wbalance = scales * 0.00015;
    bool simple_rack = false;

    ChVector3d vAh(-hbarL1 - hbarL2 * 0.5, vbarL, 0);
    ChVector3d vBh(-hbarL2 * 0.5, vbarL, 0);
    ChVector3d vCh(hbarL2 * 0.5, vbarL, 0);
    ChVector3d vDh(hbarL1 + hbarL2 * 0.5, vbarL, 0);
    ChVector3d vAl(-hbarL1 - hbarL2 * 0.5, 0, 0);
    ChVector3d vBl(-hbarL2 * 0.5, 0, 0);
    ChVector3d vCl(hbarL2 * 0.5, 0, 0);
    ChVector3d vDl(hbarL1 + hbarL2 * 0.5, 0, 0);
    ChVector3d vP(0, -Rpinion - hbarW * 0.5, 0);

    // Create a truss:
    auto body_truss = chrono_types::make_shared<ChBody>();

    body_truss->SetFixed(true);

    sys.AddBody(body_truss);

    /*
    // Attach a 'box' shape asset for visualization.
    auto mboxtruss = chrono_types::make_shared<ChVisualShapeBox>();
    mboxtruss->GetBoxGeometry().Pos  = ChVector3d(-0.01, 0,0);
    mboxtruss->GetBoxGeometry().SetLengths( ChVector3d(0.02, 0.2, 0.1) );
    body_truss->AddAsset(mboxtruss);
    */

    // Create a FEM mesh, that is a container for groups
    // of elements and their referenced nodes.
    auto my_mesh = chrono_types::make_shared<ChMesh>();

    // Create the horizontal beams

    auto msectionH = chrono_types::make_shared<ChBeamSectionAdvanced>();

    msectionH->SetDensity(7000);  //***TEST*** must be 7k
    msectionH->SetYoungModulus(200.0e9);
    msectionH->SetShearModulusFromPoisson(0.32);
    msectionH->SetAsRectangularSection(hbarW, thickZ);
    msectionH->SetRayleighDamping(0.00);

    ChBuilderBeamEuler builder;

    builder.BuildBeam(my_mesh,               // the mesh where to put the created nodes and elements
                      msectionH,             // the ChBeamSectionAdvanced to use for the ChElementBeamEuler elements
                      2,                     // the number of ChElementBeamEuler to create
                      vAh,                   // the 'Ah' point in space (beginning of beam)
                      vBh,                   // the 'Bh' point in space (end of beam)
                      ChVector3d(0, 1, 0));  // the 'Y' up direction of the section for the beam

    // After having used BuildBeam(), you can retrieve the nodes used for the beam
    std::shared_ptr<ChNodeFEAxyzrot> node_Ah = builder.GetLastBeamNodes().front();
    std::shared_ptr<ChNodeFEAxyzrot> node_Bh = builder.GetLastBeamNodes().back();

    builder.BuildBeam(my_mesh,               // the mesh where to put the created nodes and elements
                      msectionH,             // the ChBeamSectionAdvanced to use for the ChElementBeamEuler elements
                      2,                     // the number of ChElementBeamEuler to create
                      node_Bh,               // the 'Bh' point in space (beginning of beam)
                      vCh,                   // the 'Ch' point in space (end of beam)
                      ChVector3d(0, 1, 0));  // the 'Y' up direction of the section for the beam

    // After having used BuildBeam(), you can retrieve the nodes used for the beam
    std::shared_ptr<ChNodeFEAxyzrot> node_Ch = builder.GetLastBeamNodes().back();

    builder.BuildBeam(my_mesh,               // the mesh where to put the created nodes and elements
                      msectionH,             // the ChBeamSectionAdvanced to use for the ChElementBeamEuler elements
                      2,                     // the number of ChElementBeamEuler to create
                      node_Ch,               // the 'Ch' point in space (beginning of beam)
                      vDh,                   // the 'Dh' point in space (end of beam)
                      ChVector3d(0, 1, 0));  // the 'Y' up direction of the section for the beam

    // After having used BuildBeam(), you can retrieve the nodes used for the beam
    std::shared_ptr<ChNodeFEAxyzrot> node_Dh = builder.GetLastBeamNodes().back();

    // Create the vertical flexible beams

    auto msectionV = chrono_types::make_shared<ChBeamSectionAdvanced>();

    msectionV->SetDensity(7000);  //***TEST*** must be 7k
    msectionV->SetYoungModulus(200.0e9);
    msectionV->SetShearModulusFromPoisson(0.32);
    msectionV->SetAsRectangularSection(vbarW, thickZ);
    msectionV->SetRayleighDamping(0.00);

    builder.BuildBeam(my_mesh,               // the mesh where to put the created nodes and elements
                      msectionV,             // the ChBeamSectionAdvanced to use for the ChElementBeamEuler elements
                      6,                     // the number of ChElementBeamEuler to create
                      node_Ah,               // the 'Ah' point in space (beginning of beam)
                      vAl,                   // the 'Al' point in space (end of beam)
                      ChVector3d(1, 0, 0));  // the 'Y' up direction of the section for the beam

    // After having used BuildBeam(), you can retrieve the nodes used for the beam
    std::shared_ptr<ChNodeFEAxyzrot> node_Al = builder.GetLastBeamNodes().back();

    node_Al->SetFixed(true);

    builder.BuildBeam(my_mesh,               // the mesh where to put the created nodes and elements
                      msectionV,             // the ChBeamSectionAdvanced to use for the ChElementBeamEuler elements
                      6,                     // the number of ChElementBeamEuler to create
                      node_Dh,               // the 'Dh' point in space (beginning of beam)
                      vDl,                   // the 'Dl' point in space (end of beam)
                      ChVector3d(1, 0, 0));  // the 'Y' up direction of the section for the beam

    // After having used BuildBeam(), you can retrieve the nodes used for the beam
    std::shared_ptr<ChNodeFEAxyzrot> node_Dl = builder.GetLastBeamNodes().back();

    node_Dl->SetFixed(true);

    // Create the inner vertical flexible beams

    builder.BuildBeam(my_mesh,               // the mesh where to put the created nodes and elements
                      msectionV,             // the ChBeamSectionAdvanced to use for the ChElementBeamEuler elements
                      6,                     // the number of ChElementBeamEuler to create
                      node_Bh,               // the 'Bh' point in space (beginning of beam)
                      vBl,                   // the 'Bl' point in space (end of beam)
                      ChVector3d(1, 0, 0));  // the 'Y' up direction of the section for the beam

    // After having used BuildBeam(), you can retrieve the nodes used for the beam
    std::shared_ptr<ChNodeFEAxyzrot> node_Bl = builder.GetLastBeamNodes().back();

    builder.BuildBeam(my_mesh,               // the mesh where to put the created nodes and elements
                      msectionV,             // the ChBeamSectionAdvanced to use for the ChElementBeamEuler elements
                      6,                     // the number of ChElementBeamEuler to create
                      node_Ch,               // the 'Dh' point in space (beginning of beam)
                      vCl,                   // the 'Dl' point in space (end of beam)
                      ChVector3d(1, 0, 0));  // the 'Y' up direction of the section for the beam

    // After having used BuildBeam(), you can retrieve the nodes used for the beam
    std::shared_ptr<ChNodeFEAxyzrot> node_Cl = builder.GetLastBeamNodes().back();

    // Create the rack

    if (simple_rack) {
        builder.BuildBeam(my_mesh,               // the mesh where to put the created nodes and elements
                          msectionH,             // the ChBeamSectionAdvanced to use for the ChElementBeamEuler elements
                          1,                     // the number of ChElementBeamEuler to create
                          node_Bl,               // the 'Cl' point in space (beginning of beam)
                          node_Cl,               // the 'Dl' point in space (end of beam)
                          ChVector3d(0, 1, 0));  // the 'Y' up direction of the section for the beam
    }

    //
    // Final touches to mesh..
    //

    // Remember to add the mesh to the system!
    sys.Add(my_mesh);

    // ==Asset== attach a visualization of the FEM mesh.
    // This will automatically update a triangle mesh (a ChVisualShapeTriangleMesh
    // asset that is internally managed) by setting  proper
    // coordinates and vertex colours as in the FEM elements.
    // Such triangle mesh can be rendered by Irrlicht or POVray or whatever
    // postprocessor that can handle a coloured ChVisualShapeTriangleMesh).
    // Do not forget AddAsset() at the end!

    auto mvisualizebeamA = chrono_types::make_shared<ChVisualShapeFEA>();
    mvisualizebeamA->SetFEMdataType(ChVisualShapeFEA::DataType::NODE_SPEED_NORM);
    mvisualizebeamA->SetColormapRange(-30, 30);
    mvisualizebeamA->SetSmoothFaces(true);
    mvisualizebeamA->SetWireframe(false);
    my_mesh->AddVisualShapeFEA(mvisualizebeamA);

    auto mvisualizebeamC = chrono_types::make_shared<ChVisualShapeFEA>();
    mvisualizebeamC->SetFEMglyphType(ChVisualShapeFEA::GlyphType::NODE_CSYS);
    mvisualizebeamC->SetFEMdataType(ChVisualShapeFEA::DataType::NONE);
    mvisualizebeamC->SetSymbolsThickness(0.001);
    mvisualizebeamC->SetSymbolsScale(0.01);
    mvisualizebeamC->SetZbufferHide(false);
    my_mesh->AddVisualShapeFEA(mvisualizebeamC);

    //
    // The balance and the rigid rach
    //

    if (!simple_rack) {
        auto rack = chrono_types::make_shared<ChBodyEasyBox>(hbarL2, hbarW, thickZ, 7000, true, false);
        rack->SetPos(0.5 * (vBl + vCl));
        sys.Add(rack);

        auto constr_B = chrono_types::make_shared<ChLinkMateGeneric>();
        constr_B->Initialize(node_Bl, rack, false, node_Bl->Frame(), node_Bl->Frame());
        sys.Add(constr_B);

        auto constr_C = chrono_types::make_shared<ChLinkMateGeneric>();
        constr_C->Initialize(node_Cl, rack, false, node_Cl->Frame(), node_Cl->Frame());
        sys.Add(constr_C);

        auto balance = chrono_types::make_shared<ChBodyEasyCylinder>(ChAxis::Y, Rbalance, Wbalance, 7000, true, false);
        balance->SetPos(vP + ChVector3d(0, 0, -OffPin));
        balance->SetRot(QuatFromAngleAxis(CH_PI_2, VECT_X));
        for (int i = 0; i < 6; ++i) {
            double phi = CH_2PI * (i / 6.0);
            auto p1 = ChVector3d(sin(phi) * Rbalance * 0.8, Wbalance, cos(phi) * Rbalance * 0.8);
            auto p2 = p1 + ChVector3d(0, 2 * Wbalance, 0);
            ChLineSegment seg(p1, p2);
            auto vshape = chrono_types::make_shared<ChVisualShapeCylinder>(Rbalance * 0.1, seg.GetLength());
            balance->AddVisualShape(vshape, seg.GetFrame());
        }
        ChLineSegment seg(vP + ChVector3d(0, -OffPin * 10, 0), vP + ChVector3d(0, OffPin * 10, 0));
        auto vshaft = chrono_types::make_shared<ChVisualShapeCylinder>(Rpinion, seg.GetLength());
        vshaft->SetColor(ChColor(0.5f, 0.9f, 0.9f));
        balance->AddVisualShape(vshaft, seg.GetFrame());

        sys.Add(balance);

        auto revolute = chrono_types::make_shared<ChLinkLockRevolute>();
        std::shared_ptr<ChBody> mbalance = balance;
        revolute->Initialize(mbalance, body_truss, ChFrame<>(vP + ChVector3d(0, 0, -0.01)));

        sys.Add(revolute);

        auto constr_rack = chrono_types::make_shared<ChLinkMateRackPinion>();
        constr_rack->Initialize(balance, rack, false, ChFrame<>(), ChFrame<>());

        ChFrameMoving<> f_pin_abs(vP);
        ChFrameMoving<> f_rack_abs(vP + ChVector3d(0, 0.1, 0));
        ChFrameMoving<> f_pin = balance->TransformParentToLocal(f_pin_abs);
        ChFrameMoving<> f_rack = rack->TransformParentToLocal(f_rack_abs);
        constr_rack->SetPinionRadius(Rpinion);
        constr_rack->SetPinionFrame(f_pin);
        constr_rack->SetRackFrame(f_rack);

        sys.Add(constr_rack);

        balance->SetAngVelParent(ChVector3d(0, 0, 1.5));
    }
    
    // Create the Irrlicht visualization system
    auto vis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
    vis->SetWindowSize(800, 600);
    vis->SetWindowTitle("Beams and constraints");
    vis->Initialize();
    vis->AddLogo();
    vis->AddSkyBox();
    vis->AddCamera(ChVector3d(0, 0.01*scales, 0.01*scales));
    vis->AddTypicalLights();
    vis->AttachSystem(&sys);

    //
    // THE SOFT-REAL-TIME CYCLE
    //

    auto solver = chrono_types::make_shared<ChSolverMINRES>();
    solver->EnableDiagonalPreconditioner(false);
    solver->EnableWarmStart(true);
    solver->SetMaxIterations(400);
    solver->SetTolerance(1e-12);
    solver->SetVerbose(true);
    sys.SetSolver(solver);

    //***TEST***
    ChMatlabEngine matlab_engine;
    auto matlab_solver = chrono_types::make_shared<ChSolverMatlab>(matlab_engine);
    sys.SetSolver(matlab_solver);

    sys.SetGravitationalAcceleration(ChVector3d(0, 0, 0));

    std::cout << "STATIC linear solve ----\n";
    node_Cl->SetForce(ChVector3d(50, 0, 0));
    // application.GetSystem()->DoStaticLinear();
    node_Cl->SetForce(ChVector3d(0, 0, 0));

    if (simple_rack) {
        node_Cl->SetForce(ChVector3d(50, 0, 0));
        sys.DoStaticNonlinear(12);
        node_Cl->SetForce(ChVector3d(0, 0, 0));
    }

    while (vis->Run()) {
        vis->BeginScene();
        vis->Render();
        tools::drawGrid(vis.get(), 0.2, 0.2, 10, 10, ChCoordsys<>(VNULL, CH_PI_2, VECT_Z), ChColor(0.4f, 0.5f, 0.5f),
                        true);
        vis->EndScene();
        sys.DoStepDynamics(0.01);
    }

    return 0;
}

// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
//   Show how to use the OpenCASCADE features
//   implemented in the unit_CASCADE:
//
//   - load a 3D model saved in STEP format from a CAD
//   - select some sub assemblies from the STEP model
//   - make Chrono::Engine objects out of those parts
// =============================================================================

#include "chrono/core/ChRealtimeStep.h"
#include "chrono/core/ChRandom.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono_cascade/ChCascadeBodyEasy.h"
#include "chrono_cascade/ChCascadeDoc.h"
#include "chrono_cascade/ChVisualShapeCascade.h"
#include "chrono/solver/ChSolverADMM.h"

#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"

// Use the namespace of Chrono
using namespace chrono;
using namespace chrono::irrlicht;

// Use the namespace with OpenCascade stuff
using namespace cascade;

int main(int argc, char* argv[]) {
    // Set path to Chrono data directories
    SetChronoDataPath(CHRONO_DATA_DIR);

    // 1- Create a ChronoENGINE physical system: all bodies and constraints
    //    will be handled by this ChSystemNSC object.
    ChSystemNSC sys;

    // Create a surface material to be used for collisions, if any
    auto mysurfmaterial = chrono_types::make_shared<ChContactMaterialNSC>();
    mysurfmaterial->SetFriction(0.3f);
    mysurfmaterial->SetRestitution(0);

    //
    // Load a STEP file, containing a mechanism. The demo STEP file has been
    // created using a 3D CAD (in this case, SolidEdge v.18).
    //

    // Create the ChCascadeDoc, a container that loads the STEP model
    // and manages its subassemblies
    ChCascadeDoc mydoc;

    // load the STEP model using this command:
    bool load_ok = mydoc.Load_STEP(GetChronoDataFile("/cascade/IRB7600_23_500_m2000_rev1_01_decorated.stp").c_str());

    // print the contained shapes
    mydoc.Dump(std::cout);

    ChCollisionModel::SetDefaultSuggestedEnvelope(0.002);
    ChCollisionModel::SetDefaultSuggestedMargin(0.001);

    //
    // Retrieve some sub shapes from the loaded model, using
    // the GetNamedShape() function, that can use path/subpath/subsubpath/part
    // syntax and * or ? wldcards, etc.
    //

    std::shared_ptr<ChCascadeBodyEasy> body_base;
    std::shared_ptr<ChCascadeBodyEasy> body_turret;
    std::shared_ptr<ChCascadeBodyEasy> body_bicep;
    std::shared_ptr<ChCascadeBodyEasy> body_elbow;
    std::shared_ptr<ChCascadeBodyEasy> body_forearm;
    std::shared_ptr<ChCascadeBodyEasy> body_wrist;
    std::shared_ptr<ChCascadeBodyEasy> body_hand;
    std::shared_ptr<ChCascadeBodyEasy> body_cylinder;
    std::shared_ptr<ChCascadeBodyEasy> body_rod;

    // Note, In most CADs the Y axis is horizontal, but we want it vertical.
    // So define a root transformation for rotating all the imported objects.
    ChQuaternion<> rotation1;
    rotation1.SetFromAngleAxis(-CH_C_PI / 2, VECT_X);  // 1: rotate 90� on X axis
    ChQuaternion<> rotation2;
    rotation2.SetFromAngleAxis(CH_C_PI, VECT_Y);            // 2: rotate 180� on vertical Y axis
    ChQuaternion<> tot_rotation = rotation2 % rotation1;  // rotate on 1 then on 2, using quaternion product
    ChFrameMoving<> root_frame(ChVector3d(0, 0, 0), tot_rotation);

    if (load_ok) {
        TopoDS_Shape shape_base;
        if (mydoc.GetNamedShape(shape_base, "Assem10/Assem8")) {
            body_base = chrono_types::make_shared<ChCascadeBodyEasy>(shape_base, 1000, true, false);
            sys.Add(body_base);

            // The base is fixed to the ground
            body_base->SetBodyFixed(true);
            body_base->ConcatenatePreTransformation(root_frame);
        } else
            std::cout << "Warning. Desired object not found in document \n";

        TopoDS_Shape shape_turret;
        if (mydoc.GetNamedShape(shape_turret, "Assem10/Assem4")) {
            auto mbody = chrono_types::make_shared<ChCascadeBodyEasy>(shape_turret, 1000, true, false);
            body_turret = mbody;
            sys.Add(body_turret);
            // Move the body as for global displacement/rotation
            body_turret->ConcatenatePreTransformation(root_frame);
        } else
            std::cout << "Warning. Desired object not found in document \n";

        TopoDS_Shape shape_bicep;
        if (mydoc.GetNamedShape(shape_bicep, "Assem10/Assem1")) {
            body_bicep = chrono_types::make_shared<ChCascadeBodyEasy>(shape_bicep, 1000, true, false);
            sys.Add(body_bicep);
            // Move the body as for global displacement/rotation
            body_bicep->ConcatenatePreTransformation(root_frame);
        } else
            std::cout << "Warning. Desired object not found in document \n";

        TopoDS_Shape shape_elbow;
        if (mydoc.GetNamedShape(shape_elbow, "Assem10/Assem5")) {
            body_elbow = chrono_types::make_shared<ChCascadeBodyEasy>(shape_elbow, 1000, true, false);
            sys.Add(body_elbow);
            // Move the body as for global displacement/rotation
            body_elbow->ConcatenatePreTransformation(root_frame);
        } else
            std::cout << "Warning. Desired object not found in document \n";

        TopoDS_Shape shape_forearm;
        if (mydoc.GetNamedShape(shape_forearm, "Assem10/Assem7")) {
            body_forearm = chrono_types::make_shared<ChCascadeBodyEasy>(shape_forearm, 1000, true, false);
            sys.Add(body_forearm);
            // Move the body as for global displacement/rotation
            body_forearm->ConcatenatePreTransformation(root_frame);
        } else
            std::cout << "Warning. Desired object not found in document \n";

        TopoDS_Shape shape_wrist;
        if (mydoc.GetNamedShape(shape_wrist, "Assem10/Assem6")) {
            body_wrist = chrono_types::make_shared<ChCascadeBodyEasy>(shape_wrist, 1000, true, false);
            sys.Add(body_wrist);
            // Move the body as for global displacement/rotation
            body_wrist->ConcatenatePreTransformation(root_frame);
        } else
            std::cout << "Warning. Desired object not found in document \n";

        TopoDS_Shape shape_hand;
        if (mydoc.GetNamedShape(shape_hand, "Assem10/Assem9")) {
            body_hand = chrono_types::make_shared<ChCascadeBodyEasy>(shape_hand, 1000, true, false);
            sys.Add(body_hand);
            // Move the body as for global displacement/rotation
            body_hand->ConcatenatePreTransformation(root_frame);

            // also create a collision shape attached to the hand:
            auto mcube = chrono_types::make_shared<ChBodyEasyBox>(0.2, 0.08, 0.08, 1000, true, true, mysurfmaterial);
            mcube->SetCsys(ChCoordsys<>(ChVector3d(0.1, 0, 0)) >> body_hand->GetCsys());
            sys.Add(mcube);
            auto mcubelink = chrono_types::make_shared<ChLinkLockLock>();
            mcubelink->Initialize(mcube, body_hand, mcube->GetCsys());
            sys.Add(mcubelink);

        } else
            std::cout << "Warning. Desired object not found in document \n";

        TopoDS_Shape shape_cylinder;
        if (mydoc.GetNamedShape(shape_cylinder, "Assem10/Assem3")) {
            body_cylinder = chrono_types::make_shared<ChCascadeBodyEasy>(shape_cylinder, 1000, true, false);
            sys.Add(body_cylinder);
            // Move the body as for global displacement/rotation
            body_cylinder->ConcatenatePreTransformation(root_frame);
        } else
            std::cout << "Warning. Desired object not found in document \n";

        TopoDS_Shape shape_rod;
        if (mydoc.GetNamedShape(shape_rod, "Assem10/Assem2")) {
            body_rod = chrono_types::make_shared<ChCascadeBodyEasy>(shape_rod, 1000, true, false);
            sys.Add(body_rod);
            // Move the body as for global displacement/rotation
            body_rod->ConcatenatePreTransformation(root_frame);
        } else
            std::cout << "Warning. Desired object not found in document \n";

    } else
        std::cout << "Warning. Desired STEP file could not be opened/parsed \n";

    if (!body_base || !body_turret || !body_bicep || !body_elbow || !body_forearm || !body_wrist || !body_hand) {
        return 1;
    }

    // Create joints between two parts.
    // To understand where is the axis of the joint, we can exploit the fact
    // that in the STEP file that we prepared for this demo, we inserted some
    // objects called 'marker' and we placed them aligned to the shafts, so now
    // we can fetch them and get their position/rotation.

    TopoDS_Shape shape_marker;

    ChFrame<> frame_marker_base_turret;
    if (mydoc.GetNamedShape(shape_marker, "Assem10/Assem8/marker#1"))
        ChCascadeDoc::FromCascadeToChrono(shape_marker.Location(), frame_marker_base_turret);
    else
        std::cout << "Warning. Desired marker not found in document \n";
    // Transform the abs coordinates of the marker because everything was rotated/moved by 'root_frame' :
    frame_marker_base_turret >>= root_frame;

    std::shared_ptr<ChLinkLockRevolute> my_link1(new ChLinkLockRevolute);
    my_link1->Initialize(body_base, body_turret, frame_marker_base_turret.GetCsys());
    sys.AddLink(my_link1);

    ChFrame<> frame_marker_turret_bicep;
    if (mydoc.GetNamedShape(shape_marker, "Assem10/Assem4/marker#2"))
        ChCascadeDoc::FromCascadeToChrono(shape_marker.Location(), frame_marker_turret_bicep);
    else
        std::cout << "Warning. Desired marker not found in document \n";
    frame_marker_turret_bicep >>= root_frame;

    std::shared_ptr<ChLinkLockRevolute> my_link2(new ChLinkLockRevolute);
    my_link2->Initialize(body_turret, body_bicep, frame_marker_turret_bicep.GetCsys());
    sys.AddLink(my_link2);

    ChFrame<> frame_marker_bicep_elbow;
    if (mydoc.GetNamedShape(shape_marker, "Assem10/Assem1/marker#2"))
        ChCascadeDoc::FromCascadeToChrono(shape_marker.Location(), frame_marker_bicep_elbow);
    else
        std::cout << "Warning. Desired marker not found in document \n";
    frame_marker_bicep_elbow >>= root_frame;

    std::shared_ptr<ChLinkLockRevolute> my_link3(new ChLinkLockRevolute);
    my_link3->Initialize(body_bicep, body_elbow, frame_marker_bicep_elbow.GetCsys());
    sys.AddLink(my_link3);

    ChFrame<> frame_marker_elbow_forearm;
    if (mydoc.GetNamedShape(shape_marker, "Assem10/Assem5/marker#2"))
        ChCascadeDoc::FromCascadeToChrono(shape_marker.Location(), frame_marker_elbow_forearm);
    else
        std::cout << "Warning. Desired marker not found in document \n";
    frame_marker_elbow_forearm >>= root_frame;

    std::shared_ptr<ChLinkLockRevolute> my_link4(new ChLinkLockRevolute);
    my_link4->Initialize(body_elbow, body_forearm, frame_marker_elbow_forearm.GetCsys());
    sys.AddLink(my_link4);

    ChFrame<> frame_marker_forearm_wrist;
    if (mydoc.GetNamedShape(shape_marker, "Assem10/Assem7/marker#2"))
        ChCascadeDoc::FromCascadeToChrono(shape_marker.Location(), frame_marker_forearm_wrist);
    else
        std::cout << "Warning. Desired marker not found in document \n";
    frame_marker_forearm_wrist >>= root_frame;

    std::shared_ptr<ChLinkLockRevolute> my_link5(new ChLinkLockRevolute);
    my_link5->Initialize(body_forearm, body_wrist, frame_marker_forearm_wrist.GetCsys());
    sys.AddLink(my_link5);

    ChFrame<> frame_marker_wrist_hand;
    if (mydoc.GetNamedShape(shape_marker, "Assem10/Assem6/marker#2"))
        ChCascadeDoc::FromCascadeToChrono(shape_marker.Location(), frame_marker_wrist_hand);
    else
        std::cout << "Warning. Desired marker not found in document \n";
    frame_marker_wrist_hand >>= root_frame;

    std::shared_ptr<ChLinkLockRevolute> my_link6(new ChLinkLockRevolute);
    my_link6->Initialize(body_wrist, body_hand, frame_marker_wrist_hand.GetCsys());
    sys.AddLink(my_link6);

    ChFrame<> frame_marker_turret_cylinder;
    if (mydoc.GetNamedShape(shape_marker, "Assem10/Assem4/marker#3"))
        ChCascadeDoc::FromCascadeToChrono(shape_marker.Location(), frame_marker_turret_cylinder);
    else
        std::cout << "Warning. Desired marker not found in document \n";
    frame_marker_turret_cylinder >>= root_frame;

    std::shared_ptr<ChLinkLockRevolute> my_link7(new ChLinkLockRevolute);
    my_link7->Initialize(body_turret, body_cylinder, frame_marker_turret_cylinder.GetCsys());
    sys.AddLink(my_link7);

    ChFrame<> frame_marker_cylinder_rod;
    if (mydoc.GetNamedShape(shape_marker, "Assem10/Assem3/marker#2"))
        ChCascadeDoc::FromCascadeToChrono(shape_marker.Location(), frame_marker_cylinder_rod);
    else
        std::cout << "Warning. Desired marker not found in document \n";
    frame_marker_cylinder_rod >>= root_frame;

    std::shared_ptr<ChLinkLockCylindrical> my_link8(new ChLinkLockCylindrical);
    my_link8->Initialize(body_cylinder, body_rod, frame_marker_cylinder_rod.GetCsys());
    sys.AddLink(my_link8);

    ChFrame<> frame_marker_rod_bicep;
    if (mydoc.GetNamedShape(shape_marker, "Assem10/Assem2/marker#2"))
        ChCascadeDoc::FromCascadeToChrono(shape_marker.Location(), frame_marker_rod_bicep);
    else
        std::cout << "Warning. Desired marker not found in document \n";
    frame_marker_rod_bicep >>= root_frame;

    std::shared_ptr<ChLinkLockCylindrical> my_link9(new ChLinkLockCylindrical);
    my_link9->Initialize(body_rod, body_bicep, frame_marker_rod_bicep.GetCsys());
    sys.AddLink(my_link9);

    // Add a couple of markers for the 'lock' constraint between the hand and the
    // absolute reference: when we will move the marker in absolute reference, the
    // hand will follow it.

    std::shared_ptr<ChMarker> my_marker_hand(new ChMarker);
    std::shared_ptr<ChMarker> my_marker_move(new ChMarker);

    body_hand->AddMarker(my_marker_hand);
    body_base->AddMarker(my_marker_move);

    ChQuaternion<> rot_on_x;
    rot_on_x.SetFromAngleAxis(CH_C_PI / 2, VECT_X);
    ChFrame<> frame_marker_move = ChFrame<>(VNULL, rot_on_x) >> frame_marker_wrist_hand;

    my_marker_hand->Impose_Abs_Coord(frame_marker_wrist_hand.GetCsys());
    my_marker_move->Impose_Abs_Coord(frame_marker_move.GetCsys());

    std::shared_ptr<ChLinkLockLock> my_link_teacher(new ChLinkLockLock);
    my_link_teacher->Initialize(my_marker_hand, my_marker_move);
    sys.AddLink(my_link_teacher);

    // Set motions for Z and Y coordinates of the 'my_link_teacher' marker,
    // so that the hand will follow it. To do so, we create four segments of
    // motion for Z coordinate and four for Y coordinate, we join them with
    // ChFunctionSequence and we repeat sequence by ChFunctionRepeat

    std::shared_ptr<ChFunctionConstAcc> motlaw_z1(new ChFunctionConstAcc);
    motlaw_z1->Set_h(-0.7);
    motlaw_z1->Set_end(1);
    std::shared_ptr<ChFunctionConst> motlaw_z2(new ChFunctionConst);
    std::shared_ptr<ChFunctionConstAcc> motlaw_z3(new ChFunctionConstAcc);
    motlaw_z3->Set_h(0.7);
    motlaw_z3->Set_end(1);
    std::shared_ptr<ChFunctionConst> motlaw_z4(new ChFunctionConst);
    std::shared_ptr<ChFunctionSequence> motlaw_z_seq(new ChFunctionSequence);
    motlaw_z_seq->InsertFunct(motlaw_z1, 1, 1, true);
    motlaw_z_seq->InsertFunct(motlaw_z2, 1, 1, true);  // true = force c0 continuity, traslating fx
    motlaw_z_seq->InsertFunct(motlaw_z3, 1, 1, true);
    motlaw_z_seq->InsertFunct(motlaw_z4, 1, 1, true);
    auto motlaw_z = chrono_types::make_shared<ChFunctionRepeat>(motlaw_z_seq);
    motlaw_z->Set_window_length(4);

    std::shared_ptr<ChFunctionConst> motlaw_y1(new ChFunctionConst);
    std::shared_ptr<ChFunctionConstAcc> motlaw_y2(new ChFunctionConstAcc);
    motlaw_y2->Set_h(-0.6);
    motlaw_y2->Set_end(1);
    std::shared_ptr<ChFunctionConst> motlaw_y3(new ChFunctionConst);
    std::shared_ptr<ChFunctionConstAcc> motlaw_y4(new ChFunctionConstAcc);
    motlaw_y4->Set_h(0.6);
    motlaw_y4->Set_end(1);
    std::shared_ptr<ChFunctionSequence> motlaw_y_seq(new ChFunctionSequence);
    motlaw_y_seq->InsertFunct(motlaw_y1, 1, 1, true);
    motlaw_y_seq->InsertFunct(motlaw_y2, 1, 1, true);  // true = force c0 continuity, traslating fx
    motlaw_y_seq->InsertFunct(motlaw_y3, 1, 1, true);
    motlaw_y_seq->InsertFunct(motlaw_y4, 1, 1, true);
    auto motlaw_y = chrono_types::make_shared<ChFunctionRepeat>(motlaw_y_seq);
    motlaw_y->Set_window_length(4);

    my_marker_move->SetMotion_Z(motlaw_z);
    my_marker_move->SetMotion_Y(motlaw_y);

    // Create a large cube as a floor.

    std::shared_ptr<ChBodyEasyBox> mfloor(new ChBodyEasyBox(8, 1, 8, 1000, true, true, mysurfmaterial));
    mfloor->SetBodyFixed(true);
    mfloor->SetPos(ChVector3d(0, -0.5, 0));
    mfloor->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/blue.png"));
    sys.Add(mfloor);

    // Create a stack of boxes to be impacted
    if (true) {
        double brick_h = 0.3;
        for (int ix = 0; ix < 3; ++ix)
            for (int ib = 0; ib < 6; ++ib) {
                std::shared_ptr<ChBodyEasyBox> cube(
                    new ChBodyEasyBox(0.4, brick_h, 0.4, 1000, true, true, mysurfmaterial));
                cube->SetPos(ChVector3d(-1.4, (0.5 * brick_h) + ib * brick_h, -0.4 - 0.5 * ix));
                cube->SetRot(QuatFromAngleAxis(ib * 0.1, VECT_Y));
                cube->GetVisualShape(0)->SetColor(ChColor(0.5f + float(0.5 * ChRandom::Get()),  //
                                                          0.5f + float(0.5 * ChRandom::Get()),  //
                                                          0.5f + float(0.5 * ChRandom::Get())   //
                                                          ));
                sys.Add(cube);
            }
    }

    // Create the Irrlicht visualization system
    auto vis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
    vis->SetWindowSize(800, 600);
    vis->SetWindowTitle("Load a robot model from STEP file");
    vis->Initialize();
    vis->AddLogo();
    vis->AddSkyBox();
    vis->AddCamera(ChVector3d(0.2, 1.6, 3.5), ChVector3d(0, 1, 0));
    vis->AddTypicalLights();
    vis->AttachSystem(&sys);

    // Modify the settings of the solver.
    // By default, the solver might not have sufficient precision to keep the
    // robot joints 'mounted'. Expecially, the SOR, SSOR and other fixed point methods
    // cannot simulate well this robot problem because the mass of the last body in the
    // kinematic chain, i.e. the hand, is very low when compared to other bodies, so
    // the convergence of the solver would be bad when 'pulling the hand' as in this
    // 'teaching mode' IK.
    // So switch to a more precise solver, ex. BARZILAIBORWEIN

    sys.SetSolverType(ChSolver::Type::BARZILAIBORWEIN);
    sys.SetSolverMaxIterations(200);

    /*
    // Alternative: the ADMM solver offers higher precision and it can also support FEA + nonsmooth contacts
    auto solver = chrono_types::make_shared<ChSolverADMM>(); //faster, if MKL enabled:
    chrono_types::make_shared<ChSolverPardisoMKL>()); solver->EnableWarmStart(true); solver->SetMaxIterations(60);
    solver->SetRho(1);
    solver->SetSigma(1e-8);
    solver->SetStepAdjustPolicy(ChSolverADMM::AdmmStepType::BALANCED_UNSCALED);
    sys.SetSolver(solver);
    */

    // Simulation loop
    ChRealtimeStepTimer realtime_timer;
    double time_step = 0.01;

    while (vis->Run()) {
        vis->BeginScene();
        vis->Render();
        tools::drawChFunction(vis.get(), motlaw_z, 0, 10, -0.9, 0.2, 10, 400, 300, 80);
        tools::drawChFunction(vis.get(), motlaw_y, 0, 10, -0.9, 0.2, 10, 500, 300, 80);
        vis->EndScene();

        sys.DoStepDynamics(time_step);
        realtime_timer.Spin(time_step);
    }

    return 0;
}

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
// Authors: Daniel Melanz
// =============================================================================
//
// ChronoParallel test program using NSC method for frictional contact.
//
// The global reference frame has Y up.
//
// If available, OpenGL is used for run-time rendering. Otherwise, the
// simulation is carried out for a pre-defined duration and output files are
// generated for post-processing with POV-Ray.
// =============================================================================

#include <cstdio>
#include <vector>
#include <cmath>

#include "chrono/ChConfig.h"
#include "chrono/physics/ChLinkMotorRotationAngle.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsInputOutput.h"
#include "chrono/utils/ChUtilsGenerators.h"

#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_parallel/solver/ChSystemDescriptorParallel.h"

// Control use of OpenGL run-time rendering
// Note: CHRONO_OPENGL is defined in ChConfig.h
//#undef CHRONO_OPENGL

#ifdef CHRONO_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif

using namespace chrono;
using namespace chrono::collision;

const char* out_folder = "../DEMO_SUSPENSION/POVRAY";

// =============================================================================
// Generate postprocessing output with current system state.
// =============================================================================
void OutputData(ChSystemParallel* sys, int out_frame, double time) {
    char filename[100];
    sprintf(filename, "%s/data_%03d.dat", out_folder, out_frame);
    utils::WriteShapesPovray(sys, filename);
    std::cout << "time = " << time << std::flush << std::endl;
}

// =============================================================================
// First of all, define a class for the 'car' (that is, a set of
// bodies and links which are grouped within this class; so it is
// easier to manage data structures in this example).
// =============================================================================

class MySimpleCar {
  public:
    // THE DATA

    double throttle;          // actual value 0...1 of gas throttle.
    double conic_tau;         // the transmission ratio of the conic gears at the rear axle
    double gear_tau;          // the actual tau of the gear
    double max_motor_torque;  // the max torque of the motor [Nm];
    double max_motor_speed;   // the max rotation speed of the motor [rads/s]

    // The parts making the car, as 3d Irrlicht scene nodes, each containing
    // the ChBody object
    // .. truss:
    std::shared_ptr<ChBody> truss;
    // .. right front suspension:
    std::shared_ptr<ChBody> spindleRF;
    std::shared_ptr<ChBody> wheelRF;
    std::shared_ptr<ChLinkLockRevolute> link_revoluteRF;
    std::shared_ptr<ChLinkDistance> link_distRFU1;
    std::shared_ptr<ChLinkDistance> link_distRFU2;
    std::shared_ptr<ChLinkDistance> link_distRFL1;
    std::shared_ptr<ChLinkDistance> link_distRFL2;
    std::shared_ptr<ChLinkTSDA> link_springRF;
    std::shared_ptr<ChLinkDistance> link_distRSTEER;
    // .. left front suspension:
    std::shared_ptr<ChBody> spindleLF;
    std::shared_ptr<ChBody> wheelLF;
    std::shared_ptr<ChLinkLockRevolute> link_revoluteLF;
    std::shared_ptr<ChLinkDistance> link_distLFU1;
    std::shared_ptr<ChLinkDistance> link_distLFU2;
    std::shared_ptr<ChLinkDistance> link_distLFL1;
    std::shared_ptr<ChLinkDistance> link_distLFL2;
    std::shared_ptr<ChLinkTSDA> link_springLF;
    std::shared_ptr<ChLinkDistance> link_distLSTEER;
    // .. right back suspension:
    std::shared_ptr<ChBody> spindleRB;
    std::shared_ptr<ChBody> wheelRB;
    std::shared_ptr<ChLinkLockRevolute> link_revoluteRB;
    std::shared_ptr<ChLinkDistance> link_distRBU1;
    std::shared_ptr<ChLinkDistance> link_distRBU2;
    std::shared_ptr<ChLinkDistance> link_distRBL1;
    std::shared_ptr<ChLinkDistance> link_distRBL2;
    std::shared_ptr<ChLinkTSDA> link_springRB;
    std::shared_ptr<ChLinkDistance> link_distRBlat;
    std::shared_ptr<ChLinkMotorRotationAngle> link_engineL;
    // .. left back suspension:
    std::shared_ptr<ChBody> spindleLB;
    std::shared_ptr<ChBody> wheelLB;
    std::shared_ptr<ChLinkLockRevolute> link_revoluteLB;
    std::shared_ptr<ChLinkDistance> link_distLBU1;
    std::shared_ptr<ChLinkDistance> link_distLBU2;
    std::shared_ptr<ChLinkDistance> link_distLBL1;
    std::shared_ptr<ChLinkDistance> link_distLBL2;
    std::shared_ptr<ChLinkTSDA> link_springLB;
    std::shared_ptr<ChLinkDistance> link_distLBlat;
    std::shared_ptr<ChLinkMotorRotationAngle> link_engineR;
    // THE FUNCTIONS

    // Build and initialize the car, creating all bodies corresponding to
    // the various parts and adding them to the physical system - also creating
    // and adding constraints to the system.
    MySimpleCar(ChSystemParallelNSC* my_system) {
        throttle = 0;  // initially, gas throttle is 0.
        conic_tau = 0.2;
        gear_tau = 0.3;
        max_motor_torque = 800;
        max_motor_speed = 8000;
        double frontDamping = .1;
        double rearDamping = .1;
        bool useSpheres = true;
        double angularSpeed = 20;

        // Create the wheel material
        auto mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
        mat->SetFriction(1.0);

        // --- The car body ---

        truss = chrono_types::make_shared<ChBody>(chrono_types::make_shared<ChCollisionModelParallel>());
        truss->SetIdentifier(-1);
        truss->SetMass(2086.524902);
        truss->SetPos(ChVector<>(0, 0.52349, 0.055765));
        truss->SetRot(ChQuaternion<>(1, 0, 0, 0));
        truss->SetInertiaXX(ChVector<>(3570.20377, 1078.52344, 2955.66050));
        truss->GetCollisionModel()->ClearModel();
        utils::AddBoxGeometry(truss.get(), mat, ChVector<>(0.7, 0.5, 3) * 0.5, ChVector<>(0, 0, 0));
        truss->GetCollisionModel()->BuildModel();
        truss->SetCollide(false);
        truss->SetBodyFixed(false);
        my_system->AddBody(truss);

        // --- Right Front suspension ---

        // ..the car right-front spindle
        spindleRF = chrono_types::make_shared<ChBody>(chrono_types::make_shared<ChCollisionModelParallel>());
        spindleRF->SetIdentifier(-2);
        spindleRF->SetMass(14.705);
        spindleRF->SetPos(ChVector<>(0.751, -0.026, 1.648965));
        spindleRF->SetRot(ChQuaternion<>(1, 0, 0, 0));
        spindleRF->SetInertiaXX(ChVector<>(0.07352, 0.04117, 0.04117));
        spindleRF->GetCollisionModel()->ClearModel();
        utils::AddBoxGeometry(spindleRF.get(), mat, ChVector<>(0.1, 0.4, 0.4) * 0.5, ChVector<>(0, 0, 0));
        spindleRF->GetCollisionModel()->BuildModel();
        spindleRF->SetCollide(false);
        my_system->AddBody(spindleRF);

        // ..the car right-front wheel
        wheelRF = chrono_types::make_shared<ChBody>(chrono_types::make_shared<ChCollisionModelParallel>());
        wheelRF->SetIdentifier(-3);
        wheelRF->SetMass(3.0);
        wheelRF->SetPos(ChVector<>(0.91, -0.026, 1.648965));
        wheelRF->SetRot(ChQuaternion<>(1, 0, 0, 0));
        wheelRF->SetInertiaXX(ChVector<>(0.2, 0.2, 0.2));
        wheelRF->GetCollisionModel()->ClearModel();
        if (useSpheres) {
            utils::AddSphereGeometry(wheelRF.get(), mat, 0.45, ChVector<>(0, 0, 0));
        } else {
            utils::AddCylinderGeometry(wheelRF.get(), mat, 0.45, 0.2, ChVector<>(0, 0, 0), Q_from_AngZ(CH_C_PI_2));
        }
        wheelRF->GetCollisionModel()->BuildModel();
        wheelRF->SetCollide(true);
        my_system->AddBody(wheelRF);

        // .. create the revolute joint between the wheel and the spindle
        link_revoluteRF = chrono_types::make_shared<ChLinkLockRevolute>();  // right, front, upper, 1
        link_revoluteRF->Initialize(wheelRF, spindleRF, ChCoordsys<>(ChVector<>(0.91, -0.026, 1.648965),
                                                                     chrono::Q_from_AngAxis(CH_C_PI / 2, VECT_Y)));
        my_system->AddLink(link_revoluteRF);

        // .. impose distance between two parts (as a massless rod with two spherical joints at the end)
        link_distRFU1 = chrono_types::make_shared<ChLinkDistance>();  // right, front, upper, 1
        link_distRFU1->Initialize(truss, spindleRF, false, ChVector<>(0.446, 0.245, 1.640965),
                                  ChVector<>(0.716, 0.215, 1.635965));
        my_system->AddLink(link_distRFU1);

        link_distRFU2 = chrono_types::make_shared<ChLinkDistance>();  // right, front, upper, 2
        link_distRFU2->Initialize(truss, spindleRF, false, ChVector<>(0.478, 0.196, 1.420965),
                                  ChVector<>(0.716, 0.215, 1.635965));
        my_system->AddLink(link_distRFU2);

        link_distRFL1 = chrono_types::make_shared<ChLinkDistance>();  // right, front, lower, 1
        link_distRFL1->Initialize(truss, spindleRF, false, ChVector<>(0.307, 0, 1.911965),
                                  ChVector<>(0.787, -0.118, 1.652965));
        my_system->AddLink(link_distRFL1);

        link_distRFL2 = chrono_types::make_shared<ChLinkDistance>();  // right, front, lower, 2
        link_distRFL2->Initialize(truss, spindleRF, false, ChVector<>(0.307, 0, 1.465965),
                                  ChVector<>(0.787, -0.118, 1.652965));
        my_system->AddLink(link_distRFL2);

        // .. create the spring between the truss and the spindle
        link_springRF = chrono_types::make_shared<ChLinkTSDA>();
        link_springRF->Initialize(truss, spindleRF, false, ChVector<>(0.498, 0.323, 1.792965),
                                  ChVector<>(0.543, -0.047, 1.785965));
        link_springRF->SetSpringCoefficient(167062.000);
        link_springRF->SetDampingCoefficient(frontDamping);
        link_springRF->SetRestLength(0.339);
        my_system->AddLink(link_springRF);

        // .. create the rod for steering the wheel
        link_distRSTEER = chrono_types::make_shared<ChLinkDistance>();  // right steer
        link_distRSTEER->Initialize(truss, spindleRF, false, ChVector<>(0.448, 0.054, 1.438965),
                                    ChVector<>(0.821, -0.016, 1.512965));
        my_system->AddLink(link_distRSTEER);

        // --- Left Front suspension ---

        // ..the car left-front spindle
        spindleLF = chrono_types::make_shared<ChBody>(chrono_types::make_shared<ChCollisionModelParallel>());
        spindleLF->SetIdentifier(-4);
        spindleLF->SetMass(14.705);
        spindleLF->SetPos(ChVector<>(-0.751, -0.026, 1.648965));
        spindleLF->SetRot(ChQuaternion<>(1, 0, 0, 0));
        spindleLF->SetInertiaXX(ChVector<>(0.07352, 0.04117, 0.04117));
        spindleLF->GetCollisionModel()->ClearModel();
        utils::AddBoxGeometry(spindleLF.get(), mat, ChVector<>(0.1, 0.4, 0.4) * 0.5, ChVector<>(0, 0, 0));
        spindleLF->GetCollisionModel()->BuildModel();
        spindleLF->SetCollide(false);
        my_system->AddBody(spindleLF);

        // ..the car left-front wheel
        wheelLF = chrono_types::make_shared<ChBody>(chrono_types::make_shared<ChCollisionModelParallel>());
        wheelLF->SetIdentifier(-5);
        wheelLF->SetMass(3.0);
        wheelLF->SetPos(ChVector<>(-0.91, -0.026, 1.648965));
        wheelLF->SetRot(ChQuaternion<>(1, 0, 0, 0));
        wheelLF->SetInertiaXX(ChVector<>(0.2, 0.2, 0.2));
        wheelLF->GetCollisionModel()->ClearModel();
        if (useSpheres) {
            utils::AddSphereGeometry(wheelLF.get(), mat, 0.45, ChVector<>(0, 0, 0));
        } else {
            utils::AddCylinderGeometry(wheelLF.get(), mat, 0.45, 0.2, ChVector<>(0, 0, 0), Q_from_AngZ(CH_C_PI_2));
        }
        wheelLF->GetCollisionModel()->BuildModel();
        wheelLF->SetCollide(true);
        my_system->AddBody(wheelLF);

        // .. create the revolute joint between the wheel and the spindle
        link_revoluteLF = chrono_types::make_shared<ChLinkLockRevolute>();  // left, front, upper, 1
        link_revoluteLF->Initialize(wheelLF, spindleLF, ChCoordsys<>(ChVector<>(-0.91, -0.026, 1.648965),
                                                                     chrono::Q_from_AngAxis(CH_C_PI / 2, VECT_Y)));
        my_system->AddLink(link_revoluteLF);

        // .. impose distance between two parts (as a massless rod with two spherical joints at the end)
        link_distLFU1 = chrono_types::make_shared<ChLinkDistance>();  // left, front, upper, 1
        link_distLFU1->Initialize(truss, spindleLF, false, ChVector<>(-0.446, 0.245, 1.640965),
                                  ChVector<>(-0.716, 0.215, 1.635965));
        my_system->AddLink(link_distLFU1);

        link_distLFU2 = chrono_types::make_shared<ChLinkDistance>();  // left, front, upper, 2
        link_distLFU2->Initialize(truss, spindleLF, false, ChVector<>(-0.478, 0.196, 1.420965),
                                  ChVector<>(-0.716, 0.215, 1.635965));
        my_system->AddLink(link_distLFU2);

        link_distLFL1 = chrono_types::make_shared<ChLinkDistance>();  // left, front, lower, 1
        link_distLFL1->Initialize(truss, spindleLF, false, ChVector<>(-0.307, 0, 1.911965),
                                  ChVector<>(-0.787, -0.118, 1.652965));
        my_system->AddLink(link_distLFL1);

        link_distLFL2 = chrono_types::make_shared<ChLinkDistance>();  // left, front, lower, 2
        link_distLFL2->Initialize(truss, spindleLF, false, ChVector<>(-0.307, 0, 1.465965),
                                  ChVector<>(-0.787, -0.118, 1.652965));
        my_system->AddLink(link_distLFL2);

        // .. create the spring between the truss and the spindle
        link_springLF = chrono_types::make_shared<ChLinkTSDA>();
        link_springLF->Initialize(truss, spindleLF, false, ChVector<>(-0.498, 0.323, 1.792965),
                                  ChVector<>(-0.543, -0.047, 1.785965));
        link_springLF->SetSpringCoefficient(167062.000);
        link_springLF->SetDampingCoefficient(frontDamping);
        link_springLF->SetRestLength(0.339);
        my_system->AddLink(link_springLF);

        // .. create the rod for steering the wheel
        link_distLSTEER = chrono_types::make_shared<ChLinkDistance>();  // right steer
        link_distLSTEER->Initialize(truss, spindleLF, false, ChVector<>(-0.448, 0.054, 1.438965),
                                    ChVector<>(-0.821, -0.016, 1.512965));
        my_system->AddLink(link_distLSTEER);

        // --- Right Back suspension ---

        // ..the car right-back spindle
        spindleRB = chrono_types::make_shared<ChBody>(chrono_types::make_shared<ChCollisionModelParallel>());
        spindleRB->SetIdentifier(-6);
        spindleRB->SetMass(15.91);
        spindleRB->SetPos(ChVector<>(0.751, -0.026, -1.652965));
        spindleRB->SetRot(ChQuaternion<>(1, 0, 0, 0));
        spindleRB->SetInertiaXX(ChVector<>(4, 2, 2));
        spindleRB->GetCollisionModel()->ClearModel();
        utils::AddBoxGeometry(spindleRB.get(), mat, ChVector<>(0.1, 0.4, 0.4) * 0.5, ChVector<>(0, 0, 0));
        spindleRB->GetCollisionModel()->BuildModel();
        spindleRB->SetCollide(false);
        my_system->AddBody(spindleRB);

        // ..the car right-back wheel
        wheelRB = chrono_types::make_shared<ChBody>(chrono_types::make_shared<ChCollisionModelParallel>());
        wheelRB->SetIdentifier(-7);
        wheelRB->SetMass(3.0);
        wheelRB->SetPos(ChVector<>(0.91, -0.026, -1.652965));
        wheelRB->SetRot(ChQuaternion<>(1, 0, 0, 0));
        wheelRB->SetInertiaXX(ChVector<>(0.2, 0.2, 0.2));
        wheelRB->GetCollisionModel()->ClearModel();
        if (useSpheres) {
            utils::AddSphereGeometry(wheelRB.get(), mat, 0.45, ChVector<>(0, 0, 0));
        } else {
            utils::AddCylinderGeometry(wheelRB.get(), mat, 0.45, 0.2, ChVector<>(0, 0, 0), Q_from_AngZ(CH_C_PI_2));
        }
        wheelRB->GetCollisionModel()->BuildModel();
        wheelRB->SetCollide(true);
        my_system->AddBody(wheelRB);

        // .. create the revolute joint between the wheel and the spindle
        link_revoluteRB = chrono_types::make_shared<ChLinkLockRevolute>();  // right, back, upper, 1
        link_revoluteRB->Initialize(wheelRB, spindleRB, ChCoordsys<>(ChVector<>(0.91, -0.026, -1.652965),
                                                                     chrono::Q_from_AngAxis(CH_C_PI / 2, VECT_Y)));
        my_system->AddLink(link_revoluteRB);

        // .. create the motor transmission joint between the wheel and the truss (assuming small changes of alignment)
        link_engineR = chrono_types::make_shared<ChLinkMotorRotationAngle>();
        link_engineR->Initialize(wheelRB, truss, ChFrame<>(ChVector<>(0.91, -0.026, -1.652965),
                                                              chrono::Q_from_AngAxis(CH_C_PI / 2, VECT_Y)));
        link_engineR->SetSpindleConstraint(ChLinkMotorRotation::SpindleConstraint::OLDHAM);  // approx as a double Rzeppa joint
        link_engineR->SetAngleFunction(chrono_types::make_shared<ChFunction_Ramp>(0, angularSpeed));
        my_system->AddLink(link_engineR);

        // .. impose distance between two parts (as a massless rod with two spherical joints at the end)
        link_distRBU1 = chrono_types::make_shared<ChLinkDistance>();  // right, back, upper, 1
        link_distRBU1->Initialize(truss, spindleRB, false, ChVector<>(0.462, 0.228, -1.339965),
                                  ChVector<>(0.716, 0.216, -1.652965));
        my_system->AddLink(link_distRBU1);

        link_distRBU2 = chrono_types::make_shared<ChLinkDistance>();  // right, back, upper, 2
        link_distRBU2->Initialize(truss, spindleRB, false, ChVector<>(0.462, 0.224, -1.611965),
                                  ChVector<>(0.716, 0.216, -1.652965));
        my_system->AddLink(link_distRBU2);

        link_distRBL1 = chrono_types::make_shared<ChLinkDistance>();  // right, back, lower, 1
        link_distRBL1->Initialize(truss, spindleRB, false, ChVector<>(0.307, 0, -1.465965),
                                  ChVector<>(0.787, -0.118, -1.652965));
        my_system->AddLink(link_distRBL1);

        link_distRBL2 = chrono_types::make_shared<ChLinkDistance>();  // right, back, lower, 2
        link_distRBL2->Initialize(truss, spindleRB, false, ChVector<>(0.307, 0, -1.911965),
                                  ChVector<>(0.787, -0.118, -1.652965));
        my_system->AddLink(link_distRBL2);

        // .. create the spring between the truss and the spindle
        link_springRB = chrono_types::make_shared<ChLinkTSDA>();
        link_springRB->Initialize(truss, spindleRB, false, ChVector<>(0.498, 0.323, -1.79297),
                                  ChVector<>(0.544, -0.038, -1.78597));
        link_springRB->SetSpringCoefficient(369149.000);
        link_springRB->SetDampingCoefficient(rearDamping);
        link_springRB->SetRestLength(0.382);
        my_system->AddLink(link_springRB);

        // .. create the rod for avoid the steering of the wheel
        link_distRBlat = chrono_types::make_shared<ChLinkDistance>();  // right rod
        link_distRBlat->Initialize(truss, spindleRB, false, ChVector<>(0.416, 0.059, -1.465965),
                                   ChVector<>(0.821, -0.009, -1.518965));
        my_system->AddLink(link_distRBlat);

        // --- Left Back suspension ---

        // ..the car right-back spindle
        spindleLB = chrono_types::make_shared<ChBody>(chrono_types::make_shared<ChCollisionModelParallel>());
        spindleLB->SetIdentifier(-8);
        spindleLB->SetMass(15.91);
        spindleLB->SetPos(ChVector<>(-0.751, -0.026, -1.652965));
        spindleLB->SetRot(ChQuaternion<>(1, 0, 0, 0));
        spindleLB->SetInertiaXX(ChVector<>(4, 2, 2));
        spindleLB->GetCollisionModel()->ClearModel();
        utils::AddBoxGeometry(spindleLB.get(), mat, ChVector<>(0.1, 0.4, 0.4) * 0.5, ChVector<>(0, 0, 0));
        spindleLB->GetCollisionModel()->BuildModel();
        spindleLB->SetCollide(false);
        my_system->AddBody(spindleLB);

        // ..the car left-back wheel
        wheelLB = chrono_types::make_shared<ChBody>(chrono_types::make_shared<ChCollisionModelParallel>());
        wheelLB->SetIdentifier(-9);
        wheelLB->SetMass(3.0);
        wheelLB->SetPos(ChVector<>(-0.91, -0.026, -1.652965));
        wheelLB->SetRot(ChQuaternion<>(1, 0, 0, 0));
        wheelLB->SetInertiaXX(ChVector<>(0.2, 0.2, 0.2));
        wheelLB->GetCollisionModel()->ClearModel();
        if (useSpheres) {
            utils::AddSphereGeometry(wheelLB.get(), mat, 0.45, ChVector<>(0, 0, 0));
        } else {
            utils::AddCylinderGeometry(wheelLB.get(), mat, 0.45, 0.2, ChVector<>(0, 0, 0), Q_from_AngZ(CH_C_PI_2));
        }
        wheelLB->GetCollisionModel()->BuildModel();
        wheelLB->SetCollide(true);
        my_system->AddBody(wheelLB);

        // .. create the revolute joint between the wheel and the spindle
        link_revoluteLB = chrono_types::make_shared<ChLinkLockRevolute>();  // left, back, upper, 1
        link_revoluteLB->Initialize(wheelLB, spindleLB, ChCoordsys<>(ChVector<>(-0.91, -0.026, -1.652965),
                                                                     chrono::Q_from_AngAxis(CH_C_PI / 2, VECT_Y)));
        my_system->AddLink(link_revoluteLB);

        // .. create the motor transmission joint between the wheel and the truss (assuming small changes of alignment)
        link_engineL = chrono_types::make_shared<ChLinkMotorRotationAngle>();
        link_engineL->Initialize(wheelLB, truss, ChFrame<>(ChVector<>(-0.91, -0.026, -1.652965),
                                                              chrono::Q_from_AngAxis(CH_C_PI / 2, VECT_Y)));
        link_engineL->SetSpindleConstraint(ChLinkMotorRotation::SpindleConstraint::OLDHAM);  // approx as a double Rzeppa joint
        link_engineL->SetAngleFunction(chrono_types::make_shared<ChFunction_Ramp>(0, angularSpeed));
        my_system->AddLink(link_engineL);

        // .. impose distance between two parts (as a massless rod with two spherical joints at the end)
        link_distLBU1 = chrono_types::make_shared<ChLinkDistance>();  // left, front, upper, 1
        link_distLBU1->Initialize(truss, spindleLB, false, ChVector<>(-0.462, 0.228, -1.339965),
                                  ChVector<>(-0.716, 0.216, -1.652965));
        my_system->AddLink(link_distLBU1);

        link_distLBU2 = chrono_types::make_shared<ChLinkDistance>();  // left, back, upper, 2
        link_distLBU2->Initialize(truss, spindleLB, false, ChVector<>(-0.462, 0.224, -1.611965),
                                  ChVector<>(-0.716, 0.216, -1.652965));
        my_system->AddLink(link_distLBU2);

        link_distLBL1 = chrono_types::make_shared<ChLinkDistance>();  // left, back, lower, 1
        link_distLBL1->Initialize(truss, spindleLB, false, ChVector<>(-0.307, 0, -1.465965),
                                  ChVector<>(-0.787, -0.118, -1.652965));
        my_system->AddLink(link_distLBL1);

        link_distLBL2 = chrono_types::make_shared<ChLinkDistance>();  // left, back, lower, 2
        link_distLBL2->Initialize(truss, spindleLB, false, ChVector<>(-0.307, 0, -1.911965),
                                  ChVector<>(-0.787, -0.118, -1.652965));
        my_system->AddLink(link_distLBL2);

        // .. create the spring between the truss and the spindle
        link_springLB = chrono_types::make_shared<ChLinkTSDA>();
        link_springLB->Initialize(truss, spindleLB, false, ChVector<>(-0.498, 0.323, -1.79297),
                                  ChVector<>(-0.544, -0.038, -1.78597));
        link_springLB->SetSpringCoefficient(369149.000);
        link_springLB->SetDampingCoefficient(rearDamping);
        link_springLB->SetRestLength(0.382);
        my_system->AddLink(link_springLB);

        // .. create the rod for avoid the steering of the wheel
        link_distLBlat = chrono_types::make_shared<ChLinkDistance>();  // right
        link_distLBlat->Initialize(truss, spindleLB, false, ChVector<>(-0.416, 0.059, -1.465965),
                                   ChVector<>(-0.821, -0.009, -1.518965));
        my_system->AddLink(link_distLBlat);
    }

    // Delete the car object, deleting also all bodies corresponding to
    // the various parts and removing them from the physical system.  Also
    // removes constraints from the system.
    ~MySimpleCar() {
        ChSystem* mysystem = spindleRF->GetSystem();  // trick to get the system here
        // When a ChBodySceneNode is removed via ->remove() from Irrlicht 3D scene manager,
        // it is also automatically removed from the ChSystem (the ChSystem::RemoveBody() is
        // automatically called at Irrlicht node deletion - see ChBodySceneNode.h ).

        // For links, just remove them from the ChSystem using ChSystem::RemoveLink()
        mysystem->RemoveLink(link_revoluteRF);
        mysystem->RemoveLink(link_distRFU1);
        mysystem->RemoveLink(link_distRFU2);
        mysystem->RemoveLink(link_distRFL1);
        mysystem->RemoveLink(link_distRFL2);
        mysystem->RemoveLink(link_springRF);
        mysystem->RemoveLink(link_distRSTEER);

        mysystem->RemoveLink(link_revoluteLF);
        mysystem->RemoveLink(link_distLFU1);
        mysystem->RemoveLink(link_distLFU2);
        mysystem->RemoveLink(link_distLFL1);
        mysystem->RemoveLink(link_distLFL2);
        mysystem->RemoveLink(link_springLF);
        mysystem->RemoveLink(link_distLSTEER);

        mysystem->RemoveLink(link_revoluteRB);
        mysystem->RemoveLink(link_distRBU1);
        mysystem->RemoveLink(link_distRBU2);
        mysystem->RemoveLink(link_distRBL1);
        mysystem->RemoveLink(link_distRBL2);
        mysystem->RemoveLink(link_springRB);
        mysystem->RemoveLink(link_distRBlat);
        mysystem->RemoveLink(link_engineR);

        mysystem->RemoveLink(link_revoluteLB);
        mysystem->RemoveLink(link_distLBU1);
        mysystem->RemoveLink(link_distLBU2);
        mysystem->RemoveLink(link_distLBL1);
        mysystem->RemoveLink(link_distLBL2);
        mysystem->RemoveLink(link_springLB);
        mysystem->RemoveLink(link_distLBlat);
        mysystem->RemoveLink(link_engineL);
    }

    // This can be used, at each time step, to compute the actual value of torque
    // transmitted to the wheels, according to gas throttle / speed / gear value.
    // The following is a very simplified model (the torque curve of the motor is linear
    // and no latency or inertial or clutch effects in gear train are considered.)
    double ComputeWheelTorque() {
        // Assume clutch is never used. Given the kinematics of differential,
        // the speed of the engine transmission shaft is the average of the two wheel speeds,
        // multiplied the conic gear transmission ratio inversed:
        double shaftspeed = (1.0 / this->conic_tau) * 0.5 *
                            (this->link_engineL->GetMotorRot_dt() + this->link_engineR->GetMotorRot_dt());
        // The motorspeed is the shaft speed multiplied by gear ratio inversed:
        double motorspeed = (1.0 / this->gear_tau) * shaftspeed;
        // The torque depends on speed-torque curve of the motor: here we assume a
        // very simplified model a bit like in DC motors:
        double motortorque = max_motor_torque - motorspeed * (max_motor_torque / max_motor_speed);
        // Motor torque is linearly modulated by throttle gas value:
        motortorque = motortorque * this->throttle;
        // The torque at motor shaft:
        double shafttorque = motortorque * (1.0 / this->gear_tau);
        // The torque at wheels - for each wheel, given the differential transmission,
        // it is half of the shaft torque  (multiplied the conic gear transmission ratio)
        double singlewheeltorque = 0.5 * shafttorque * (1.0 / this->conic_tau);
        // debug:print infos on screen:
        // GetLog() << "motor torque="<< motortorque<< "  speed=" <<
        // motorspeed << "  wheel torque = " << singlewheeltorque << "\n";
        // If needed, return also the value of wheel torque:
        return singlewheeltorque;
    }
};

// =============================================================================
// Create the granular material
//
// Granular material consisting of identical spheres with specified radius and
// material properties; the spheres are generated in a number of vertical
// layers with locations within each layer obtained using Poisson Disk sampling,
// thus ensuring that no two spheres are closer than twice the radius.
// =============================================================================

int CreateGranularMaterial(ChSystemParallel* system,
                           double pitLocation_z,
                           double pitDepth,
                           double pitLength,
                           double groundWidth) {
    // Parameters for the granular material
    int Id_g = 1;           // start body ID for particles
    double r_g = .04;       // [cm] radius of granular sphers
    double rho_g = 2500.0;  // [g/cm^3] density of granules

    double desiredBulkDensity = 1.3894;  // [g/cm^3] desired bulk density

    float mu_g = 0.5f;

    // -------------------------------------------
    // Create a material for the granular material
    // -------------------------------------------

    auto mat_g = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    mat_g->SetFriction(mu_g);

    // ---------------------------------------------
    // Create a mixture entirely made out of spheres
    // ---------------------------------------------

    // Create the particle generator with a mixture of 100% spheres
    utils::Generator gen(system);

    std::shared_ptr<utils::MixtureIngredient> m1 = gen.AddMixtureIngredient(utils::MixtureType::SPHERE, 1.0);
    m1->setDefaultMaterial(mat_g);
    m1->setDefaultDensity(rho_g);
    m1->setDefaultSize(r_g);

    // Ensure that all generated particle bodies will have positive IDs.
    gen.setBodyIdentifier(Id_g);

    // ----------------------
    // Generate the particles
    // ----------------------

    double r = 1.01 * r_g;
    ChVector<> hdims(groundWidth - 2 * r, 2 * (pitDepth - 2 * r), pitLength - 2 * r);
    ChVector<> center(0, 0, pitLocation_z + pitLength / 3 - 2 * r);

    // while (center.z() < pitDepth)
    //{
    gen.createObjectsBox(utils::SamplingType::POISSON_DISK, 2 * r, center, hdims / 2);
    //  center.z() += 2 * r;
    //}

    // Return the number of generated particles.
    return gen.getTotalNumBodies();
}

// =============================================================================
// Create a bin consisting of five boxes attached to the ground and a mixer
// blade attached through a revolute joint to ground. The mixer is constrained
// to rotate at constant angular velocity.
// =============================================================================
void AddGround(ChSystemParallelNSC* sys) {
    // IDs for the two bodies
    int groundId = -200;

    // Create a common material
    auto mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    mat->SetFriction(1.0f);

    // Create the containing bin (2 x 2 x 1)
    ChVector<> pos(0, -.6, -2);
    double pitLocation_z = 5;
    double pitDepth = .5;
    double pitLength = 8;
    double groundLength = 40;
    double groundWidth = 5;
    double wallHeight = 6;
    double thickness = 0.1;
    auto ground = chrono_types::make_shared<ChBody>(chrono_types::make_shared<ChCollisionModelParallel>());
    ground->SetIdentifier(groundId);
    ground->SetMass(1);
    ground->SetPos(pos);
    ground->SetRot(ChQuaternion<>(1, 0, 0, 0));
    ground->SetCollide(true);
    ground->SetBodyFixed(true);

    ground->GetCollisionModel()->ClearModel();
    utils::AddBoxGeometry(ground.get(), mat, ChVector<>(groundWidth / 2, thickness, 10 * pitLocation_z / 2),
                          ChVector<>(0, -thickness, 0));
    ground->GetCollisionModel()->BuildModel();

    sys->AddBody(ground);

    // CreateGranularMaterial(sys, pitLocation_z, pitDepth, pitLength, groundWidth);
}

// -----------------------------------------------------------------------------
// Create the system, specify simulation parameters, and run simulation loop.
// -----------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    int threads = 16;

    // Simulation parameters
    // ---------------------

    double gravity = 9.81;
    double time_step = 1e-3;
    double time_end = 10;

    double out_fps = 60;

    // Solver settings
    int max_iteration_normal = 0;
    int max_iteration_sliding = 10000;
    int max_iteration_spinning = 0;
    int max_iteration_bilateral = 0;
    double contact_recovery_speed = 1;
    bool clamp_bilaterals = false;
    double bilateral_clamp_speed = 10e30;
    double tolerance = 1e-2;
    bool thread_tuning = false;

    // Create system
    // -------------

    ChSystemParallelNSC msystem;

    // Set gravitational acceleration
    msystem.Set_G_acc(ChVector<>(0, -gravity, 0));

    // Set number of threads.
    int max_threads = omp_get_num_procs();
    if (threads > max_threads)
        threads = max_threads;
    omp_set_num_threads(threads);
    std::cout << "Using " << threads << " threads" << std::endl;

    msystem.GetSettings()->max_threads = threads;
    msystem.GetSettings()->perform_thread_tuning = thread_tuning;

    // Edit system settings
    msystem.GetSettings()->solver.tolerance = tolerance;
    msystem.GetSettings()->solver.max_iteration_bilateral = max_iteration_bilateral;
    msystem.GetSettings()->solver.clamp_bilaterals = clamp_bilaterals;
    msystem.GetSettings()->solver.bilateral_clamp_speed = bilateral_clamp_speed;
    msystem.GetSettings()->solver.solver_mode = SolverMode::SLIDING;
    msystem.GetSettings()->solver.max_iteration_normal = max_iteration_normal;
    msystem.GetSettings()->solver.max_iteration_sliding = max_iteration_sliding;
    msystem.GetSettings()->solver.max_iteration_spinning = max_iteration_spinning;
    msystem.GetSettings()->solver.alpha = 0;
    msystem.GetSettings()->solver.contact_recovery_speed = contact_recovery_speed;
    msystem.SetMaxPenetrationRecoverySpeed(contact_recovery_speed);
    msystem.ChangeSolverType(SolverType::APGDREF);

    msystem.GetSettings()->collision.collision_envelope = 0.01;
    msystem.GetSettings()->collision.bins_per_axis = vec3(10, 10, 10);

    // Create the fixed and moving bodies
    // ----------------------------------

    AddGround(&msystem);
    MySimpleCar* mycar = new MySimpleCar(&msystem);

// Perform the simulation
// ----------------------

#ifdef CHRONO_OPENGL
    opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
    gl_window.Initialize(1280, 720, "Suspension", &msystem);
    gl_window.SetCamera(ChVector<>(0, 0, -10), ChVector<>(0, 0, 0), ChVector<>(0, 1, 0));
#endif

    // Run simulation for specified time
    int num_steps = (int)std::ceil(time_end / time_step);
    int out_steps = (int)std::ceil((1 / time_step) / out_fps);
    double time = 0;
    int out_frame = 0;

    for (int i = 0; i < num_steps; i++) {
        if (i % out_steps == 0) {
            OutputData(&msystem, out_frame, time);
            out_frame++;
        }

#ifdef CHRONO_OPENGL
        if (gl_window.Active()) {
            gl_window.DoStepDynamics(time_step);
            gl_window.Render();
        } else
            break;
#else
        msystem.DoStepDynamics(time_step);
#endif

        time += time_step;
    }

    return 0;
}

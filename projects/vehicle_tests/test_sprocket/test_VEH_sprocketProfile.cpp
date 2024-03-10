#include <cmath>

#include "chrono/physics/ChBody.h"
#include "chrono/physics/ChSystemNSC.h"

#include "chrono/geometry/ChLineArc.h"
#include "chrono/geometry/ChLinePath.h"
#include "chrono/geometry/ChLineSegment.h"

#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"

using namespace chrono;
using namespace chrono::irrlicht;

// -----------------------------------------------------------------------------
std::shared_ptr<ChLinePath> CreateProfile(int num_teeth, double R_T, double R_C, double R) {
    auto profile = chrono_types::make_shared<ChLinePath>();

    double beta = CH_C_2PI / num_teeth;
    double sbeta = std::sin(beta / 2);
    double cbeta = std::cos(beta / 2);
    double y = (R_T * R_T + R_C * R_C - R * R) / (2 * R_C);
    double x = std::sqrt(R_T * R_T - y * y);
    double gamma = std::asin(x / R);

    for (int i = 0; i < num_teeth; ++i) {
        double alpha = -i * beta;
        ChVector3d p0(0, R_C, 0);
        ChVector3d p1(-R_T * sbeta, R_T * cbeta, 0);
        ChVector3d p2(-x, y, 0);
        ChVector3d p3(x, y, 0);
        ChVector3d p4(R_T * sbeta, R_T * cbeta, 0);
        ChQuaternion<> quat;
        quat.SetFromAngleZ(alpha);
        ChMatrix33<> rot(quat);
        p0 = rot * p0;
        p1 = rot * p1;
        p2 = rot * p2;
        p3 = rot * p3;
        p4 = rot * p4;
        ChLineSegment seg1(p1, p2);
        double angle1 = alpha + 1.5 * CH_C_PI - gamma;
        double angle2 = alpha + 1.5 * CH_C_PI + gamma;
        ChLineArc arc(ChCoordsys<>(p0), R, angle1, angle2, true);
        ChLineSegment seg2(p3, p4);
        profile->AddSubLine(seg1);
        profile->AddSubLine(arc);
        profile->AddSubLine(seg2);
    }

    return profile;
}

// -----------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    // Set path to Chrono data directory
    SetChronoDataPath(CHRONO_DATA_DIR);

    ChSystemNSC system;
    system.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);

    // Create a shared material (default properties)
    auto mat = chrono_types::make_shared<ChContactMaterialNSC>();

    // ----------------------
    // Create the ground body
    // ----------------------

    auto ground = chrono_types::make_shared<ChBody>();
    ground->SetIdentifier(-1);
    ground->SetFixed(true);
    system.AddBody(ground);

    // -----------------------------
    // Create the rotating gear body
    // -----------------------------

    ChVector3d gear_loc(0, 0, 0);

    auto gear = chrono_types::make_shared<ChBody>();
    gear->SetIdentifier(0);
    gear->SetPos(gear_loc);
    // gear->SetWvel_loc(ChVector3d(0, 0, 0.1));
    system.AddBody(gear);

    int n_teeth = 10;
    double R_T = 0.2605;
    double R_C = 0.3;
    double R = 0.089;
    double separation = 0.25;

    std::shared_ptr<ChLinePath> gear_profile = CreateProfile(n_teeth, R_T, R_C, R);

    // Add the collision shape to gear
    gear->EnableCollision(true);
    auto gear_ct_shape = chrono_types::make_shared<ChCollisionShapePath2D>(mat, gear_profile);
    gear->AddCollisionShape(gear_ct_shape, ChFrame<>(ChVector3d(0, 0, separation / 2), QUNIT));
    gear->AddCollisionShape(gear_ct_shape, ChFrame<>(ChVector3d(0, 0, -separation / 2), QUNIT));
    gear->GetCollisionModel()->SetSafeMargin(0.02);

    // Add ChVisualShapeLine visualization asset to gear
    auto gear_profile_plus = chrono_types::make_shared<ChVisualShapeLine>();
    gear_profile_plus->SetLineGeometry(gear_profile);
    gear_profile_plus->SetColor(ChColor(1, 0, 0));
    gear->AddVisualShape(gear_profile_plus, ChFrame<>(ChVector3d(0, 0, separation / 2)));

    auto gear_profile_minus = chrono_types::make_shared<ChVisualShapeLine>();
    gear_profile_minus->SetLineGeometry(gear_profile);
    gear_profile_minus->SetColor(ChColor(0, 1, 0));
    gear->AddVisualShape(gear_profile_minus, ChFrame<>(ChVector3d(0, 0, -separation / 2)));

    auto gear_axle = chrono_types::make_shared<ChVisualShapeCylinder>(0.1 * R, separation);
    gear->AddVisualShape(gear_axle);

    // Revolute constraint (gear-ground)
    auto revolute = chrono_types::make_shared<ChLinkLockRevolute>();
    revolute->Initialize(gear, ground, ChFrame<>(gear_loc));
    system.Add(revolute);

    // ---------------------
    // Create a falling body
    // ---------------------

    double pin_hlen = 0.6 * separation;
    double pin_radius = 0.3 * R;
    ChVector3d pin_loc = gear_loc + ChVector3d(0, R_T, 0);

    auto pin = chrono_types::make_shared<ChBody>();
    pin->SetIdentifier(1);
    pin->SetPos(pin_loc);
    // pin->SetFixed(true);
    system.AddBody(pin);

    auto pin_profile = chrono_types::make_shared<ChLinePath>();
    ChLineArc pin_circle(ChCoordsys<>(), pin_radius, CH_C_2PI, 0, false);
    pin_profile->AddSubLine(pin_circle);

    // Add collision shapes to pin
    pin->EnableCollision(true);
    auto pin_ct_shape = chrono_types::make_shared<ChCollisionShapePath2D>(mat, pin_profile);
    pin->AddCollisionShape(pin_ct_shape, ChFrame<>(ChVector3d(0, 0, separation / 2), QUNIT));
    pin->AddCollisionShape(pin_ct_shape, ChFrame<>(ChVector3d(0, 0, -separation / 2), QUNIT));
    pin->GetCollisionModel()->SetSafeMargin(0.02);

    // Add pin visualization
    auto pin_cyl = chrono_types::make_shared<ChVisualShapeCylinder>(pin_radius, 2 * pin_hlen);
    pin_cyl->SetColor(ChColor(0.2f, 0.5f, 0.8f));
    pin->AddVisualShape(pin_cyl);

    // Create the Irrlicht visualization
    auto vis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
    vis->SetWindowSize(800, 600);
    vis->SetWindowTitle("Paths");
    vis->Initialize();
    vis->AddLogo();
    vis->AddSkyBox();
    vis->AddCamera(ChVector3d(0.5, 0.5, -1));
    vis->AddTypicalLights();
    vis->SetSymbolScale(1e-4);
    vis->EnableContactDrawing(ContactsDrawMode::CONTACT_FORCES);
    vis->AttachSystem(&system);

    // ---------------
    // Simulation loop
    // ---------------

    while (vis->Run()) {
        vis->BeginScene();
        vis->Render();
        vis->EndScene();
        system.DoStepDynamics(0.01);
    }

    return 0;
}

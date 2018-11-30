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
// Test for mesh-box collision detection, using a single ANCF toroidal tire.
//
// The coordinate frame respects the ISO standard adopted in Chrono::Vehicle:
// right-handed frame with X pointing towards the front, Y to the left, and Z up
//
// =============================================================================
////#include <float.h>
////unsigned int fp_control_state = _controlfp(_EM_INEXACT, _MCW_EM);

#include <cmath>
#include <cstdio>
#include <vector>

#include "chrono/assets/ChGlyphs.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#include "chrono/fea/ChContactSurfaceNodeCloud.h"

#include "chrono_vehicle/ChConfigVehicle.h"
#include "chrono_vehicle/ChVehicleModelData.h"
#include "chrono_vehicle/terrain/RigidTerrain.h"
#include "chrono_vehicle/wheeled_vehicle/tire/ANCFToroidalTire.h"

#include "chrono_thirdparty/filesystem/path.h"

#include "../BaseTest.h"

using namespace chrono;
using namespace chrono::vehicle;

// =============================================================================
// Global definitions
// =============================================================================

// Wheel location
ChVector<> wheel_pos(0, 0, 0);

// Wheel orientation
ChQuaternion<> wheeel_rot = Q_from_AngX(CH_C_DEG_TO_RAD * 30);

// Offset of the lowest tire vertex above terrain
// (a negative value results in contact)
double tire_offset = 0;

// Radius of contact node (associated with mesh vertices)
double node_radius = 0.01;

// Rigid terrain dimensions
double terrain_length = 100.0;  // size in X direction
double terrain_width = 2.0;     // size in Y direction

// Render contact points (from default collider or from custom collider)
enum ColliderType { DEFAULT_COLLIDER, CUSTOM_COLLIDER };
ColliderType render_which = CUSTOM_COLLIDER;
double point_size = 0.003;

// Tire mesh wireframe only?
bool tire_mesh_wireframe = false;

// =============================================================================

// Test class
class toroidalTireTest : public BaseTest {
  public:
    toroidalTireTest(const std::string& testName, const std::string& testProjectName)
        : BaseTest(testName, testProjectName), m_execTime(0) {}

    ~toroidalTireTest() {}

    // Override corresponding functions in BaseTest
    virtual bool execute() override;
    virtual double getExecutionTime() const override { return m_execTime; }

  private:
    double m_execTime;
};

// =============================================================================
// Contact reporter class
// =============================================================================

class MyContactReporter : public ChReportContactCallback {
  public:
    MyContactReporter(std::shared_ptr<ChBody> ground) : m_ground(ground), m_num_contacts(0) {}

    unsigned int GetNumAddedContacts() const { return m_num_contacts; }
    const std::vector<ChVector<>>& GetTerrainPoints() const { return m_terrain_points; }
    const std::vector<ChVector<>>& GetNodePoints() const { return m_node_points; }

  private:
    virtual bool ReportContactCallback(const ChVector<>& pA,
                                       const ChVector<>& pB,
                                       const ChMatrix33<>& plane_coord,
                                       const double& distance,
                                       const ChVector<>& react_forces,
                                       const ChVector<>& react_torques,
                                       ChContactable* objA,
                                       ChContactable* objB) override {
        m_num_contacts++;
        ChVector<> normal = plane_coord.Get_A_Xaxis();
        ChVector<> pointT = (objA == m_ground.get()) ? pA : pB;
        ChVector<> pointN = (objA == m_ground.get()) ? pB : pA;
        printf("%3d | ", m_num_contacts);
        printf("%+10.5e | ", distance);
        printf("%+10.5e  %+10.5e  %+10.5e | ", pointT.x, pointT.y, pointT.z);
        printf("%+10.5e  %+10.5e  %+10.5e | ", pointN.x, pointN.y, pointN.z);
        printf("%+10.5e  %+10.5e  %+10.5e \n", normal.x, normal.y, normal.z);

        // Cache the contact points
        m_terrain_points.push_back(pointT);
        m_node_points.push_back(pointN);

        return true;
    }

    std::shared_ptr<ChBody> m_ground;
    int m_num_contacts;
    std::vector<ChVector<>> m_node_points;
    std::vector<ChVector<>> m_terrain_points;
};

// =============================================================================
// Custom collision detection class
// =============================================================================

class TireTestCollisionManager : public ChSystem::ChCustomComputeCollisionCallback {
  public:
    TireTestCollisionManager(std::shared_ptr<fea::ChContactSurfaceNodeCloud> surface,
                             std::shared_ptr<RigidTerrain> terrain,
                             double radius)
        : m_surface(surface), m_terrain(terrain), m_radius(radius), m_num_contacts(0) {}

    unsigned int GetNumAddedContacts() const { return m_num_contacts; }
    const std::vector<ChVector<>>& GetTerrainPoints() const { return m_terrain_points; }
    const std::vector<ChVector<>>& GetNodePoints() const { return m_node_points; }

  private:
    virtual void PerformCustomCollision(ChSystem* system) override {
        printf("\n>>>> Custom collision detection <<<\n\n");

        for (unsigned int in = 0; in < m_surface->GetNnodes(); in++) {
            // Represent the contact node as a sphere (P, m_radius)
            auto contact_node = std::static_pointer_cast<fea::ChContactNodeXYZsphere>(m_surface->GetNode(in));
            const ChVector<>& P = contact_node->GetNode()->GetPos();

            // Represent the terrain as a plane (Q, normal)
            ChVector<> normal = m_terrain->GetNormal(P.x, P.y);
            ChVector<> Q(P.x, P.y, m_terrain->GetHeight(P.x, P.y));

            // Calculate signed height of sphere center above plane
            double height = Vdot(normal, P - Q);

            // No collision if the sphere center is above plane by more than radius
            if (height >= m_radius)
                continue;

            m_num_contacts++;

            // Create a collision info structure:
            //    modelA: terrain collision model
            //    modelB: node collision model
            //    vN: normal (from A to B)
            //    vpA: contact point on terrain
            //    vpB: contact point on node
            //    distance: penetration (negative)
            collision::ChCollisionInfo contact;
            contact.modelA = m_terrain->GetGroundBody()->GetCollisionModel();
            contact.modelB = contact_node->GetCollisionModel();
            contact.vN = normal;
            contact.vpA = P - height * normal;
            contact.vpB = P - m_radius * normal;
            contact.distance = height - m_radius;

            printf("%3d | ", m_num_contacts);
            printf("%+10.5e | ", contact.distance);
            printf("%+10.5e  %+10.5e  %+10.5e | ", contact.vpA.x, contact.vpA.y, contact.vpA.z);
            printf("%+10.5e  %+10.5e  %+10.5e | ", contact.vpB.x, contact.vpB.y, contact.vpB.z);
            printf("%+10.5e  %+10.5e  %+10.5e \n", contact.vN.x, contact.vN.y, contact.vN.z);

            // Note: do not register the new contact
            ////system->GetContactContainer()->AddContact(contact);

            // Cache the contact points (on terrain and nodes)
            m_terrain_points.push_back(contact.vpA);
            m_node_points.push_back(contact.vpB);
        }

        printf("\nFound: %d\n\n", m_num_contacts);
    }

    std::shared_ptr<fea::ChContactSurfaceNodeCloud> m_surface;
    std::shared_ptr<RigidTerrain> m_terrain;
    double m_radius;
    unsigned int m_num_contacts;
    std::vector<ChVector<>> m_node_points;
    std::vector<ChVector<>> m_terrain_points;
};

bool toroidalTireTest::execute() {
    // Create the mechanical system
    // ----------------------------

    ChSystemDEM system;
    system.Set_G_acc(ChVector<>(0.0, 0.0, 0.0));

    // Create the wheel (rim)
    // ----------------------

    auto wheel = std::make_shared<ChBody>(ChMaterialSurfaceBase::DEM);
    system.AddBody(wheel);
    wheel->SetIdentifier(2);
    wheel->SetName("wheel");
    wheel->SetBodyFixed(false);
    wheel->SetCollide(false);
    wheel->SetMass(10);
    wheel->SetInertiaXX(ChVector<>(1, 1, 1));
    wheel->SetPos(wheel_pos);
    wheel->SetRot(wheeel_rot);

    // Create the tire
    // ---------------

    auto tire = std::make_shared<ANCFToroidalTire>("ANCF_Tire");

    tire->EnablePressure(true);
    tire->EnableRimConnection(true);
    tire->EnableContact(true);

    tire->SetContactSurfaceType(ChANCFTire::NODE_CLOUD);
    tire->SetContactNodeRadius(node_radius);

    tire->Initialize(wheel, LEFT);

    tire->SetVisualizationType(VisualizationType::NONE);

    double tire_radius = tire->GetRadius();
    double rim_radius = tire->GetRimRadius();
    double tire_width = tire->GetWidth();

    // Find lowest mesh node
    auto tire_mesh = tire->GetMesh();
    double z_min = 0;
    for (unsigned int in = 0; in < tire_mesh->GetNnodes(); in++) {
        auto node = std::dynamic_pointer_cast<fea::ChNodeFEAxyz>(tire_mesh->GetNode(in));
        if (node->GetPos().z < z_min)
            z_min = node->GetPos().z;
    }

    // Create the terrain
    // ------------------

    auto terrain = std::make_shared<RigidTerrain>(&system);
    terrain->SetContactFrictionCoefficient(0.9f);
    terrain->SetContactRestitutionCoefficient(0.01f);
    terrain->SetContactMaterialProperties(2e7f, 0.3f);
    ////terrain->SetTexture(vehicle::GetDataFile("terrain/textures/tile4.jpg"), 200, 4);
    terrain->Initialize(z_min - tire_offset, terrain_length, terrain_width);

    // Create the custom collision detector
    // ------------------------------------

    // Extract the contact surface from the tire mesh
    auto surface = std::dynamic_pointer_cast<fea::ChContactSurfaceNodeCloud>(tire_mesh->GetContactSurface(0));

    // Add custom collision callback
    TireTestCollisionManager collider(surface, terrain, tire->GetContactNodeRadius());
    system.SetCustomComputeCollisionCallback(&collider);

    // Complete system setup
    // ---------------------

    system.SetupInitial();

    // Perform collision detection
    // ---------------------------
    system.Update(true);
    system.ComputeCollisions();

    // Report tire-terrain contacts
    MyContactReporter reporter(terrain->GetGroundBody());
    printf("\n>>>> Default collision detection <<<\n\n");
    system.GetContactContainer()->ReportAllContacts(&reporter);
    printf("\nFound: %d\n\n", system.GetContactContainer()->GetNcontacts());

    // Asset for rendering contact points
    unsigned int npoints = collider.GetNumAddedContacts();
    std::vector<ChVector<>> terrain_points;
    std::vector<ChVector<>> node_points;

    switch (render_which) {
        case DEFAULT_COLLIDER:
            printf("\nRender contact points from DEFAULT COLLIDER (Bullet)\n");
            terrain_points = reporter.GetTerrainPoints();
            node_points = reporter.GetNodePoints();
            break;
        case CUSTOM_COLLIDER:
            printf("\nRender contact points from CUSTOM COLLIDER (exact)\n");
            terrain_points = collider.GetTerrainPoints();
            node_points = collider.GetNodePoints();
            break;
    }

    auto glyph_points = std::make_shared<ChGlyphs>();
    glyph_points->Reserve(2 * npoints);
    glyph_points->SetGlyphsSize(point_size);
    for (unsigned int ip = 0; ip < npoints; ip++) {
        glyph_points->SetGlyphPoint(ip, wheel->TransformPointParentToLocal(terrain_points[ip]), ChColor(1, 0, 0));
        glyph_points->SetGlyphPoint(npoints + ip, wheel->TransformPointParentToLocal(node_points[ip]),
                                    ChColor(0, 1, 0));
    }

    wheel->AddAsset(glyph_points);

    // Rendering loop
    // --------------
    int num_steps = 100;  // Is this enough?
    double step_size = 1e-3;
    double time_total = 0;
    for (int istep = 0; istep < num_steps; istep++) {
        system.DoStepDynamics(step_size);
        time_total += system.GetTimerStep();
    }
    m_execTime = time_total;
    return true;
}
// =============================================================================
// Main driver program
// =============================================================================

int main(int argc, char* argv[]) {
    std::string out_dir = "../METRICS";
    if (!filesystem::create_directory(filesystem::path(out_dir))) {
        std::cout << "Error creating directory " << out_dir << std::endl;
        return 1;
    }

    toroidalTireTest test("metrics_VEH_collisionToroidalTire", "Chrono::VEH");
    test.setOutDir(out_dir);
    test.setVerbose(true);
    bool passed = test.run();
    test.print();
    return 0;
}

//
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
// Test for co-simulation between Chrono::FEA and Chrono::Multicore
//
// =============================================================================

#include "chrono/geometry/ChTriangleMeshConnected.h"
#include "chrono/physics/ChLoadBodyMesh.h"
#include "chrono/physics/ChLoadContainer.h"
#include "chrono/physics/ChLoaderUV.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/solver/ChIterativeSolverLS.h"

#include "chrono/fea/ChLoadContactSurfaceMesh.h"
#include "chrono/fea/ChMesh.h"
#include "chrono/fea/ChMeshFileLoader.h"
#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"

#include "chrono_multicore/physics/ChSystemMulticore.h"

#include "chrono/ChConfig.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsGenerators.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#ifdef CHRONO_OPENGL
#include "chrono_opengl/ChVisualSystemOpenGL.h"
#endif

using namespace chrono;
using namespace chrono::fea;
using namespace chrono::irrlicht;

// Utility to draw some triangles that are affected by cosimulation.
// Also plot forces as vectors.
// Mostly for debugging.
void draw_affected_triangles(ChVisualSystemIrrlicht& vis,
                             std::vector<ChVector3d>& vert_pos,
                             std::vector<ChVector3i>& triangles,
                             std::vector<int>& vert_indexes,
                             std::vector<ChVector3d>& vert_forces,
                             double forcescale = 0.01) {
    for (int it = 0; it < triangles.size(); ++it) {
        bool vert_hit = false;
        for (int io = 0; io < vert_indexes.size(); ++io) {
            if (triangles[it].x() == vert_indexes[io] || triangles[it].y() == vert_indexes[io] ||
                triangles[it].z() == vert_indexes[io])
                vert_hit = true;
        }
        if (vert_hit == true) {
            std::vector<chrono::ChVector3d> fourpoints = {vert_pos[triangles[it].x()], vert_pos[triangles[it].y()],
                                                          vert_pos[triangles[it].z()], vert_pos[triangles[it].x()]};
            tools::drawPolyline(&vis, fourpoints, ChColor(0.9f, 0.7f, 0), true);
        }
    }
    if (forcescale > 0)
        for (int io = 0; io < vert_indexes.size(); ++io) {
            std::vector<chrono::ChVector3d> forceline = {vert_pos[vert_indexes[io]],
                                                         vert_pos[vert_indexes[io]] + vert_forces[io] * forcescale};
            tools::drawPolyline(&vis, forceline, ChColor(0.9f, 0, 0), true);
        }
}

int main(int argc, char* argv[]) {
    // Set path to Chrono data directory
    SetChronoDataPath(CHRONO_DATA_DIR);

    double time_step = 0.005;
    int max_iteration = 30;
    double tolerance = 1e-3;
    double time_end = 2.0;

    // Frequency for visualization output
    bool saveData = true;
    int out_fps = 60;

    // Global parameter for tire:
    double tire_rad = 0.8;
    double tire_vel_z0 = -3;
    ChVector3d tire_center(0, 1 + 0.02 + tire_rad, 0);
    ChMatrix33<> tire_alignment(QuatFromAngleAxis(CH_PI, VECT_Y));  // create rotated 180� on y

    double tire_w0 = tire_vel_z0 / tire_rad;

    // Create a Chrono::Engine physical system
    ChSystemSMC sys;
    sys.SetCollisionSystemType(ChCollisionSystem::Type::MULTICORE);

    //
    // CREATE A FINITE ELEMENT MESH
    //

    // Create the surface material, containing information
    // about friction etc.

    auto mysurfmaterial = chrono_types::make_shared<ChContactMaterialSMC>();
    mysurfmaterial->SetYoungModulus(10e4);
    mysurfmaterial->SetFriction(0.3f);
    mysurfmaterial->SetRestitution(0.2f);
    mysurfmaterial->SetAdhesion(0);

    // Create a mesh, that is a container for groups
    // of FEA elements and their referenced nodes.
    auto my_mesh = chrono_types::make_shared<ChMesh>();
    sys.Add(my_mesh);

    // Create a material, that must be assigned to each solid element in the mesh,
    // and set its parameters
    auto mmaterial = chrono_types::make_shared<ChContinuumElastic>();
    mmaterial->SetYoungModulus(0.003e9);  // rubber 0.01e9, steel 200e9
    mmaterial->SetPoissonRatio(0.4);
    mmaterial->SetRayleighDampingBeta(0.004);
    mmaterial->SetDensity(1000);

    // Load an ABAQUS .INP tetahedron mesh file from disk, defining a tetahedron mesh.
    // Note that not all features of INP files are supported. Also, quadratic tetahedrons are promoted to linear.
    // This is much easier than creating all nodes and elements via C++ programming.
    // Ex. you can generate these .INP files using Abaqus or exporting from the SolidWorks simulation tool.
    std::map<std::string, std::vector<std::shared_ptr<ChNodeFEAbase>>> node_sets;
    try {
        ChMeshFileLoader::FromAbaqusFile(my_mesh, GetChronoDataFile("fea/tractor_wheel_coarse.INP").c_str(), mmaterial,
                                         node_sets, tire_center, tire_alignment);
    } catch (std::exception& myerr) {
        std::cout << myerr.what();
        return 0;
    }

    // Create the contact surface(s).
    // Use the AddFacesFromBoundary() to select automatically the outer skin of the tetrahedron mesh:
    auto mcontactsurf = chrono_types::make_shared<ChContactSurfaceMesh>(mysurfmaterial);
    my_mesh->AddContactSurface(mcontactsurf);
    mcontactsurf->AddFacesFromBoundary(*my_mesh);

    /// Create a mesh load for cosimulation, acting on the contact surface above
    /// (forces on nodes will be computed by an external procedure)
    auto mloadcontainer = chrono_types::make_shared<ChLoadContainer>();
    sys.Add(mloadcontainer);
    auto mrigidmeshload = chrono_types::make_shared<ChLoadContactSurfaceMesh>(mcontactsurf);
    mloadcontainer->Add(mrigidmeshload);

    // ==Asset== attach a visualization of the FEM mesh.
    // This will automatically update a triangle mesh (a ChVisualShapeTriangleMesh
    // asset that is internally managed) by setting  proper
    // coordinates and vertex colours as in the FEM elements.
    auto mvisualizemesh = chrono_types::make_shared<ChVisualShapeFEA>();
    mvisualizemesh->SetFEMdataType(ChVisualShapeFEA::DataType::NODE_SPEED_NORM);
    mvisualizemesh->SetColormapRange(0.0, 10);
    mvisualizemesh->SetSmoothFaces(true);
    my_mesh->AddVisualShapeFEA(mvisualizemesh);

    //
    // END CREATE A FINITE ELEMENT MESH
    //

    /*
    //
    // CREATE A RIGID BODY WITH A MESH
    //

    // Create also a rigid body with a rigid mesh that will be used for the cosimulation,
    // this time the ChLoadContactSurfaceMesh cannot be used as in the FEA case, so we
    // will use the ChLoadBodyMesh class:

    auto mrigidbody = chrono_types::make_shared<ChBody>();
    sys.Add(mrigidbody);
    mrigidbody->SetMass(200);
    mrigidbody->SetInertiaXX(ChVector3d(20,20,20));
    mrigidbody->SetPos(tire_center);

    auto mrigidmesh = chrono_types::make_shared<ChVisualShapeTriangleMesh>();
    mrigidmesh->GetMesh().LoadWavefrontMesh(GetChronoDataFile("models/tractor_wheel_fine.obj"));
    mrigidmesh->GetMesh().Transform(VNULL, QuatFromAngleAxis(CH_PI, VECT_Y) );
    mrigidbody->AddAsset(mrigidmesh);

    auto mcol = chrono_types::make_shared<ChColorAsset>();
    mcol->SetColor(ChColor(0.3f, 0.3f, 0.3f));
    mrigidbody->AddAsset(mcol);

    /// Create a mesh load for cosimulation, acting on the contact surface above
    /// (forces on nodes will be computed by an external procedure)

    auto mloadcontainer = chrono_types::make_shared<ChLoadContainer>();
    sys.Add(mloadcontainer);

    // this is used to use the mesh in cosimulation!
    auto mrigidmeshload = chrono_types::make_shared<ChLoadBodyMesh>(mrigidbody, mrigidmesh->GetMesh());
    mloadcontainer->Add(mrigidmeshload);

    //
    // END CREATE A RIGID BODY WITH A MESH
    //
     *
    */

#ifndef CHRONO_OPENGL
    auto vis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
    vis->SetWindowSize(800, 600);
    vis->SetWindowTitle("FEA contacts");
    vis->Initialize();
    vis->AddLogo();
    vis->AddSkyBox();
    vis->AddCamera(ChVector3d(3, 1.4, -3.2));
    vis->AddTypicalLights();
    vis->AttachSystem(&sys);
#endif

    //
    // THE SOFT-REAL-TIME CYCLE
    //

    // Change solver to embedded MINRES
    // NOTE! it is strongly advised that you compile the optional PardisoMKL module
    // if you need higher precision, and switch to its Pardiso solver - see demos for FEA & PardisoMKL.
    auto solver = chrono_types::make_shared<ChSolverMINRES>();
    solver->EnableWarmStart(true);
    solver->SetMaxIterations(40);
    solver->SetTolerance(1e-12);
    sys.SetSolver(solver);

    // Change type of integrator:
    sys.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT_LINEARIZED);  // fast, less precise

    // BEGIN MULTICORE SYSTEM INITIALIZATION
    ChSystemMulticoreNSC* systemG = new ChSystemMulticoreNSC();

    // Set gravitational acceleration
    systemG->SetGravitationalAcceleration(sys.GetGravitationalAcceleration());

    // Set solver parameters
    systemG->GetSettings()->solver.solver_mode = SolverMode::SLIDING;
    systemG->GetSettings()->solver.max_iteration_normal = max_iteration / 3;
    systemG->GetSettings()->solver.max_iteration_sliding = max_iteration / 3;
    systemG->GetSettings()->solver.max_iteration_spinning = 0;
    systemG->GetSettings()->solver.max_iteration_bilateral = max_iteration / 3;
    systemG->GetSettings()->solver.tolerance = tolerance;
    systemG->GetSettings()->solver.alpha = 0;
    systemG->GetSettings()->solver.contact_recovery_speed = 10000;
    systemG->ChangeSolverType(SolverType::APGD);
    systemG->GetSettings()->collision.narrowphase_algorithm = ChNarrowphase::Algorithm::HYBRID;

    systemG->GetSettings()->collision.collision_envelope = 0.01;
    systemG->GetSettings()->collision.bins_per_axis = vec3(10, 10, 10);

    auto triMat = chrono_types::make_shared<ChContactMaterialNSC>();
    triMat->SetFriction(0.4f);

    // Create the triangles for the tire geometry
    ChVector3d pos(0, 0, 0);
    ChVector3d vel(0, 0, 0);

    std::vector<ChVector3d> vert_pos;
    std::vector<ChVector3d> vert_vel;
    std::vector<ChVector3i> triangles;
    std::vector<ChVector3d> vert_forces;
    std::vector<int> vert_indexes;
    vert_forces.clear();
    vert_indexes.clear();
    std::vector<ChVector3d> vert_forcesVisualization;
    std::vector<int> vert_indexesVisualization;
    vert_forcesVisualization.clear();
    vert_indexesVisualization.clear();
    mrigidmeshload->OutputSimpleMesh(vert_pos, vert_vel, triangles);
    for (int i = 0; i < vert_pos.size(); i++) {
        vert_forces.push_back(ChVector3d(0, 0, 0));
        vert_indexes.push_back(i);
    }

    double mass = 2;  // mrigidbody->GetMass()/((double) triangles.size());
    double radius = 0.005;
    ChVector3d inertia = (2.0 / 5.0) * mass * radius * radius * ChVector3d(1, 1, 1);

    int tri_tag = 0;
    std::vector<std::shared_ptr<ChBody>> tri_bodies;
    for (int i = 0; i < triangles.size(); i++) {
        auto triangle = chrono_types::make_shared<ChBody>();
        triangle->SetTag(tri_tag++);
        triangle->SetMass(mass);
        triangle->SetInertiaXX(inertia);
        pos = (vert_pos[triangles[i].x()] + vert_pos[triangles[i].y()] + vert_pos[triangles[i].z()]) / 3.0;
        vel = (vert_vel[triangles[i].x()] + vert_vel[triangles[i].y()] + vert_vel[triangles[i].z()]) / 3.0;
        triangle->SetPos(pos);
        triangle->SetLinVel(vel);
        triangle->SetRot(ChQuaternion<>(1, 0, 0, 0));
        triangle->EnableCollision(true);
        triangle->SetFixed(true);

        std::string name = "tri" + std::to_string(tri_tag);
        chrono::utils::AddTriangleGeometry(triangle.get(), triMat, vert_pos[triangles[i].x()] - pos,
                                           vert_pos[triangles[i].y()] - pos, vert_pos[triangles[i].z()] - pos, name);
        triangle->GetCollisionModel()->SetFamily(1);
        triangle->GetCollisionModel()->DisallowCollisionsWith(1);

        tri_bodies.push_back(triangle);
        systemG->AddBody(triangle);
    }

    // Add the terrain, MUST BE ADDED AFTER TIRE GEOMETRY (for index assumptions)
    chrono::utils::CreateBoxContainer(systemG, triMat, ChVector3d(2, 2, 2), 0.2, ChVector3d(0, -1, 0), QUNIT, true,
                                      true, false);

    double r = 0.1;  // 0.02;//
    double shapeRatio = 0.4;
    chrono::utils::ChPDSampler<double> sampler(2 * r);
    chrono::utils::ChGenerator gen(systemG);
    auto m1 = gen.AddMixtureIngredient(chrono::utils::MixtureType::ELLIPSOID, 1.0);
    m1->SetDefaultMaterial(triMat);
    m1->SetDefaultDensity(2500);
    m1->SetDefaultSize(ChVector3d(r, r * shapeRatio, r));

    ChVector3d hdims(1 - r * 1.01, 0.5, 1 - r * 1.01);
    ChVector3d center(0, 0, 0);
    gen.CreateObjectsBox(sampler, center, hdims);

#ifdef CHRONO_OPENGL
    // Initialize OpenGL
    opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
    gl_window.Initialize(1280, 720, "DEMO TTI", systemG);
    gl_window.SetCamera(ChVector3d(1, 1.4, -1.2), ChVector3d(0, tire_rad, 0), ChVector3d(0, 1, 0));
    gl_window.SetRenderMode(opengl::WIREFRAME);
#endif
    // END MULTICORE SYSTEM INITIALIZATION

    // Begin time loop
    int out_steps = (int)std::ceil((1.0 / time_step) / out_fps);
    int timeIndex = 0;
    double time = 0;
    int frameIndex = 0;
#ifdef CHRONO_OPENGL
    while (true) {
#else
    while (vis->Run()) {
#endif
// while (time<time_end) {

// STEP 1: ADVANCE DYNAMICS OF GRANULAR SYSTEM
#ifdef CHRONO_OPENGL
        if (gl_window.Active()) {
            gl_window.DoStepDynamics(time_step);
            gl_window.Render();
        } else
            break;
#else
        systemG->DoStepDynamics(time_step);

        if (timeIndex % out_steps == 0 && saveData) {
            char filename[100];
            sprintf(filename, "../POVRAY/data_%d.dat", frameIndex);

            chrono::utils::WriteVisualizationAssets(systemG, filename, false);
            std::string delim = ",";
            chrono::utils::ChWriterCSV csv(delim);
            csv << triangles.size() << std::endl;
            for (int i = 0; i < triangles.size(); i++) {
                csv << systemG->GetBodies().at(i)->GetPos() << vert_pos[triangles[i].x()]
                    << vert_pos[triangles[i].y()] << vert_pos[triangles[i].z()] << std::endl;
            }
            sprintf(filename, "../POVRAY/triangles_%d.dat", frameIndex);
            csv.WriteToFile(filename);
        }
#endif
        // END STEP 1

        // STEP 2: APPLY CONTACT FORCES FROM GRANULAR TO TIRE SYSTEM
        real3 force(0, 0, 0);
        real3 torque(0, 0, 0);
        systemG->CalculateContactForces();

        vert_forces.clear();
        for (int i = 0; i < vert_pos.size(); i++) {
            vert_forces.push_back(ChVector3d(0, 0, 0));
        }

        for (int i = 0; i < triangles.size(); i++) {
            force = systemG->GetBodyContactForce(tri_bodies[i]);
            torque = systemG->GetBodyContactTorque(tri_bodies[i]);

            // TODO: Calculate force based on the position in the triangle
            vert_forces[triangles[i].x()] += ChVector3d(force.x, force.y, force.z) / 3;
            vert_forces[triangles[i].y()] += ChVector3d(force.x, force.y, force.z) / 3;
            vert_forces[triangles[i].z()] += ChVector3d(force.x, force.y, force.z) / 3;
        }
        mrigidmeshload->InputSimpleForces(vert_forces, vert_indexes);
// END STEP 2

// STEP 3: ADVANCE DYNAMICS OF TIRE SYSTEM
#ifdef CHRONO_OPENGL
        sys.DoStepDynamics(time_step);
#else
        vis->BeginScene();
        vis->Render();

        sys.DoStepDynamics(time_step);

        if (timeIndex % out_steps == 0 && saveData) {
            // takeScreenshot(application.GetDevice(),frameIndex);
            frameIndex++;
        }
#endif
        // END STEP 3

        // STEP 4: UPDATE THE POSITION/VELOCITY OF THE TIRE GEOMETRY IN GRANULAR SYSTEM
        vert_pos.clear();
        vert_vel.clear();
        triangles.clear();

        mrigidmeshload->OutputSimpleMesh(vert_pos, vert_vel, triangles);

        for (int i = 0; i < triangles.size(); i++) {
            std::shared_ptr<ChBody> triBody = systemG->GetBodies().at(i);
            pos = (vert_pos[triangles[i].x()] + vert_pos[triangles[i].y()] + vert_pos[triangles[i].z()]) / 3.0;
            triBody->SetPos(pos);
            vel = (vert_vel[triangles[i].x()] + vert_vel[triangles[i].y()] + vert_vel[triangles[i].z()]) / 3.0;
            triBody->SetLinVel(vel);

            //            // Update visual assets TODO: chrono_opengl cannot handle dynamic meshes yet
            //            for (int j = 0; j < triBody->GetAssets().size(); j++) {
            //              std::shared_ptr<ChAsset> asset = triBody->GetAssets()[j];
            //              if (std::dynamic_pointer_cast<ChVisualShapeTriangleMesh>(asset)) {
            //                //std::cout << j << std::endl;
            //                //std::cout << vert_pos[triangles[i].x()].x() << " " << vert_pos[triangles[i].x()].y() <<
            //                " " << vert_pos[triangles[i].x()].z() << std::endl;
            //                ((ChVisualShapeTriangleMesh*)(asset.get()))->GetMesh().m_vertices[0] =
            //                vert_pos[triangles[i].x()];
            //                ((ChVisualShapeTriangleMesh*)(asset.get()))->GetMesh().m_vertices[1] =
            //                vert_pos[triangles[i].y()];
            //                ((ChVisualShapeTriangleMesh*)(asset.get()))->GetMesh().m_vertices[2] =
            //                ChVector3d(0);//vert_pos[triangles[i].z()];
            //              }
            //            }

            // Update collision information
            systemG->data_manager->cd_data->shape_data.triangle_rigid[3 * i + 0] =
                real3(vert_pos[triangles[i].x()].x() - pos.x(), vert_pos[triangles[i].x()].y() - pos.y(),
                      vert_pos[triangles[i].x()].z() - pos.z());
            systemG->data_manager->cd_data->shape_data.triangle_rigid[3 * i + 1] =
                real3(vert_pos[triangles[i].y()].x() - pos.x(), vert_pos[triangles[i].y()].y() - pos.y(),
                      vert_pos[triangles[i].y()].z() - pos.z());
            systemG->data_manager->cd_data->shape_data.triangle_rigid[3 * i + 2] =
                real3(vert_pos[triangles[i].z()].x() - pos.x(), vert_pos[triangles[i].z()].y() - pos.y(),
                      vert_pos[triangles[i].z()].z() - pos.z());
        }
        // END STEP 4

#ifndef CHRONO_OPENGL
        // now, just for debugging and some fun, draw some triangles
        // (only those that have a vertex that has a force applied):
        vert_forcesVisualization.clear();
        vert_indexesVisualization.clear();
        for (int i = 0; i < vert_forces.size(); i++) {
            if (vert_forces[i].Length() > 1e-5) {
                vert_forcesVisualization.push_back(vert_forces[i]);
                vert_indexesVisualization.push_back(vert_indexes[i]);
            }
        }
        draw_affected_triangles(*vis, vert_pos, triangles, vert_indexesVisualization, vert_forcesVisualization,
                                0.01);

        // End of cosimulation block
        // -------------------------------------------------------------------------

        tools::drawGrid(vis.get(), 0.1, 0.1, 20, 20, ChCoordsys<>(VNULL, CH_PI_2, VECT_X), ChColor(0.4f, 0.4f, 0.4f),
                        true);

        vis->EndScene();
#endif
        timeIndex++;
        time += time_step;
        std::cout << time << std::endl;
    }

    return 0;
}

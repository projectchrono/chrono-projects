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

///////////////////////////////////////////////////
//
//   Demo code about
//
//     - hard-sphere (DVI) direct shear box test
//
//     - written by Jonathan Fleischmann 01/30/2015
//
//   CHRONO-PARALLEL
//   ------
//   Multibody dynamics engine
//
///////////////////////////////////////////////////

#include <iostream>
#include <vector>
#include <valarray>
#include <string>
#include <sstream>
#include <cmath>

#include "core/ChFileutils.h"
#include "core/ChStream.h"

#include "chrono_utils/ChUtilsCreators.h"
#include "chrono_utils/ChUtilsGenerators.h"
#include "chrono_utils/ChUtilsInputOutput.h"

#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_parallel/lcp/ChLcpSystemDescriptorParallel.h"

// Control use of OpenGL run-time rendering
//#undef CHRONO_PARALLEL_HAS_OPENGL

#ifdef CHRONO_PARALLEL_HAS_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif

using namespace chrono;
using namespace chrono::collision;

using std::cout;
using std::flush;
using std::endl;

// Utility for adding (visible or invisible) walls

void AddWall(ChSharedBodyPtr& body, const ChVector<>& dim, const ChVector<>& loc, bool visible) {
  body->GetCollisionModel()->AddBox(dim.x, dim.y, dim.z, loc);

  if (visible == true) {
    ChSharedPtr<ChBoxShape> box(new ChBoxShape);
    box->GetBoxGeometry().Size = dim;
    //    box->GetBoxGeometry().Pos = loc;		// for Chrono
    box->Pos = loc;  // for Chrono-Parallel
    box->SetColor(ChColor(1, 0, 0));
    box->SetFading(0.6f);
    body->AddAsset(box);
  }
}

// Output

const std::string out_dir = "./SHEAR_DVI";
const std::string pov_dir = out_dir + "/POVRAY";
const std::string shear_file = out_dir + "/shear.dat";

int main(int argc, char* argv[]) {
  // Create output directories

  if (ChFileutils::MakeDirectory(out_dir.c_str()) < 0) {
    cout << "Error creating directory " << out_dir << endl;
    return 1;
  }
  if (ChFileutils::MakeDirectory(pov_dir.c_str()) < 0) {
    cout << "Error creating directory " << pov_dir << endl;
    return 1;
  }

  // Simulation parameters

  double gravity = 9.81;
  double time_step = 1e-3;
  double data_out_step = 1e1 * time_step;
  double visual_out_step = 1e2 * time_step;
  double settling_time = 0.2;
  double begin_shear_time = 10.0;
  double end_simulation_time = 20.0;
  double normal_pressure = 1e3;  // Pa
  double shear_speed = 0.001;    // m/s

  bool write_povray_data = true;

  // Parameters for the balls

  int ballId = 1;  // first ball id

  const int a = 50;
  const int b = 6;
  const int c = 6;
  int numballs = a * b * c;  // number of falling balls = (a X b X c)

  bool dense = true;

  double radius = 0.0025;  // m
  double density = 2500;   // kg/m^3
  double mass = density * (4.0 / 3) * CH_C_PI * radius * radius * radius;
  double Y = 8.0e6;  // Pa
  double nu = 0.3;
  double COR = 0.9;
  double mu = 0.5;

  // Parameters for containing bin, shear box, and load plate

  int groundId = 0;
  int binId = -1;
  int boxId = -2;
  int plateId = -3;
  double width = 0.06;
  double length = 0.06;
  double height = 0.06;
  double thickness = 0.01;

  double shear_Area;
  double shear_Height;
  double shear_Disp;

  ChVector<> pos(0, 0, 0);
  ChQuaternion<> rot(1, 0, 0, 0);
  ChVector<> vel(0, 0, 0);
  real3 force(0, 0, 0);

  // Define two quaternions representing:
  // - a rotation of -90 degrees around x (z2y)
  // - a rotation of +90 degrees around y (z2x)

  ChQuaternion<> z2y;
  ChQuaternion<> z2x;
  z2y.Q_from_AngAxis(-CH_C_PI / 2, ChVector<>(1, 0, 0));
  z2x.Q_from_AngAxis(CH_C_PI / 2, ChVector<>(0, 1, 0));

  // Create the system

  ChSystemParallelDVI* my_system = new ChSystemParallelDVI();
  my_system->GetSettings()->solver.tolerance = 0.01;
  my_system->GetSettings()->solver.max_iteration_bilateral = 100;
  my_system->GetSettings()->solver.clamp_bilaterals = false;
  my_system->GetSettings()->solver.bilateral_clamp_speed = 0.1;

  my_system->GetSettings()->solver.solver_mode = SLIDING;
  my_system->GetSettings()->solver.max_iteration_normal = 0;
  my_system->GetSettings()->solver.max_iteration_sliding = 1000;
  my_system->GetSettings()->solver.max_iteration_spinning = 0;
  my_system->GetSettings()->solver.alpha = 0;
  my_system->GetSettings()->solver.contact_recovery_speed = 0.1;
  my_system->ChangeSolverType(APGDREF);

  my_system->GetSettings()->collision.collision_envelope = 0.05 * radius;
  my_system->GetSettings()->collision.narrowphase_algorithm = NARROWPHASE_HYBRID_MPR;

  my_system->Set_G_acc(ChVector<>(0, -gravity, 0));

#ifdef CHRONO_PARALLEL_HAS_OPENGL
  // Create the OpenGL visualization

  opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
  gl_window.Initialize(800, 600, "hard-sphere (DVI) direct shear box test", my_system);
  gl_window.SetCamera(ChVector<>(3 * width, 0, 0), ChVector<>(0, 0, 0), ChVector<>(0, 1, 0), radius, radius);
  gl_window.SetRenderMode(opengl::SOLID);
#endif

  // Create a material (will be used by all objects)

  ChSharedPtr<ChMaterialSurface> material;
  material = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
  material->SetRestitution(COR);
  material->SetFriction(mu);

  // Create lower bin

  ChSharedPtr<ChBody> bin(new ChBody(new ChCollisionModelParallel));

  bin->SetIdentifier(binId);
  bin->SetMass(1);
  bin->SetPos(ChVector<>(0, -height / 2, 0));
  bin->SetBodyFixed(true);
  bin->SetCollide(true);

  bin->SetMaterialSurface(material);

  bin->GetCollisionModel()->ClearModel();
  AddWall(bin, ChVector<>(width / 2, thickness / 2, length / 2), ChVector<>(0, 0, 0), true);
  AddWall(bin,
          ChVector<>(thickness / 2, height / 2, length / 2 + thickness),
          ChVector<>(-width / 2 - thickness / 2, 0, 0),
          true);
  AddWall(bin,
          ChVector<>(thickness / 2, height / 2, length / 2 + thickness),
          ChVector<>(width / 2 + thickness / 2, 0, 0),
          false);
  AddWall(bin,
          ChVector<>(width / 2 + thickness, height / 2, thickness / 2),
          ChVector<>(0, 0, -length / 2 - thickness / 2),
          true);
  AddWall(bin,
          ChVector<>(width / 2 + thickness, height / 2, thickness / 2),
          ChVector<>(0, 0, length / 2 + thickness / 2),
          true);
  bin->GetCollisionModel()->SetFamily(1);
  bin->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(2);
  bin->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(3);
  bin->GetCollisionModel()->SetFamilyMaskDoCollisionWithFamily(4);
  bin->GetCollisionModel()->BuildModel();

  my_system->AddBody(bin);

  // Create upper shear box

  ChSharedPtr<ChBody> box(new ChBody(new ChCollisionModelParallel));

  box->SetIdentifier(boxId);
  box->SetMass(1);
  box->SetPos(ChVector<>(0, height / 2 + radius, 0));
  box->SetBodyFixed(true);
  box->SetCollide(true);

  box->SetMaterialSurface(material);

  box->GetCollisionModel()->ClearModel();
  AddWall(box,
          ChVector<>(thickness / 2, height / 2, length / 2 + thickness),
          ChVector<>(-width / 2 - thickness / 2, 0, 0),
          true);
  AddWall(box,
          ChVector<>(thickness / 2, height / 2, length / 2 + thickness),
          ChVector<>(width / 2 + thickness / 2, 0, 0),
          false);
  AddWall(box,
          ChVector<>(width / 2 + thickness, height / 2, thickness / 2),
          ChVector<>(0, 0, -length / 2 - thickness / 2),
          true);
  AddWall(box,
          ChVector<>(width / 2 + thickness, height / 2, thickness / 2),
          ChVector<>(0, 0, length / 2 + thickness / 2),
          true);
  box->GetCollisionModel()->SetFamily(2);
  box->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
  box->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(3);
  box->GetCollisionModel()->SetFamilyMaskDoCollisionWithFamily(4);
  box->GetCollisionModel()->BuildModel();

  my_system->AddBody(box);

  // Create upper load plate

  ChSharedPtr<ChBody> plate(new ChBody(new ChCollisionModelParallel));

  shear_Area = width * length;

  plate->SetIdentifier(binId);
  plate->SetMass(normal_pressure * shear_Area / gravity);
  plate->SetPos(ChVector<>(0, 2.0 * radius * float(a) + thickness, 0));
  plate->SetBodyFixed(true);
  plate->SetCollide(true);

  plate->SetMaterialSurface(material);

  plate->GetCollisionModel()->ClearModel();
  AddWall(plate, ChVector<>(width / 2, thickness / 2, length / 2), ChVector<>(0, 0, 0), true);
  plate->GetCollisionModel()->SetFamily(3);
  plate->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
  plate->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(2);
  plate->GetCollisionModel()->SetFamilyMaskDoCollisionWithFamily(4);
  plate->GetCollisionModel()->BuildModel();

  my_system->AddBody(plate);

  // Create (a X b X c) many falling balls

  int i, j, k;
  double ball_x, ball_y, ball_z;

  for (i = 0; i < a; i++) {
    for (j = 0; j < b; j++) {
      for (k = 0; k < c; k++) {
        ball_y = 2.0 * radius * float(i);

        ball_x = 4.0 * radius * (float(j - b / 2) + 0.5) + 0.99 * radius * (float(rand() % 100) / 50 - 1.0);
        ball_z = 4.0 * radius * (float(k - c / 2) + 0.5) + 0.99 * radius * (float(rand() % 100) / 50 - 1.0);

        ChSharedPtr<ChBody> ball(new ChBody(new ChCollisionModelParallel));

        ball->SetIdentifier(ballId + 6 * 6 * i + 6 * j + k);
        ball->SetMass(mass);
        ball->SetInertiaXX((2.0 / 5.0) * mass * radius * radius * ChVector<>(1, 1, 1));
        ball->SetPos(ChVector<>(ball_x, ball_y, ball_z));
        ball->SetBodyFixed(false);
        ball->SetCollide(true);

        ball->SetMaterialSurface(material);

        ball->GetCollisionModel()->ClearModel();
        ball->GetCollisionModel()->AddSphere(radius);
        ball->GetCollisionModel()->SetFamily(4);
        ball->GetCollisionModel()->BuildModel();

        ChSharedPtr<ChSphereShape> sphere(new ChSphereShape);

        sphere->GetSphereGeometry().rad = radius;
        sphere->SetColor(ChColor(1, 0, 1));
        ball->AddAsset(sphere);

        my_system->AddBody(ball);
      }
    }
  }

  // Create prismatic (translational) joint between load plate and shear box.
  // The translational axis of a prismatic joint is along the Z axis of the
  // specified joint coordinate system.  Here, we apply the 'z2y' rotation to
  // align it with the Y axis of the global reference frame.

  ChSharedPtr<ChLinkLockPrismatic> prismatic_plate_box(new ChLinkLockPrismatic);
  prismatic_plate_box->SetName("prismatic_plate_box");
  prismatic_plate_box->Initialize(plate, box, ChCoordsys<>(ChVector<>(0, 0, 0), z2y));
  my_system->AddLink(prismatic_plate_box);

  // Setup output

  ChStreamOutAsciiFile shearStream(shear_file.c_str());
  shearStream.SetNumFormat("%16.4e");

  // Begin simulation

  bool settling = true;
  bool shearing = false;

  int data_out_frame = 0;
  int visual_out_frame = 0;

  while (my_system->GetChTime() < end_simulation_time) {
    if (my_system->GetChTime() > settling_time && settling == true) {
      if (dense == true)
        material->SetFriction(0.05);
      plate->SetPos(ChVector<>(0, height, 0));
      plate->SetBodyFixed(false);
      settling = false;
    }

    if (my_system->GetChTime() > begin_shear_time && shearing == false) {
      if (dense == true)
        material->SetFriction(mu);
      shear_Height = plate->GetPos().y;
      shear_Disp = bin->GetPos().z;
      shearing = true;
    }

    if (shearing == true) {
      bin->SetPos(ChVector<>(0, -height / 2, -shear_speed * begin_shear_time + shear_speed * my_system->GetChTime()));
      bin->SetRot(QUNIT);
    } else {
      bin->SetPos(ChVector<>(0, -height / 2, 0));
      bin->SetRot(QUNIT);
    }

    //  Do time step
    // Advance dynamics.
#ifdef CHRONO_PARALLEL_HAS_OPENGL
    if (gl_window.Active()) {
		gl_window.DoStepDynamics(time_step);
		gl_window.Render();
    } else
    	break;
#else
    my_system->DoStepDynamics(time_step);
#endif

    //  Output to screen

    if (my_system->GetChTime() >= data_out_frame * data_out_step) {
      my_system->CalculateContactForces();
      force = my_system->GetBodyContactForce(0);

      cout << my_system->GetChTime() << "	" << plate->GetPos().y - bin->GetPos().y << "	" << bin->GetPos().x
           << "	" << bin->GetPos().y << "	" << bin->GetPos().z << "	" << force.x << "	" << force.y
           << "	" << force.z << "\n";

      //  Output to shear data file

      if (shearing == true) {
        shearStream << (bin->GetPos().z - shear_Disp) / (2 * radius) << "	";
        shearStream << -force.z / (shear_Area * normal_pressure) << "	";
        shearStream << (plate->GetPos().y - shear_Height) / (2 * radius) << "\n";

        cout << (bin->GetPos().z - shear_Disp) / (2 * radius) << "	";
        cout << -force.z / (shear_Area * normal_pressure) << "	";
        cout << (plate->GetPos().y - shear_Height) / (2 * radius) << "\n";
      }

      data_out_frame++;
    }

    //  Output to POV-Ray

    if (my_system->GetChTime() >= visual_out_frame * visual_out_step) {
      if (write_povray_data) {
        char filename[100];
        sprintf(filename, "%s/data_%03d.dat", pov_dir.c_str(), visual_out_frame + 1);
        utils::WriteShapesPovray(my_system, filename, false);
      }

      visual_out_frame++;
    }

  }

  return 0;
}

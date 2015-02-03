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
//     - soft-sphere (DEM) direct shear box test
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

#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_parallel/lcp/ChLcpSystemDescriptorParallel.h"

#include "chrono_opengl/ChOpenGLWindow.h"

#include "core/ChFileutils.h"
#include "core/ChStream.h"

#include "chrono_utils/ChUtilsCreators.h"

using namespace chrono;
using namespace chrono::collision;

using std::cout;
using std::flush;
using std::endl;

// Utility for adding (visible or invisible) walls

void AddWall(ChSharedBodyDEMPtr& body, const ChVector<>& dim, const ChVector<>& loc, bool visible) {
  body->GetCollisionModel()->AddBox(dim.x, dim.y, dim.z, loc);

  if (visible == true) {
    ChSharedPtr<ChBoxShape> box(new ChBoxShape);
    box->GetBoxGeometry().Size = dim;
//    box->GetBoxGeometry().Pos = loc;		// for Chrono
    box->Pos = loc;							// for Chrono-Parallel
    box->SetColor(ChColor(1, 0, 0));
    box->SetFading(0.6f);
    body->AddAsset(box);
  }
}

// Output

const std::string shear_file = "shearDEM.dat";

int main(int argc, char* argv[]) {

  // Simulation parameters

  double gravity = 9.81;
  double time_step = 1e-5;
  double data_out_step = 1e3 * time_step;
  double visual_out_step = 1e4 * time_step;
  double settling_time = 1.0;
  double begin_shear_time = 10.0;
  double end_simulation_time = 20.0;
  double normal_pressure = 1e3;		// Pa
  double shear_speed = 0.001;		// m/s

  // Parameters for the balls

  int ballId = 1;
  int numballs = 1800;  // right now 1800 balls is hard-wired
  bool dense = true;
  double radius = 0.0025;  // m
  double density = 2500;   // kg/m^3
  double mass = density * (4.0 / 3) * CH_C_PI * radius * radius * radius;
  double Y = 8.0e6;  // Pa
  double nu = 0.3;
  double COR = 0.9;
  double alpha = -log(COR) * 6.66264 / (3.85238 + log(COR));  // Hu et al. (2011)
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
  ChVector<> force(0, 0, 0);

  // Define two quaternions representing:
  // - a rotation of -90 degrees around x (z2y)
  // - a rotation of +90 degrees around y (z2x)

  ChQuaternion<> z2y;
  ChQuaternion<> z2x;
  z2y.Q_from_AngAxis(-CH_C_PI / 2, ChVector<>(1, 0, 0));
  z2x.Q_from_AngAxis(CH_C_PI / 2, ChVector<>(0, 1, 0));

  // Create the system

  ChSystemParallelDEM* my_system = new ChSystemParallelDEM();
  my_system->GetSettings()->solver.tolerance = 0.01;
  my_system->GetSettings()->solver.max_iteration_bilateral = 100;

  my_system->Set_G_acc(ChVector<>(0, -gravity, 0));

  // Create the OpenGL visualization

  opengl::ChOpenGLWindow &gl_window = opengl::ChOpenGLWindow::getInstance();
  gl_window.Initialize(800, 600, "soft-sphere (DEM) direct shear box test", my_system);
  gl_window.SetCamera(ChVector<>(3*width,0,0), ChVector<>(0,0,0),ChVector<>(0,1,0), 0.01, 0.01);
  gl_window.SetRenderMode(opengl::SOLID);

  // Create a material (will be used by all objects)

  ChSharedPtr<ChMaterialSurfaceDEM> material;
  material = ChSharedPtr<ChMaterialSurfaceDEM>(new ChMaterialSurfaceDEM);
  material->SetYoungModulus(Y);
  material->SetPoissonRatio(nu);
  material->SetDissipationFactor(alpha);
  material->SetFriction(mu);

  // Create lower bin

  ChSharedPtr<ChBodyDEM> bin(new ChBodyDEM(new ChCollisionModelParallel));

  bin->SetIdentifier(binId);
  bin->SetMass(1);
  bin->SetPos(ChVector<>(0, -height / 2, 0));
  bin->SetBodyFixed(false);
  bin->SetCollide(true);

  bin->SetMaterialSurfaceDEM(material);

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

  ChSharedPtr<ChBodyDEM> box(new ChBodyDEM(new ChCollisionModelParallel));

  box->SetIdentifier(boxId);
  box->SetMass(1);
  box->SetPos(ChVector<>(0, height / 2 + radius, 0));
  box->SetBodyFixed(true);
  box->SetCollide(true);

  box->SetMaterialSurfaceDEM(material);

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

  ChSharedPtr<ChBodyDEM> plate(new ChBodyDEM(new ChCollisionModelParallel));

  shear_Area = width * length;

  plate->SetIdentifier(binId);
  plate->SetMass(normal_pressure * shear_Area / gravity);
  plate->SetPos(ChVector<>(0, height, 0));
  plate->SetBodyFixed(true);
  plate->SetCollide(false);

  plate->SetMaterialSurfaceDEM(material);

  plate->GetCollisionModel()->ClearModel();
  AddWall(plate, ChVector<>(width / 2, thickness / 2, length / 2), ChVector<>(0, 0, 0), true);
  plate->GetCollisionModel()->SetFamily(3);
  plate->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
  plate->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(2);
  plate->GetCollisionModel()->SetFamilyMaskDoCollisionWithFamily(4);
  plate->GetCollisionModel()->BuildModel();

  my_system->AddBody(plate);

  // Create 1800 falling balls

  int i, j, k;
  double ball_x, ball_y, ball_z;

  for (i = 0; i < 50; i++) {
    for (j = 0; j < 6; j++) {
      for (k = 0; k < 6; k++) {
        ball_y = 2 * radius * float(i);

        ball_x = -0.025 + 0.01 * float(j) + 0.99 * radius * (float(rand() % 100) / 50 - 1.0);
        ball_z = -0.025 + 0.01 * float(k) + 0.99 * radius * (float(rand() % 100) / 50 - 1.0);

        ChSharedPtr<ChBodyDEM> ball(new ChBodyDEM(new ChCollisionModelParallel));

        ball->SetIdentifier(ballId + 6 * 6 * i + 6 * j + k);
        ball->SetMass(mass);
        ball->SetInertiaXX((2.0 / 5.0) * mass * radius * radius * ChVector<>(1, 1, 1));
        ball->SetPos(ChVector<>(ball_x, ball_y, ball_z));
        ball->SetBodyFixed(false);
        ball->SetCollide(true);

        ball->SetMaterialSurfaceDEM(material);

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

  // Create locked (6 DOF fixed) joint between lower bin and shear box.

  ChSharedPtr<ChLinkLockLock> lock_bin_box(new ChLinkLockLock);
  lock_bin_box->SetName("lock_bin_box");
  lock_bin_box->Initialize(bin, box, ChCoordsys<>(ChVector<>(0, 0, 0), QUNIT));
  my_system->AddLink(lock_bin_box);

  // Setup output

  ChStreamOutAsciiFile shearStream(shear_file.c_str());
  shearStream.SetNumFormat("%16.4e");

  // The soft-real-time cycle

  double time = 0.0;
  double data_out_time = 0.0;
  double visual_out_time = 0.0;

  // Begin simulation

  bool settling = true;
  bool shearing = false;
  while (gl_window.Active() && time < end_simulation_time) {

    gl_window.Render();

    while (time < visual_out_time) {
      while (time < data_out_time) {
        if (time > settling_time && settling == true) {
          if (dense == true)
            material->SetFriction(0.05);
          plate->SetBodyFixed(false);
          plate->SetCollide(true);
          settling = false;
        }

        if (time > begin_shear_time && shearing == false) {
          if (dense == true)
            material->SetFriction(mu);
          shear_Height = plate->GetPos().y;
          shear_Disp = bin->GetPos().z;
          shearing = true;
        }

        if(shearing == true)
        {
        	bin->SetPos(ChVector<>(0, -height / 2, -shear_speed * begin_shear_time + shear_speed*time));
//        	bin->SetPos_dt(ChVector<>(0, 0, shear_speed));
//        	bin->SetPos_dtdt(ChVector<>(0, 0, 0));
        } else
        {
        	bin->SetPos(ChVector<>(0, -height / 2, 0));
//        	bin->SetPos_dt(ChVector<>(0, 0, 0));
//        	bin->SetPos_dtdt(ChVector<>(0, 0, 0));
        }

        my_system->DoStepDynamics(time_step);
        time += time_step;
      }

      cout << time << "	" << plate->GetPos().y - bin->GetPos().y << "	"
    		  << bin->GetPos().x << "	"
    		  << bin->GetPos().y << "	"
    		  << bin->GetPos().z << "	"
    		  << lock_bin_box->Get_react_force().x << "	"
    		  << lock_bin_box->Get_react_force().y << "	"
    		  << lock_bin_box->Get_react_force().z << "\n";

      if (shearing == true) {
        shearStream << (bin->GetPos().z - shear_Disp) / (2 * radius) << "	";
        shearStream << lock_bin_box->Get_react_force().z / (shear_Area * normal_pressure) << "	";
        shearStream << (plate->GetPos().y - shear_Height) / (2 * radius) << "\n";

        cout << (bin->GetPos().z - shear_Disp) / (2 * radius) << "	";
        cout << lock_bin_box->Get_react_force().z / (shear_Area * normal_pressure) << "	";
        cout << (plate->GetPos().y - shear_Height) / (2 * radius) << "\n";
      }

      data_out_time += data_out_step;
    }

    visual_out_time += visual_out_step;
  }

  return 0;
}

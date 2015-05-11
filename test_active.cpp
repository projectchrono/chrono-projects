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
// Author: Radu Serban
// =============================================================================
//
// Simple test for the ability to deactivate a body that exits a predefined
// box (axis-aligned).
//
// =============================================================================

#include <stdio.h>
#include <vector>
#include <cmath>

#include "core/ChFileutils.h"
#include "core/ChStream.h"

#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_parallel/lcp/ChLcpSystemDescriptorParallel.h"

#include "chrono_utils/ChUtilsCreators.h"
#include "chrono_utils/ChUtilsInputOutput.h"

// Control use of OpenGL run-time rendering

//#undef CHRONO_PARALLEL_HAS_OPENGL

#ifdef CHRONO_PARALLEL_HAS_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif

using namespace chrono;
using namespace chrono::collision;


int main(int argc, char* argv[]) {
#ifndef CHRONO_PARALLEL_HAS_OPENGL
  std::cout << "This test requires OpenGL support." << std::endl;
  return 1;
#endif

  // Create system
  ChSystemParallelDEM* system = new ChSystemParallelDEM();

  // Set number of threads.
  int threads = 1;
  int max_threads = omp_get_num_procs();
  if (threads > max_threads)
    threads = max_threads;
  system->SetParallelThreadNumber(threads);
  omp_set_num_threads(threads);

  // Set gravitational acceleration
  system->Set_G_acc(ChVector<>(0, 0, 0));

  // Edit system settings
  system->GetSettings()->solver.tolerance = 1e-3;
  system->GetSettings()->solver.max_iteration_bilateral = 20;
  system->GetSettings()->collision.bins_per_axis = I3(10, 10, 10);
  system->GetSettings()->collision.narrowphase_algorithm = NARROWPHASE_R;

  // Specify containing box
  system->GetSettings()->collision.use_aabb_active = true;
  system->GetSettings()->collision.aabb_min = R3(-2, -2, -2);
  system->GetSettings()->collision.aabb_max = R3(2, 2, 2);

  // Create a ground object (just to attach a visual asset)
  ChSharedPtr<ChBody> ground(new ChBody(new ChCollisionModelParallel, ChBody::DEM));
  ground->SetCollide(false);
  ground->SetBodyFixed(true);

  ChSharedPtr<ChBoxShape> box(new ChBoxShape);
  box->GetBoxGeometry().Size = ChVector<>(2, 2, 2);
  box->Pos = ChVector<>(0, 0, 0);
  ground->AddAsset(box);

  system->AddBody(ground);

  // Create a body at origin with some initial velocity
  ChSharedPtr<ChBody> ball(new ChBody(new ChCollisionModelParallel, ChBody::DEM));
  ball->SetCollide(true);
  ball->SetBodyFixed(false);
  ball->SetMass(1);
  ball->SetPos(ChVector<>(0, 0, 0));
  ball->SetPos_dt(ChVector<>(1, 0, 1));

  ball->GetCollisionModel()->ClearModel();
  utils::AddSphereGeometry(ball.get_ptr(), 1);
  ball->GetCollisionModel()->BuildModel();

  system->AddBody(ball);

  // Initialize OpenGL
#ifdef CHRONO_PARALLEL_HAS_OPENGL
  opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
  gl_window.Initialize(1280, 720, "ACTIVE BOX TEST", system);
  gl_window.SetCamera(ChVector<>(0, -20, 0), ChVector<>(0, 0, 0), ChVector<>(0, 0, 1));
  gl_window.SetRenderMode(opengl::WIREFRAME);
#endif
  // Perform the simulation
  double time_step = 1e-3;

  while (true) {
#ifdef CHRONO_PARALLEL_HAS_OPENGL
    if (gl_window.Active()) {
      gl_window.DoStepDynamics(time_step);
      gl_window.Render();
    } else {
      break;
    }
#else
    system->DoStepDynamics(time_step);
#endif

  }

  return 0;
}

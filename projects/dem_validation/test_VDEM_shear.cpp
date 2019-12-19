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
// Author: Daniel Melanz
// =============================================================================
//
// ChronoParallel demo program for testing the filling process.
//
// The global reference frame has Z up.
// All units SI (CGS, i.e., centimeter - gram - second)
//
// =============================================================================

#include <iostream>
#include <vector>
#include <valarray>
#include <string>
#include <sstream>
#include <cmath>

#include "chrono/core/ChStream.h"

#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_parallel/lcp/ChLcpSystemDescriptorParallel.h"

#include "chrono_utils/ChUtilsCreators.h"
#include "chrono_utils/ChUtilsGenerators.h"
#include "chrono_utils/ChUtilsInputOutput.h"

#include "chrono_thirdparty/filesystem/path.h"

// Control use of OpenGL run-time rendering
//#undef CHRONO_PARALLEL_HAS_OPENGL

#ifdef CHRONO_PARALLEL_HAS_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif

#include "demo_utils.h"

using namespace chrono;
using namespace chrono::collision;

using std::cout;
using std::flush;
using std::endl;

// -----------------------------------------------------------------------------
// Global problem definitions
// -----------------------------------------------------------------------------

// Desired number of OpenMP threads (will be clamped to maximum available)
int threads = 20;

// Perform dynamic tuning of number of threads?
bool thread_tuning = false;

// Save PovRay post-processing data?
bool write_povray_data = true;

// Output
std::string out_dir = "../TEST_SHEAR";

std::string pov_dir = out_dir + "/POVRAY";
std::string fill_file = out_dir + "/filling.dat";
std::string stats_file = out_dir + "/stats.dat";
std::string settled_ckpnt_file = out_dir + "/settled.dat";

// Frequency for visualization output
int out_fps = 60;

// Frequency for writing results to output file
int write_fps = 1000;

// Simulation frame at which detailed timing information is printed
int timing_frame = -1;

// =============================================================================
// Create the granular material
//
// Granular material consisting of identical spheres with specified radius and
// material properties; the spheres are generated in a number of vertical
// layers with locations within each layer obtained using Poisson Disk sampling,
// thus ensuring that no two spheres are closer than twice the radius.
// =============================================================================

int CreateGranularMaterial(ChSystem* system, double mu_g, double r_g, double rho_g, double hdimX, double hdimY, double hdimZ)
{
  // -------------------------------------------
  // Create a material for the granular material
  // -------------------------------------------

  ChSharedPtr<ChMaterialSurface> mat_g(new ChMaterialSurface);
  mat_g->SetFriction(mu_g);

  // ---------------------------------------------
  // Create a mixture entirely made out of spheres
  // ---------------------------------------------

  // Create the particle generator with a mixture of 100% spheres
  utils::Generator gen(system);

  utils::MixtureIngredientPtr& m1 = gen.AddMixtureIngredient(utils::MixtureType::SPHERE, 1.0);
  m1->setDefaultMaterialDVI(mat_g);
  m1->setDefaultDensity(rho_g);
  m1->setDefaultSize(r_g);

  // Ensure that all generated particle bodies will have positive IDs.
  gen.setBodyIdentifier(0);

  // ----------------------
  // Generate the particles
  // ----------------------

  double r = 1.01 * r_g;
  ChVector<> hdims(hdimX - r, 0, hdimZ - r);
  ChVector<> center(0, r, 0);

  while (center.y < 2 * hdimY - r)
  {
    gen.createObjectsBox(utils::SamplingType::POISSON_DISK, 2 * r, center, hdims);
    center.y += 2 * r;
  }

  // Return the number of generated particles.
  return gen.getTotalNumBodies();
}

ChSharedPtr<ChBody> createShearPlate(ChSystem* mphysicalSystem, ChVector<> size, double TH, ChVector<> position, ChQuaternion<> rotation, ChColor color)
{
  double L = size.x;
  double H = size.y;
  double W = size.z;

  // create the body
  ChSharedPtr<ChBody> shearPlate(new ChBody(new ChCollisionModelParallel));
  shearPlate->SetPos(position);
  shearPlate->SetRot(rotation);
  mphysicalSystem->Add(shearPlate);

  // Define the collision shape
  shearPlate->GetCollisionModel()->ClearModel();
  utils::AddBoxGeometry(shearPlate.get_ptr(), ChVector<>((L+2*TH)*.5,H*.5,TH*.5), ChVector<>(0,0,.5*W+.5*TH));
  utils::AddBoxGeometry(shearPlate.get_ptr(), ChVector<>((L+2*TH)*.5,H*.5,TH*.5), ChVector<>(0,0,-.5*W-.5*TH));
  utils::AddBoxGeometry(shearPlate.get_ptr(), ChVector<>(TH*.5,H*.5,W*.5), ChVector<>(.5*L+.5*TH,0,0));
  utils::AddBoxGeometry(shearPlate.get_ptr(), ChVector<>(TH*.5,H*.5,W*.5), ChVector<>(-.5*L-.5*TH,0,0));
  shearPlate->GetCollisionModel()->BuildModel();
  shearPlate->SetCollide(true);

  return shearPlate;
}

ChSharedPtr<ChBody> createShearPlateBottom(ChSystem* mphysicalSystem, ChVector<> size, double TH, ChVector<> position, ChQuaternion<> rotation, ChColor color)
{
  double L = size.x;
  double H = size.y;
  double W = size.z;

  // create the body
  ChSharedPtr<ChBody> shearPlate(new ChBody(new ChCollisionModelParallel));
  shearPlate->SetPos(position);
  shearPlate->SetRot(rotation);
  mphysicalSystem->Add(shearPlate);

  // Define the collision shape
  shearPlate->GetCollisionModel()->ClearModel();
  utils::AddBoxGeometry(shearPlate.get_ptr(), ChVector<>((L+2*TH)*.5,H*.5,TH*.5), ChVector<>(0,0,.5*W+.5*TH));
  utils::AddBoxGeometry(shearPlate.get_ptr(), ChVector<>((L+2*TH)*.5,H*.5,TH*.5), ChVector<>(0,0,-.5*W-.5*TH));
  utils::AddBoxGeometry(shearPlate.get_ptr(), ChVector<>(TH*.5,H*.5,W*.5), ChVector<>(.5*L+.5*TH,0,0));
  utils::AddBoxGeometry(shearPlate.get_ptr(), ChVector<>(TH*.5,H*.5,W*.5), ChVector<>(-.5*L-.5*TH,0,0));
  utils::AddBoxGeometry(shearPlate.get_ptr(), ChVector<>((L+2*TH)*.5,H*.5,(W+2*TH)*.5), ChVector<>(0,-H,0));
  shearPlate->GetCollisionModel()->BuildModel();
  shearPlate->SetCollide(true);

  return shearPlate;
}

ChSharedPtr<ChBody> createBox(ChSystem* mphysicalSystem, ChVector<> size, ChVector<> position, ChQuaternion<> rotation, ChColor color)
{
  double L = size.x;
  double H = size.y;
  double W = size.z;

  // create the body
  ChSharedPtr<ChBody> box(new ChBody(new ChCollisionModelParallel));
  box->SetPos(position);
  box->SetRot(rotation);
  mphysicalSystem->Add(box);

  // Define the collision shape
  box->GetCollisionModel()->ClearModel();
  utils::AddBoxGeometry(box.get_ptr(), ChVector<>(L*.5,H*.5,W*.5), ChVector<>(0,0,0));
  box->GetCollisionModel()->BuildModel();
  box->SetCollide(true);

  return box;
}

class ShearBox {
public:
  // THE DATA
  double L;
  double H;
  double HTOP;
  double W;
  double TH;
  double desiredVelocity;
  double normalPressure;
  double mu_particles;
  double particleRadius;

  double sumCielingVelocity;
  double nCielingVelocity;
  double nextFillTime;
  double sphereVolume;

  double percentLoad;
  int stage2Index;

  // bodies
  ChSharedPtr<ChBody> bottom;
  ChSharedPtr<ChBody> top;
  ChSharedPtr<ChBody> cieling;

  ChSystem* sys;

  // THE FUNCTIONS

  ShearBox(ChSystem*  mphysicalSystem, double gravity,  ///< the chrono::engine physical system
    double L, double H, double HTOP, double W, double TH,
    double desiredVelocity, double normalPressure, double mu, double mu_particles, double particleDensity, double r, double dist, int useSpheres, std::string data_folder)
  {
    this->L = L;
    this->H = H;
    this->HTOP = HTOP;
    this->W = W;
    this->TH = TH;
    this->desiredVelocity = desiredVelocity;
    this->normalPressure = normalPressure;
    this->particleRadius = r;
    this->mu_particles = mu_particles;
    sumCielingVelocity = 10;
    nCielingVelocity = 1;
    nextFillTime = 0;
    sphereVolume = 0;
    sys = mphysicalSystem;
    percentLoad = 0;
    stage2Index = 0;

    ChSharedPtr<ChMaterialSurface> mmaterial(new ChMaterialSurface);
    mmaterial->SetFriction(mu_particles); // Friction coefficient of steel

    int ground_coll_fam = 2;
    int bottom_coll_fam = 3;
    int top_coll_fam = 4;
    int cieling_coll_fam = 5;

    double normalLoad = normalPressure*L*W;
    double mass = 1.0;

    // Create bottom
    bottom = createShearPlateBottom(mphysicalSystem, ChVector<>(L,H,W), TH, ChVector<>(0, 0.5*H, 0), ChQuaternion<>(1,0,0,0), ChColor(0.6,0.3,0.6));
    bottom->SetMaterialSurface(mmaterial);
    bottom->SetMass(mass);
    bottom->GetCollisionModel()->SetFamily(bottom_coll_fam);
    bottom->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(ground_coll_fam);
    bottom->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(top_coll_fam);
    bottom->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(cieling_coll_fam);
    bottom->SetBodyFixed(true);
    bottom->SetIdentifier(-2);

    // Create top
    top = createShearPlate(mphysicalSystem, ChVector<>(L,HTOP+TH,W), TH, ChVector<>(0,H+0.5*(HTOP+TH), 0), ChQuaternion<>(1,0,0,0), ChColor(0.3,0.3,0.6));
    top->SetMaterialSurface(mmaterial);
    top->SetMass(mass);
    top->GetCollisionModel()->SetFamily(top_coll_fam);
    top->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(ground_coll_fam);
    top->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(bottom_coll_fam);
    top->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(cieling_coll_fam);
    top->SetBodyFixed(true);
    top->SetIdentifier(-3);

    // Create cieling
    cieling = createBox(mphysicalSystem, ChVector<>(.95*L, TH, .95*W), ChVector<>(0,H+HTOP+0.5*TH,0), ChQuaternion<>(1,0,0,0), ChColor(0.3,0.5,0.3));
    cieling->SetMaterialSurface(mmaterial);
    cieling->SetMass(mass);
    cieling->GetCollisionModel()->SetFamily(cieling_coll_fam);
    cieling->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(ground_coll_fam);
    cieling->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(bottom_coll_fam);
    cieling->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(top_coll_fam);
    cieling->SetIdentifier(-4);

    // create the translational joint between the truss and weight load
    ChSharedPtr<ChLinkLockPrismatic> translationalVert(new ChLinkLockPrismatic);
    translationalVert->Initialize(cieling, top, ChCoordsys<>(top->GetPos(), chrono::Q_from_AngAxis(CH_C_PI/2,VECT_X)) );
    mphysicalSystem->AddLink(translationalVert);

    CreateGranularMaterial(mphysicalSystem, mu_particles, r, particleDensity, 0.5*L, 0.5*(H+HTOP), 0.5*W);
  }

  int applyNormalPressure(double normalPressure)
  {
    cieling->Empty_forces_accumulators();
    cieling->Accumulate_force(ChVector<>(0, -normalPressure*L*W, 0), cieling->GetPos(), false);

    return 0;
  }

  double setParticleDensity(double rho)
  {
    double totalMass = 0;
    std::vector<ChBody*>::iterator abody = sys->Get_bodylist()->begin();
    while (abody != sys->Get_bodylist()->end()) {
      ChBody* bpointer = (*abody);

      double volume = 4.0*CH_C_PI*particleRadius*particleRadius*particleRadius/3.0;
      double mass = rho*volume;
      bpointer->SetMass(mass);
      totalMass+=mass;

      abody++;
    }

    return totalMass;
  }

  int setFriction(double mu)
  {
	ChSharedPtr<ChMaterialSurface> mmaterial(new ChMaterialSurface);
	mmaterial->SetFriction(mu); // Friction coefficient of steel

    std::vector<ChBody*>::iterator abody = sys->Get_bodylist()->begin();
    while (abody != sys->Get_bodylist()->end()) {
      ChBody* bpointer = (*abody);
      bpointer->SetMaterialSurface(mmaterial);
      abody++;
    }
	mmaterial->SetFriction(0); // Friction coefficient of steel
    //ground->SetMaterialSurface(mmaterial);
    top->SetMaterialSurface(mmaterial);
    bottom->SetMaterialSurface(mmaterial);
    cieling->SetMaterialSurface(mmaterial);

    return 0;
  }
};

bool checkFilled(ChSystem* mphysicalSystem, double minTime, double maxTime)
{
  bool filled = false;
  double averageVelocity = 0;
  std::vector<ChBody*>::iterator abody = mphysicalSystem->Get_bodylist()->begin();
  while (abody != mphysicalSystem->Get_bodylist()->end())
  {
    ChVector<> v = (*abody)->GetPos_dt();
    averageVelocity += sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
    abody++;
  }
  averageVelocity = averageVelocity/mphysicalSystem->GetNbodiesTotal();
  printf("  Average velocity: %f\n",averageVelocity);

  if((averageVelocity < 0.1 && mphysicalSystem->GetChTime() > minTime) || mphysicalSystem->GetChTime() > maxTime) filled = true;

  return filled;
}

int fill(ChSystem* mphysicalSystem, ShearBox* shearBox, double desiredBulkDensity, double gravity)
{
  printf("Time: %.4f (Stage 1): Num. Bodies = %d\n", mphysicalSystem->GetChTime(), mphysicalSystem->GetNbodiesTotal());

  // Check to see if box is already filled
  bool filled = checkFilled(mphysicalSystem, 0.1, 0.3);

  int nextStage = 1;
  if(filled) {
    nextStage = 2;
  }

  return nextStage;
}

int compress(ChSystem* mphysicalSystem, ShearBox* shearBox, double muParticles, double gravity)
{
  //printf("Time: %.3f (Stage 2): Ave. vertical speed: %.3f m/s, Ground Pressure: (%.3f, %.3f, %.3f)\n", mphysicalSystem.GetChTime(), shearBox->sumCielingVelocity/shearBox->nCielingVelocity,forceOnGround.x/(L*W),forceOnGround.y/(L*W),forceOnGround.z/(L*W));
  printf("Time: %.4f (Stage 2) Percent loaded: %f\n",mphysicalSystem->GetChTime(), shearBox->percentLoad);

  if(shearBox->stage2Index%5==0) {
    shearBox->cieling->SetMass(shearBox->normalPressure*shearBox->L*shearBox->W*shearBox->percentLoad/gravity);
    //shearBox->cieling->Empty_forces_accumulators();
    //shearBox->cieling->Accumulate_force(ChVector<>(0,-shearBox->normalPressure*shearBox->L*shearBox->W,0)*shearBox->percentLoad,shearBox->cieling->GetPos(),false);
    shearBox->percentLoad+=0.005;
  }
  shearBox->stage2Index++;

  int nextStage = 2;
  if(shearBox->percentLoad>=1.0) {
    shearBox->setFriction(muParticles);
    nextStage = 3;
  }

  return nextStage;
}

int shear(ChSystem* mphysicalSystem, ShearBox* shearBox, double lengthToRun, double hh)
{
  shearBox->bottom->SetPos_dt(ChVector<>(shearBox->desiredVelocity,0,0));
  shearBox->bottom->SetPos(shearBox->bottom->GetPos()+hh*ChVector<>(shearBox->desiredVelocity,0,0));

  printf("Time: %.4f (Stage 3): Speed: %.3f m/s, Pos: (%.7f, %.4f, %.4f)\n", mphysicalSystem->GetChTime(), shearBox->bottom->GetPos_dt().x, shearBox->bottom->GetPos().x, shearBox->bottom->GetPos().y, shearBox->bottom->GetPos().z);

  int nextStage = 3;
  if(fabs(shearBox->bottom->GetPos().x)>fabs(lengthToRun)) nextStage = 4;
  return nextStage;
}

// =============================================================================

int main(int argc, char* argv[])
{
  bool visualize = 1;
  double particleRadius = 0.3; // [cm]
  double disturbance = 0.1; // How how much should I disturb the particles by?
  double normalPressure = 125000; // apply normal force to cieling
  double desiredBulkDensity = 1.3894; // 1.389475 [g/cm^3]
  double muParticles = 0.18;
  double muParticlesRolling = 0;
  double muParticlesSpinning = 0;
  double cohesion = 0;
  double hh = 1e-4; // [s]
  double tolerance = 5000;
  int maxiterations = 5000;
  int useSpheres = 0;
  int useWarmStarting = 1;
  std::string data_folder = "../TEST_SHEAR/";
  pov_dir = data_folder + "/POVRAY";
  stats_file = data_folder + "/stats.dat";
  settled_ckpnt_file = data_folder + "/settled.dat";

  if (argc > 1)
  {
    visualize      = atoi(argv[1]);
    particleRadius = atof(argv[2]);   // [cm]
    normalPressure = atof(argv[3]);   // apply normal force to cieling// 16,888.1 Pa // 44,127.0 Pa// 71,365.9 Pa
    muParticles    = atof(argv[4]);
    hh             = atof(argv[5]);
    tolerance      = atof(argv[6]);
    useSpheres     = atoi(argv[7]);

    std::stringstream dataFolderStream;
    dataFolderStream << "../dst_" << visualize << "_" << particleRadius << "_" << normalPressure << "_" << muParticles << "_" << hh << "_" << tolerance << "_" << useSpheres << "/";
    data_folder = dataFolderStream.str();

    pov_dir = data_folder + "/POVRAY";
    stats_file = data_folder + "/stats.dat";
    settled_ckpnt_file = data_folder + "/settled.dat";
  }

  // Create output directories.
  if (!filesystem::create_directory(filesystem::path(data_folder))) {
    cout << "Error creating directory " << data_folder << endl;
    return 1;
  }

  if (!filesystem::create_directory(filesystem::path(pov_dir))) {
    cout << "Error creating directory " << pov_dir << endl;
    return 1;
  }

  // working in [cm, g, s]
  double gravity = -981; // [cm/s^2];
  double L = 6; // [cm]
  double H = 3; // [cm]
  double HTOP = 10; // [cm]
  double W = 6; // [cm]
  double TH = 3; // [cm]
  double desiredVelocity = 1;//0.166; // [cm/s]
  double particleDensity = 2.55;//atof(argv[4]); // [g/cm^3]
  double lengthToRun = 0.5; // [cm]
  double muWalls = 0;
  bool debug = false;

  // Initialize the shear stress vs. displacement output file
  std::stringstream shearStream;
  shearStream << data_folder << "/dst_" << visualize << "_" << particleRadius << "_" << normalPressure << "_" << muParticles << "_" << hh << "_" << tolerance << "_" << useSpheres << ".csv";
  ChStreamOutAsciiFile shearData(shearStream.str().c_str());

  // -------------
  // Create system
  // -------------

  //cout << "Create DVI system" << endl;
  ChSystemParallelDVI* msystem = new ChSystemParallelDVI();
  msystem->Set_G_acc(ChVector<>(0, gravity, 0));

  // Create the shear box (with no shearing motion and the cieling locked)
  ShearBox* shearBox = new ShearBox(msystem, gravity, L, H, HTOP, W, TH, desiredVelocity, normalPressure, muWalls, muParticles, particleDensity, particleRadius, disturbance, useSpheres, data_folder);
  shearBox->setFriction(0);

  // Set number of threads.
  int max_threads = omp_get_num_procs();
  if (threads > max_threads) threads = max_threads;
  msystem->SetParallelThreadNumber(threads);
  omp_set_num_threads(threads);
  cout << "Using " << threads << " threads" << endl;

  msystem->GetSettings()->max_threads = threads;
  msystem->GetSettings()->perform_thread_tuning = thread_tuning;

  // Edit system settings
  msystem->GetSettings()->solver.tolerance = tolerance;
  msystem->GetSettings()->solver.max_iteration_bilateral = 0;
  msystem->GetSettings()->solver.clamp_bilaterals = false;
  msystem->GetSettings()->solver.bilateral_clamp_speed = 10e30;
  msystem->GetSettings()->solver.solver_mode = SLIDING;
  msystem->GetSettings()->solver.max_iteration_normal = 0;
  msystem->GetSettings()->solver.max_iteration_sliding = maxiterations;
  msystem->GetSettings()->solver.max_iteration_spinning = 0;
  msystem->GetSettings()->solver.alpha = 0;
  msystem->GetSettings()->solver.contact_recovery_speed = 10e30;
  msystem->SetMaxPenetrationRecoverySpeed(10e30);
  msystem->ChangeSolverType(APGDREF);

  msystem->GetSettings()->collision.collision_envelope = 0.05 * particleRadius;

  msystem->GetSettings()->collision.bins_per_axis = I3(10, 10, 10);
  msystem->GetSettings()->collision.narrowphase_algorithm = NARROWPHASE_HYBRID_MPR;

  // --------------
  // Problem set up
  // --------------

  // Depending on problem type:
  // - Select end simulation time
  // - Select output FPS
  // - Create / load objects

  // ----------------------
  // Perform the simulation
  // ----------------------

  // Set number of simulation steps and steps between successive output
  int out_steps = (int) std::ceil((1.0 / hh) / out_fps);
  int write_steps = (int) std::ceil((1.0 / hh) / write_fps);

  // Initialize counters
  double time = 0;
  int sim_frame = 0;
  int out_frame = 0;
  int next_out_frame = 0;
  double exec_time = 0;
  int num_contacts = 0;
  int currentStage = 1;

  // Create output files
  ChStreamOutAsciiFile statsStream(stats_file.c_str());
  shearData.SetNumFormat("%16.4e");
  statsStream.SetNumFormat("%16.4e");

#ifdef CHRONO_PARALLEL_HAS_OPENGL
  opengl::ChOpenGLWindow &gl_window = opengl::ChOpenGLWindow::getInstance();
  gl_window.Initialize(1280, 720, "Shearing Test", msystem);
  gl_window.SetCamera(ChVector<>(0,3,-20), ChVector<>(0,3,0),ChVector<>(0,1,0));
  gl_window.SetRenderMode(opengl::WIREFRAME);
#endif

  // Loop until reaching the end time...
  while (currentStage<=3) {

    // Switch between the direct shear stages
    switch(currentStage)
    {
    case 1:
      // Stage 1: Fill
      currentStage = fill(msystem,shearBox,desiredBulkDensity,gravity);
      break;
    case 2:
      // Stage 2: Compress
      currentStage = compress(msystem,shearBox,muParticles, abs(gravity));
      break;
    case 3:
      // Stage 3: Shear
      msystem->GetSettings()->solver.max_iteration_sliding = 50000;
      currentStage = shear(msystem,shearBox,lengthToRun,hh);
      break;
    }

    // If at an output frame, write PovRay file and print info
    if (sim_frame == next_out_frame) {
      cout << "------------ Output frame:   " << out_frame + 1 << endl;
      cout << "             Sim frame:      " << sim_frame << endl;
      cout << "             Time:           " << time << endl;
      cout << "             Execution time: " << exec_time << endl;

      // Save PovRay post-processing data.
      if (write_povray_data) {
        char filename[100];
        sprintf(filename, "%s/data_%03d.dat", pov_dir.c_str(), out_frame + 1);
        utils::WriteShapesPovray(msystem, filename, false);
      }

      // Create a checkpoint from the current state.
      cout << "             Write checkpoint data " << flush;
      utils::WriteCheckpoint(msystem, settled_ckpnt_file);
      cout << msystem->Get_bodylist()->size() << " bodies" << endl;

      // Increment counters
      out_frame++;
      next_out_frame += out_steps;
    }

    // Advance simulation by one step
#ifdef CHRONO_PARALLEL_HAS_OPENGL
    if (gl_window.Active()) {
      gl_window.DoStepDynamics(hh);
      gl_window.Render();
    }
#else
    msystem->DoStepDynamics(hh);
#endif
    std::vector<double> history = ((ChLcpIterativeSolver*) (msystem->GetLcpSolverSpeed()))->GetViolationHistory();
    int numIters = history.size();
    double residual = 0;
    if(numIters) residual = history[numIters-1];
    cout << "  Residual: " << residual << ", Iterations: " << numIters << endl;

    // Record stats about the simulation
    if (sim_frame % write_steps == 0) {
      // Compute contact force on container (container should be the first body)
      real3 force(0, 0, 0);
      double shearForce = 0;
      if(currentStage == 3 && msystem->GetNcontacts())
      {
        msystem->CalculateContactForces();
    	force = msystem->GetBodyContactForce(0);
        shearForce = force.x;
        cout << "  Shear force: " << shearForce << endl;
      }

      // write fill info
      shearData << time << ", " << numIters << ", " << residual << ", " << shearForce << ", " << shearBox->bottom->GetPos().x << ", " << shearBox->top->GetPos().y << ", " << shearBox->top->GetPos().z << ", " << shearBox->top->GetPos_dt().x << ", " << shearBox->top->GetPos_dt().y << ", " << shearBox->top->GetPos_dt().z << ", \n";
      shearData.GetFstream().flush();

      // write stat info
      statsStream << time << ", " << exec_time << ", " << num_contacts/write_steps << ", " << numIters << ", " << residual << ", \n";
      statsStream.GetFstream().flush();

      num_contacts = 0;
    }

    // Increment counters
    time += hh;
    sim_frame++;
    exec_time += msystem->GetTimerStep();
    num_contacts += msystem->GetNcontacts();

    // If requested, output detailed timing information for this step
    if (sim_frame == timing_frame)
      msystem->PrintStepStats();
  }

  // ----------------
  // Final processing
  // ----------------

  // Create a checkpoint from the last state
  cout << "             Write checkpoint data " << flush;
  utils::WriteCheckpoint(msystem, settled_ckpnt_file);
  cout << msystem->Get_bodylist()->size() << " bodies" << endl;

  // Final stats
  cout << "==================================" << endl;
  cout << "Number of bodies:  " << msystem->Get_bodylist()->size() << endl;
  cout << "Simulation time:   " << exec_time << endl;
  cout << "Number of threads: " << threads << endl;

  return 0;
}

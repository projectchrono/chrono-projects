#include <stdio.h>
#include <vector>
#include <cmath>

#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_parallel/lcp/ChLcpSystemDescriptorParallel.h"

#include "chrono_utils/ChUtilsInputOutput.h"
#include "chrono_utils/ChUtilsSamplers.h"
#include "chrono_utils/ChUtilsGenerators.h"

using namespace chrono;
using namespace geometry;

// =======================================================================
// Global problem definitions
const std::string data_folder = "../TEST/";

void writePoints(const std::string filename, const utils::PointVectorD& points) {
  utils::CSV_writer csv(" ");
  for (int i = 0; i < points.size(); i++) {
    csv << points[i] << std::endl;
  }
  char file[100];
  csv.write_to_file(data_folder + filename);
  std::cout << "Wrote " << points.size() << " points to " << data_folder + filename << std::endl;
}

void checkSamplers() {
  utils::PointVectorD points;

  // Grid sampler
  utils::GridSampler<> gs(ChVector<>(0.1, 0.2, 0.3));

  points = gs.SampleBox(ChVector<>(2, 0, 0), ChVector<>(0, 2, 1));
  writePoints("GS_rectangleX.out", points);

  points = gs.SampleBox(ChVector<>(0, 3, 0), ChVector<>(1, 0, 2));
  writePoints("GS_rectangleY.out", points);

  points = gs.SampleBox(ChVector<>(0, 0, 2), ChVector<>(1, 2, 0));
  writePoints("GS_rectangleZ.out", points);

  points = gs.SampleBox(ChVector<>(-1, -2, -3), ChVector<>(1, 1.5, 0.5));
  writePoints("GS_box.out", points);

  points = gs.SampleSphere(ChVector<>(1, 1, 1), 1);
  writePoints("GS_sphere.out", points);

  points = gs.SampleCylinderX(ChVector<>(1, 0, -1), 2, 0);
  writePoints("GS_circleX.out", points);

  points = gs.SampleCylinderY(ChVector<>(1, 0, -1), 2, 0);
  writePoints("GS_circleY.out", points);

  points = gs.SampleCylinderZ(ChVector<>(1, 0, -1), 2, 0);
  writePoints("GS_circleZ.out", points);

  points = gs.SampleCylinderZ(ChVector<>(1, 0, -1), 2, 0.5);
  writePoints("GS_cylinderZ.out", points);

  // Poisson Disk sampler
  utils::PDSampler<> pd(0.2);

  points = pd.SampleBox(ChVector<>(2, 0, 0), ChVector<>(0, 2, 1));
  writePoints("PD_rectangleX.out", points);

  points = pd.SampleBox(ChVector<>(0, 3, 0), ChVector<>(1, 0, 2));
  writePoints("PD_rectangleY.out", points);

  points = pd.SampleBox(ChVector<>(0, 0, 2), ChVector<>(1, 2, 0));
  writePoints("PD_rectangleZ.out", points);

  points = pd.SampleCylinderX(ChVector<>(1, 0, -1), 2, 0);
  writePoints("PD_circleX.out", points);

  points = pd.SampleCylinderY(ChVector<>(1, 0, -1), 2, 0);
  writePoints("PD_circleY.out", points);

  points = pd.SampleCylinderZ(ChVector<>(1, 0, -1), 2, 0);
  writePoints("PD_circleZ.out", points);

  points = pd.SampleBox(ChVector<>(-1, -2, -3), ChVector<>(1, 1.5, 1));
  writePoints("PD_box.out", points);

  points = pd.SampleSphere(ChVector<>(1, 1, 1), 1);
  writePoints("PD_sphere.out", points);

  points = pd.SampleCylinderZ(ChVector<>(1, 0, -1), 2, 0.5);
  writePoints("PD_cylinderZ.out", points);
}

void checkGenerators() {
  char filename[100];

  // Parameters for the falling ball
  int ballId = 100;
  double radius = .1;
  double density = 2000;
  double volume = (4.0 / 3) * CH_PI * radius * radius * radius;
  double mass = density * volume;
  ChVector<> inertia = 0.4 * mass * radius * radius * ChVector<>(1, 1, 1);
  ChVector<> init_vel(0, 0, 0);

  // Create system
  ChSystemParallelDEM* msystem = new ChSystemParallelDEM();

  // Create a material for the balls
  ChSharedPtr<ChMaterialSurfaceDEM> ballMat;
  ballMat = ChSharedPtr<ChMaterialSurfaceDEM>(new ChMaterialSurfaceDEM);
  ballMat->SetYoungModulus(5e4f);
  ballMat->SetFriction(0.4f);
  ballMat->SetRestitution(0.2f);

  // CHECK GENERATOR
  utils::Generator gen(msystem);

  utils::MixtureIngredientPtr& m1 = gen.AddMixtureIngredient(utils::MixtureType::SPHERE, 0.1);
  m1->setDefaultSize(radius);
  m1->setDistributionDensity(1000, 600, 500, 3000);

  utils::MixtureIngredientPtr& m2 = gen.AddMixtureIngredient(utils::MixtureType::BOX, 0.25);
  m2->setDefaultSize(radius);
  m2->setDefaultDensity(1000);

  utils::MixtureIngredientPtr& m3 = gen.AddMixtureIngredient(utils::MixtureType::ELLIPSOID, 0.15);
  m3->setDefaultSize(radius);
  m3->setDefaultDensity(2000);

  gen.createObjectsBox(utils::SamplingType::POISSON_DISK, radius, ChVector<>(0, 0, 2), ChVector<>(3, 3, 0));

  std::cout << "Total mass: " << gen.getTotalMass() << "  Total volume: " << gen.getTotalVolume() << std::endl;

  gen.writeObjectInfo(data_folder + "body_info.dat");
}

int main(int argc, char* argv[]) {
  checkSamplers();
  // checkGenerators();

  return 0;
}

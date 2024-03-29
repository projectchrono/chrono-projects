#include "chrono_distributed/physics/ChSystemDistributed.h"
#include "chrono_distributed/collision/ChCollisionModelDistributed.h"
#include "chrono/physics/ChBody.h"

#include <mpi.h>
#include <memory>

using namespace chrono;

double dt = 0.001;

// To be run on 2 MPI ranks
int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    ChSystemDistributed sys(MPI_COMM_WORLD, 1.0, 10000);
    sys.GetDomain()->SetSimDomain(ChVector<>(0, 0, 0), ChVector<>(10, 10, 20));

    sys.Set_G_acc(ChVector<double>(0, 0, -9.8));

    auto ball1 = chrono_types::make_shared<ChBody>(chrono_types::make_shared<ChCollisionModelDistributed>());
    auto ball2 = chrono_types::make_shared<ChBody>(chrono_types::make_shared<ChCollisionModelDistributed>());

    auto material = chrono_types::make_shared<ChMaterialSurfaceSMC>();

    ChVector<double> pos1(5, 5, 0.01);
    ball1->SetPos(pos1);
    ball1->GetCollisionModel()->Clear();
    auto ct_shape1 = chrono_types::make_shared<ChCollisionShapeSphere>(material, 0.15);
    ball1->GetCollisionModel()->AddShape(ct_shape1, ChFrame<>(pos1, QUNIT));
    ball1->GetCollisionModel()->Build();
    ball1->SetCollide(true);

    ChVector<double> pos2(5, 5, 11.02);
    ball2->SetPos(pos2);
    ball2->GetCollisionModel()->Clear();
    auto ct_shape2 = chrono_types::make_shared<ChCollisionShapeSphere>(material, 0.15);
    ball2->GetCollisionModel()->AddShape(ct_shape2, ChFrame<>(pos2, QUNIT));
    ball2->GetCollisionModel()->Build();
    ball2->SetCollide(true);

    sys.GetDomain()->PrintDomain();
    sys.AddBody(ball1);
    sys.AddBody(ball2);

    for (int i = 0; i < 10000; i++) {
        if (i % 10 == 0) {
            // sys.PrintBodyStatus();
        }
        sys.DoStepDynamics(dt);
    }

    MPI_Finalize();
    return 0;
}

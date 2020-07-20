#pragma once
#include "chrono_models/vehicle/m113/M113_Chassis.h"
#include "chrono_models/vehicle/m113/M113_SimpleDriveline.h"
#include "chrono_models/vehicle/m113/M113_SimplePowertrain.h"
#include "chrono_models/vehicle/m113/M113_TrackAssemblyDoublePin.h"
#include "chrono_models/vehicle/m113/M113_TrackAssemblySinglePin.h"
#include "chrono_models/vehicle/m113/M113_Vehicle.h"
using namespace chrono::vehicle;
using namespace chrono::vehicle::m113;

class M113_Chassis_Custom : public M113_Chassis {
  public:
    M113_Chassis_Custom(const std::string& name, bool fixed) : M113_Chassis(name, fixed) {}

    virtual void Initialize(ChSystem* system,                ///< [in] containing system
                            const ChCoordsys<>& chassisPos,  ///< [in] absolute chassis position
                            double chassisFwdVel,            ///< [in] initial chassis forward velocity
                            int collision_family = 0         ///< [in] collision family
    ) override {
        m_body = std::shared_ptr<ChBodyAuxRef>(system->NewBodyAuxRef());
        m_body->SetIdentifier(0);
        m_body->SetName("chassis");
        m_body->SetMass(GetMass());
        m_body->SetFrame_COG_to_REF(ChFrame<>(GetLocalPosCOM(), ChQuaternion<>(1, 0, 0, 0)));
        m_body->SetInertia(GetInertia());

        m_body->SetFrame_REF_to_abs(ChFrame<>(chassisPos));
        m_body->SetPos_dt(chassisFwdVel * chassisPos.TransformDirectionLocalToParent(ChVector<>(1, 0, 0)));

        auto chassis_mesh = std::make_shared<chrono::geometry::ChTriangleMeshConnected>();
        std::vector<std::vector<ChVector<double> > > chassis_hulls;
        utils::LoadConvexHulls(vehicle::GetDataFile("M113/Chassis_Hulls.obj"), *chassis_mesh, chassis_hulls);

        auto mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();

        m_body->GetAssets().clear();
        m_body->GetCollisionModel()->ClearModel();
        utils::AddConvexCollisionModel(m_body, mat, chassis_mesh, chassis_hulls, VNULL, QUNIT);

        m_body->GetCollisionModel()->BuildModel();

        m_body->GetCollisionModel()->SetFamily(4);
        m_body->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(4);

        m_body->SetCollide(true);

        system->Add(m_body);

        m_body->SetBodyFixed(m_fixed);
    }
};

class M113_Vehicle_Custom : public ChTrackedVehicle {
  public:
    M113_Vehicle_Custom(bool fixed, TrackShoeType shoe_type)
        : ChTrackedVehicle("M113 Vehicle", ChContactMethod::NSC), m_type(shoe_type) {
        Create(fixed);
    }

    M113_Vehicle_Custom(bool fixed, TrackShoeType shoe_type, ChSystem* system)
        : ChTrackedVehicle("M113 Vehicle", system), m_type(shoe_type) {
        assert(system->GetContactMethod() == ChContactMethod::NSC);
        Create(fixed);
    }

    ~M113_Vehicle_Custom() {}

    virtual void Initialize(const ChCoordsys<>& chassisPos, double chassisFwdVel = 0) override {
        m_chassis->Initialize(m_system, chassisPos, chassisFwdVel);

        // Initialize the left and right track assemblies.
        double track_offset = 1.0795;

        if (ChSystemParallelNSC* system_dvi = dynamic_cast<ChSystemParallelNSC*>(m_system)) {
            printf("BodiesA: %d\n", system_dvi->data_manager->num_rigid_shapes);
            m_tracks[0]->Initialize(m_chassis, ChVector<>(0, track_offset, 0));

            printf("BodiesB: %d\n", system_dvi->data_manager->num_rigid_shapes);

            m_tracks[1]->Initialize(m_chassis, ChVector<>(0, -track_offset, 0));

            printf("BodiesC: %d\n", system_dvi->data_manager->num_rigid_shapes);

            for (int i = 25; i < 558; i++) {
                system_dvi->data_manager->shape_data.fam_rigid[i].x = 1 << 4;
            }
        }

        // Initialize the driveline subsystem
        m_driveline->Initialize(m_chassis, m_tracks[0], m_tracks[1]);
    }

  protected:
    void Create(bool fixed) {
        // Create the chassis subsystem
        m_chassis = std::make_shared<M113_Chassis_Custom>("Chassis", fixed);

        // Create the track assembly subsystems
        switch (m_type) {
            case TrackShoeType::SINGLE_PIN:
                m_tracks[0] = std::make_shared<M113_TrackAssemblySinglePin>(LEFT, BrakeType::SIMPLE);
                m_tracks[1] = std::make_shared<M113_TrackAssemblySinglePin>(RIGHT, BrakeType::SIMPLE);
                break;
            case TrackShoeType::DOUBLE_PIN:
                m_tracks[0] = std::make_shared<M113_TrackAssemblyDoublePin>(LEFT, BrakeType::SIMPLE);
                m_tracks[1] = std::make_shared<M113_TrackAssemblyDoublePin>(RIGHT, BrakeType::SIMPLE);
                break;
        }

        // Create the driveline
        m_driveline = std::make_shared<M113_SimpleDriveline>();
        ////m_driveline = std::make_shared<M113_DrivelineBDS>();

        GetLog() << "M113 custom vehicle mass = " << GetVehicleMass() << " kg.\n";
    }

    TrackShoeType m_type;  ///< type of track assembly (SINGLE_PIN or DOUBLE_PIN)
};

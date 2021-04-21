// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Radu Serban, Justin Madsen, Conlain Kelly, Aaron Young
// =============================================================================
//
// Interactive driver for a vehicle. This class implements the
// functionality required by its base ChDriver class using keyboard or joystick
// inputs. If a joystick is present it will use that as an input; it will
// otherwise default to a keyboard input.
//
// =============================================================================

#ifndef CH_LIDAR_WAYPOINT_DRIVER_H
#define CH_LIDAR_WAYPOINT_DRIVER_H

#include <string>

#include "chrono_sensor/ChLidarSensor.h"
#include "chrono_vehicle/ChApiVehicle.h"
#include "chrono_vehicle/ChDriver.h"
#include "chrono_vehicle/driver/ChDataDriver.h"
#include "chrono_vehicle/driver/ChPathFollowerACCDriver.h"

using namespace chrono::vehicle;
using namespace chrono::sensor;

namespace chrono {
namespace synchrono {

class CH_VEHICLE_API ChLidarWaypointDriver : public ChDriver {
  public:
    /// Construct an interactive driver.
    ChLidarWaypointDriver(ChVehicle& vehicle,
                          std::shared_ptr<ChLidarSensor> lidar,
                          std::shared_ptr<ChBezierCurve> path,
                          const std::string& path_name,
                          double target_speed,           ///< constant target speed
                          double target_following_time,  ///< seconds of following time
                          double target_min_distance,    ///< min following distance
                          double current_distance,       ///< current distance to the vehicle in front
                          bool isClosedPath              ///< Treat the path as a closed loop
    );

    virtual ~ChLidarWaypointDriver() {}

    /// Update the state of this driver system at the specified time.
    virtual void Synchronize(double time) override;

    /// Advance the state of this driver system by the specified time step.
    virtual void Advance(double step) override;

    void SetCurrentDistance(double dist) {
        m_current_distance = std::min(m_current_distance, dist);
        next_dist_reset_time += 0.1;
    }

    /// Set gains for internal dynamics.
    void SetGains(double lookahead,
                  double p_steer,
                  double i_steer,
                  double d_steer,
                  double p_acc,
                  double i_acc,
                  double d_acc);

  private:
    void MinDistFromLidar();

    std::shared_ptr<ChLidarSensor> m_lidar;
    double m_current_distance;
    double m_last_lidar_time = 0;
    double next_dist_reset_time = 0;
    double m_target_speed;
    std::shared_ptr<ChBezierCurve> m_path;

    std::shared_ptr<ChPathFollowerACCDriver> m_acc_driver;  ///< underlying acc driver
};

}  // namespace synchrono
}  // namespace chrono

#endif

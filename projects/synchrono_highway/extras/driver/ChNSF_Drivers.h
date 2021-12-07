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
// Authors: Simone Benatti
// =============================================================================
//
// Custom drivers for the NSF project. 
// Both are specialization of ChPathFollowerDriver
// The leader will adjust its behavior depending on the traveled distance
// The follower will adjust the speed to reach a target gap with the leader
//
// =============================================================================

#ifndef CH_NSF_DRIVER_H
#define CH_NSF_DRIVER_H

#include <string>

#include "chrono_sensor/sensors/ChLidarSensor.h"
#include "chrono_vehicle/ChApiVehicle.h"
#include "chrono_vehicle/ChDriver.h"
#include "chrono_vehicle/driver/ChDataDriver.h"
#include "chrono_vehicle/driver/ChPathFollowerACCDriver.h"
#include "chrono_vehicle/driver/ChPathFollowerDriver.h"

using namespace chrono::vehicle;
using namespace chrono::sensor;

#define MS_TO_MPH 2.23694
#define MPH_TO_MS 0.44704
#define AUDI_LENGTH 4.86

namespace chrono {
namespace synchrono {

class CH_VEHICLE_API ChNSFLeaderDriver : public ChPathFollowerDriver {
  public:
    /// Construct an interactive driver.
    ChNSFLeaderDriver(ChVehicle& vehicle,                                      ///< associated vehicle
                         const std::string& steering_filename,                 ///< JSON file with steering controller specification
                         const std::string& speed_filename,                    ///< JSON file with speed controller specification
                         std::shared_ptr<ChBezierCurve> path,                  ///< Bezier curve with target path
                         const std::string& path_name,                         ///< name of the path curve
                         double target_speed,                                  ///< constant target speed
                         std::vector<std::vector<double>> behavior,            ///< piecewise directives
                         bool isClosedPath = false                             ///< Treat the path as a closed loop
    );

    virtual ~ChNSFLeaderDriver() {}

    void Synchronize(double time);

  private:
    // starting pos to compare with to obtain traveled dist 
    ChVector<> startingPos;
    //JSON file containing the instruction for target speed
     std::vector<std::vector<double>> behavior_data;
};

}  // namespace synchrono
}  // namespace chrono

#endif

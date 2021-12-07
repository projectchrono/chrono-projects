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

#include "ChNSF_Drivers.h"
#include "chrono_vehicle/wheeled_vehicle/vehicle/WheeledVehicle.h"

#include <algorithm>
#include <climits>
#include <fstream>
#include <iostream>
#include <sstream>

namespace chrono {
namespace synchrono {

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
ChNSFLeaderDriver::ChNSFLeaderDriver(ChVehicle& vehicle,             ///< associated vehicle
                         const std::string& steering_filename,       ///< JSON file with steering controller specification
                         const std::string& speed_filename,          ///< JSON file with speed controller specification
                         std::shared_ptr<ChBezierCurve> path,        ///< Bezier curve with target path
                         const std::string& path_name,               ///< name of the path curve
                         double target_speed,                        ///< constant target speed
                         std::vector<std::vector<double>> behavior,  ///< JSON file with piecewise directives
                         bool isClosedPath)                          ///< Treat the path as a closed loop
    : ChPathFollowerDriver(vehicle, steering_filename, speed_filename,  path, path_name, target_speed, isClosedPath), behavior_data(behavior) {
        startingPos = vehicle.GetChassis()->GetPos();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void ChNSFLeaderDriver::Synchronize(double time) {
    double dist = (m_vehicle.GetChassis()->GetPos() - startingPos).Length();

    ChPathFollowerDriver::Synchronize(time);
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
/*void ChLidarWaypointDriver::Advance(double step) {
    // calculate a new current 
    int a = 0;
}
*/
}  // namespace synchrono
}  // end namespace chrono

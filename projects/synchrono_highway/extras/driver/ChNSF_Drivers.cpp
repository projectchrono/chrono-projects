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

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
ChNSFLeaderDriver::ChNSFLeaderDriver(
    ChVehicle& vehicle,                         ///< associated vehicle
    const std::string& steering_filename,       ///< JSON file with steering controller specification
    const std::string& speed_filename,          ///< JSON file with speed controller specification
    std::shared_ptr<ChBezierCurve> path,        ///< Bezier curve with target path
    const std::string& path_name,               ///< name of the path curve
    double target_speed,                        ///< constant target speed
    std::vector<std::vector<double>> behavior,  ///< JSON file with piecewise directives
    bool isClosedPath)                          ///< Treat the path as a closed loop
    : ChPathFollowerDriver(vehicle, steering_filename, speed_filename, path, path_name, target_speed, isClosedPath),
      behavior_data(behavior),
      cruise_speed(target_speed) {
    previousPos = vehicle.GetChassis()->GetPos();
    dist = 0;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void ChNSFLeaderDriver::Synchronize(double time) {
    // In this portion we adjust the target speed according to custom piece-wise sinusoidal defined in behavior_data:
    // start [miles], end [miles], offset [mph], amplitude [mph], period [miles], phase [miles]
    dist += (m_vehicle.GetChassis()->GetPos() - previousPos).Length() * M_TO_MILE;
    previousPos = m_vehicle.GetChassis()->GetPos();
    SetDesiredSpeed(cruise_speed);
    /*
    for (auto piece_data : behavior_data) {
        // if the traveled dist is > max, we inspect the following piece
        if (dist > piece_data[1])
            continue;
        else {
            // if the new piece has not been reached yet, we keep cruise speed:
            if (dist < piece_data[0])
                SetDesiredSpeed(cruise_speed);
            else {
                double v_mph =
                    piece_data[2] + piece_data[3] * sin(CH_C_2PI * (1 / piece_data[4]) * dist + piece_data[5]);
                // double des_speed = v_mph * MPH_TO_MS ;
                // std::cout<< "Miles traveled:" << dist << "   Desired Speed is :" << des_speed << "\n";
                SetDesiredSpeed(v_mph * MPH_TO_MS);
            }
        }
    }*/

    ChPathFollowerDriver::Synchronize(time);
}

ChNSFFollowererDriver::ChNSFFollowererDriver(
    ChVehicle& vehicle,                       ///< associated vehicle
    const std::string& steering_filename,     ///< JSON file with steering controller specification
    const std::string& speed_filename,        ///< JSON file with speed controller specification
    std::shared_ptr<ChBezierCurve> path,      ///< Bezier curve with target path
    const std::string& path_name,             ///< name of the path curve
    double target_speed,                      ///< constant target speed
    const ChVehicle& lead_vehicle,            ///< followed_vehicle
    std::vector<std::vector<double>> params,  ///< JSON file with piecewise params
    bool isClosedPath)                        ///< Treat the path as a closed loop
    : ChPathFollowerDriver(vehicle, steering_filename, speed_filename, path, path_name, target_speed, isClosedPath),
      behavior_data(params),
      cruise_speed(target_speed),
      leader(lead_vehicle) {
    previousPos = vehicle.GetChassis()->GetPos();
    dist = 0;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void ChNSFFollowererDriver::Synchronize(double time, double step) {
    // In this portion we adjust the target speed according to custom piece-wise sinusoidal defined in behavior_data.
    // We use the driver model explained here, using a desired speed instead:
    // https://traffic-simulation.de/info/info_IDM.html the parameters are: start [miles], end [miles], v0 desired v
    // [m/s], T desired time headway [s], desired space headway [m], a: accel reate a [m/s^2], b: comfort decel [m/s^2],
    // delta: accel exponent
    dist += (m_vehicle.GetChassis()->GetPos() - previousPos).Length() * M_TO_MILE;
    previousPos = m_vehicle.GetChassis()->GetPos();
    for (auto piece_data : behavior_data) {
        // if the traveled dist is > max, we inspect the following piece
        if (dist > piece_data[1])
            continue;
        else {
            // if the new piece has not been reached yet, we keep cruise speed:
            if (dist < piece_data[0])
                SetDesiredSpeed(cruise_speed);
            else {
                double s = (m_vehicle.GetChassis()->GetPos() - leader.GetChassis()->GetPos()).Length() - AUDI_LENGTH;
                double v = m_vehicle.GetChassis()->GetSpeed();
                double delta_v = v - leader.GetChassis()->GetSpeed();
                double s_star =
                    piece_data[4] +
                    ChMax(0.0, v * piece_data[3] + (v * delta_v) / (2 * sqrt(piece_data[5] * piece_data[6])));
                double dv_dt = piece_data[5] * (1 - pow(v / piece_data[2], piece_data[7]) - pow(s_star / s, 2));
                double v_ms = ChMax(0.0, v + dv_dt * step);
                // double des_speed = v_mph * MPH_TO_MS ;
                // std::cout<< "Miles traveled:" << dist << "   Desired Speed is :" << des_speed << "\n";
                SetDesiredSpeed(v_ms);
            }
        }
    }

    ChPathFollowerDriver::Synchronize(time);
}

void ChNSFFollowererDriver::SetCruiseSpeed(double speed) {
    cruise_speed = speed;
}

void ChNSFLeaderDriver::SetCruiseSpeed(double speed) {
    cruise_speed = speed;
}

}  // end namespace chrono

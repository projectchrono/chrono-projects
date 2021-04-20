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

#include "ChLidarWaypointDriver.h"

#include <algorithm>
#include <climits>
#include <fstream>
#include <iostream>
#include <sstream>

namespace chrono {
namespace synchrono {

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
ChLidarWaypointDriver::ChLidarWaypointDriver(ChVehicle& vehicle,
                                             std::shared_ptr<ChLidarSensor> lidar,
                                             std::shared_ptr<ChBezierCurve> path,
                                             const std::string& path_name,
                                             double target_speed,           ///< constant target speed
                                             double target_following_time,  ///< seconds of following time
                                             double target_min_distance,    ///< min following distance
                                             double current_distance,  ///< current distance to the vehicle in front
                                             bool isClosedPath         ///< Treat the path as a closed loop
                                             )
    : ChDriver(vehicle), m_lidar(lidar), m_target_speed(target_speed), m_path(path), m_current_distance(100.0) {
    m_acc_driver = chrono_types::make_shared<ChPathFollowerACCDriver>(vehicle, path, path_name, target_speed,
                                                                      target_following_time, target_min_distance,
                                                                      current_distance, isClosedPath);
    m_acc_driver->GetSpeedController().SetGains(0.5, 0, 0);
    m_acc_driver->GetSteeringController().SetGains(0.5, 0, 0);
    m_acc_driver->GetSteeringController().SetLookAheadDistance(8.0);
    m_acc_driver->Initialize();
    m_acc_driver->Reset();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

void ChLidarWaypointDriver::SetGains(double lookahead,
                                     double p_steer,
                                     double i_steer,
                                     double d_steer,
                                     double p_acc,
                                     double i_acc,
                                     double d_acc) {
    m_acc_driver->GetSpeedController().SetGains(p_acc, i_acc, d_acc);
    m_acc_driver->GetSteeringController().SetGains(p_steer, i_steer, d_steer);
    m_acc_driver->GetSteeringController().SetLookAheadDistance(lookahead);
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void ChLidarWaypointDriver::Synchronize(double time) {
    if (time > next_dist_reset_time && !m_lidar) {
        next_dist_reset_time = time + 0.01;
        m_current_distance = 100;
    }

    m_acc_driver->Synchronize(time);
}

void ChLidarWaypointDriver::MinDistFromLidar() {
    // only do this at 20 Hz
    double update_period = 0.05;
    double min_distance = 100;

    // box of interest
    double x_min = .1;
    double x_max = 100;
    double y_min = -2.0;
    double y_max = 2.0;
    double z_min = -.1;
    double z_max = 1.5;

    UserXYZIBufferPtr xyzi_buffer = m_lidar->GetMostRecentBuffer<UserXYZIBufferPtr>();
    if (xyzi_buffer->TimeStamp > m_last_lidar_time + update_period) {
        m_last_lidar_time = xyzi_buffer->TimeStamp;
        for (int i = 0; i < xyzi_buffer->Height; i++) {
            for (int j = 0; j < xyzi_buffer->Width; j++) {
                PixelXYZI xyzi = xyzi_buffer->Buffer[i * xyzi_buffer->Width + j];
                // intensity threshold
                if (xyzi.intensity > 0) {
                    // x threshold
                    if (xyzi.x < std::min(x_max, min_distance) && xyzi.x > x_min) {
                        if (xyzi.y < y_max && xyzi.y > y_min) {
                            if (xyzi.z < z_max && xyzi.z > z_min) {
                                min_distance = xyzi.x;
                            }
                        }
                    }
                }
            }
        }
        m_current_distance = min_distance;
    }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void ChLidarWaypointDriver::Advance(double step) {
    // calculate a new current distance using the lidar
    if (m_lidar) {
        MinDistFromLidar();
    }

    // calculate a new current speed using the path curvature
    double curve_location_const = 4.0;

    ChVector<double> curvature_location =
        m_vehicle.GetVehiclePos() +
        curve_location_const * (m_acc_driver->GetSteeringController().GetTargetLocation() - m_vehicle.GetVehiclePos());

    double t;
    double t_actual;
    double min_dist = 1e6;
    int knot_id = 0;
    for (int i = 0; i < m_path->getNumPoints(); i++) {
        ChVector<> loc = m_path->calcClosestPoint(curvature_location, i, t);
        double tmp_min_dist = (loc - curvature_location).Length();
        if (tmp_min_dist < min_dist) {
            min_dist = tmp_min_dist;
            knot_id = i;
            t_actual = t;
        }
    }

    ChVector<> d = m_path->evalD(knot_id, t_actual);
    d.Normalize();
    ChVector<> heading = m_vehicle.GetVehicleRot().Rotate({1, 0, 0});
    double dotangle = d.Dot(heading);

    m_acc_driver->SetDesiredSpeed(m_target_speed * std::max(0.3, dotangle * 1.5 - .5));
    // std::cout << "Speed: " << m_target_speed * std::max(0.3, dotangle * 2.0 - 1.0) << std::endl;
    m_acc_driver->SetCurrentDistance(m_current_distance);
    m_acc_driver->Advance(step);

    double max_dt = 0.01;
    m_throttle = ChClamp(m_acc_driver->GetThrottle(), m_throttle - max_dt, m_throttle + max_dt);
    m_steering = ChClamp(m_acc_driver->GetSteering(), m_steering - max_dt, m_steering + max_dt);
    m_braking = ChClamp(m_acc_driver->GetBraking(), m_braking - max_dt, m_braking + max_dt);
}

}  // namespace synchrono
}  // end namespace chrono

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
// Interactive driver for a vehicle. This class implements the
// functionality required by its base ChDriver class using keyboard or joystick
// inputs. If a joystick is present it will use that as an input; it will
// otherwise default to a keyboard input.
//
// =============================================================================

#ifndef CH_SOUND_ENGINE_H
#define CH_SOUND_ENGINE_H
#include <string>

#include "chrono/physics/ChSystem.h"
#include "chrono/utils/ChUtilsChaseCamera.h"
#include "chrono_vehicle/ChApiVehicle.h"
#include "chrono_vehicle/ChVehicle.h"
#include "chrono_vehicle/wheeled_vehicle/ChWheeledVehicleVisualSystemIrrlicht.h"

#ifdef CHRONO_IRRKLANG
#include "irrKlang.h"

using namespace chrono::vehicle;

namespace chrono {

// Sound effect tools for the CSL simulator
class ChCSLSoundEngine {
  public:
    /// Construct the sound reproduction engine
    ChCSLSoundEngine(ChVehicle* vehicle);

    ~ChCSLSoundEngine();

    /// Updates sound engine
    void Synchronize(double time);

  private:
    ChVehicle* thisvehicle;
    irrklang::ISound* car_sound;
    // irrklang::ISound* motor_sound;
    irrklang::ISoundEngine* sound_engine;
    // std::vector<std::string> motor_soundfiles;
    // std::vector<irrklang::ISound*> motor_sounds;
    double last_time_played = 0;
    int last_threshold = 0;
};

}  // end namespace chrono

#endif
#endif

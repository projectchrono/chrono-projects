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
#include "ChCSLSoundEngine.h"
#include <string>

#include "chrono_vehicle/ChApiVehicle.h"
#ifdef CHRONO_IRRKLANG
namespace chrono {


ChCSLSoundEngine::ChCSLSoundEngine(ChVehicle* vehicle){
  thisvehicle = vehicle;
  sound_engine = irrklang::createIrrKlangDevice();
  /*motor_soundfiles = {"/600rpm_noload", "/1000_rpm_manifoldNorm", "/1500_rpm_manifoldNorm", 
                      "/2000_rpm_manifoldNorm", "/2500_rpm_manifoldNorm", "/3000_rpm_manifoldNorm", "/3500_rpm_manifoldNorm"};*/
  /*for(auto file : motor_soundfiles){
    irrklang::ISound* motor_sound = sound_engine->play2D((std::string(STRINGIFY(HIGHWAY_DATA_DIR)) + 
                                                          "/Environments/Iowa/Sounds" + file + ".wav").c_str(), true, false, true);
            motor_sound->setIsPaused(true);
            motor_sounds.push_back(motor_sound);
        }*/
    car_sound = sound_engine->play2D((std::string(STRINGIFY(HIGHWAY_DATA_DIR)) + 
                                                          "/Environments/Iowa/Sounds/audi.ogg").c_str(), true, false, true);
    car_sound->setIsPaused(true);

}

ChCSLSoundEngine::~ChCSLSoundEngine(){
  delete sound_engine;
}

void ChCSLSoundEngine::Synchronize(double time){
  //update every 0.01 sec
  if(time-last_time_played>0.5){
      last_time_played = time;
      double rpm = thisvehicle->GetPowertrainAssembly()->GetEngine()->GetMotorSpeed() * 60 / CH_C_2PI;
      double soundspeed = rpm / (2000.);  // denominator: to guess
      if (soundspeed < 0.1)
          soundspeed = 0.1;
      if (car_sound->getIsPaused())
          car_sound->setIsPaused(false);
      car_sound->setPlaybackSpeed((irrklang::ik_f32)soundspeed);
  }
    
    // we round it to change only every 500 rpm
    /*if(rpm > 500){
	    uint32_t new_threshold = (rpm - 500)/500;
	    if(new_threshold != last_threshold){
	      motor_sounds[last_threshold]->setIsPaused(true);
	      motor_sounds[new_threshold]->setIsPaused(false);
	      motor_sounds[new_threshold]->setVolume(1);
	      std::cout << "Playing new audio" << "\n";
	      last_threshold = new_threshold;
	      }
	    if(new_threshold > 0){
	      double soundspeed = rpm/(rpm-500);
	      std::cout << "Sound Playback Speed:  " << soundspeed << "\n";
	      motor_sounds[last_threshold]->setPlaybackSpeed((irrklang::ik_f32)soundspeed);
	      }*/
	}
}


#endif

// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2019 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Conlain Kelly, Ruochun Zhang
// =============================================================================
// Chrono::Granular simulation of up to one million spherical particles
// settling in a box in order to measure run time.
// =============================================================================

#include <fstream>
#include <iostream>
#include <string>

#include "chrono/core/ChTimer.h"
#include "chrono/utils/ChUtilsSamplers.h"
#include "chrono_gpu/physics/ChSystemGpu.h"
#include "chrono_thirdparty/filesystem/path.h"

#include "../utils.h"

using namespace chrono;
using namespace chrono::gpu;

std::string output_prefix = "../results";

// Default values
float ballRadius = 1.f;
float ballDensity = 2.50f;
float timeEnd = .25f;
float grav_acceleration = -980.f;
float normStiffness_S2S = 1e7f;
float normStiffness_S2W = 1e7f;

float normalDampS2S = 1000;
float normalDampS2W = 1000;
float adhesion_ratio_s2w = 0.0f;
float timestep = 1e-4f;

CHGPU_OUTPUT_MODE write_mode = CHGPU_OUTPUT_MODE::BINARY;
CHGPU_VERBOSITY verbose = CHGPU_VERBOSITY::INFO;
float cohesion_ratio = 0;

// -----------------------------------------------------------------------------
// Run a wavetank for a monodisperse collection of spheres in a rectangular box, undergoing a wave motion
// -----------------------------------------------------------------------------
double run_test(float box_size_X, float box_size_Y, float box_size_Z) {
    // Setup simulation
    ChSystemGpu gpu_system(ballRadius, ballDensity, ChVector3f(box_size_X, box_size_Y, box_size_Z));
    gpu_system.SetKn_SPH2SPH(normStiffness_S2S);
    gpu_system.SetKn_SPH2WALL(normStiffness_S2W);
    gpu_system.SetGn_SPH2SPH(normalDampS2S);
    gpu_system.SetGn_SPH2WALL(normalDampS2W);

    gpu_system.SetCohesionRatio(cohesion_ratio);
    gpu_system.SetAdhesionRatio_SPH2WALL(adhesion_ratio_s2w);
    gpu_system.SetGravitationalAcceleration(ChVector3d(0.f, 0.f, grav_acceleration));
    gpu_system.SetParticleOutputMode(write_mode);

    // Fill the bottom half with material
    chrono::utils::ChHCPSampler<float> sampler(2.4f * ballRadius);  // Add epsilon
    ChVector3f center(0, 0, -0.25f * box_size_Z);
    ChVector3f hdims(box_size_X / 2 - ballRadius, box_size_X / 2 - ballRadius, box_size_Z / 4 - ballRadius);
    std::vector<ChVector3f> body_points = sampler.SampleBox(center, hdims);
    gpu_system.SetParticles(body_points);

    filesystem::create_directory(filesystem::path(output_prefix));

    gpu_system.SetBDFixed(true);
    gpu_system.SetFrictionMode(CHGPU_FRICTION_MODE::FRICTIONLESS);
    gpu_system.SetTimeIntegrator(CHGPU_TIME_INTEGRATOR::EXTENDED_TAYLOR);
    gpu_system.SetVerbosity(verbose);

    gpu_system.SetFixedStepSize(timestep);
    ChTimer timer;

    // Run wavetank experiment and time it
    timer.start();
    gpu_system.Initialize();
    int fps = 50;
    // assume we run for at least one frame
    float frame_step = 1.0f / fps;
    float curr_time = 0;
    int currframe = 0;

    // Run settling experiments
    while (curr_time < timeEnd) {
        gpu_system.AdvanceSimulation(frame_step);
        curr_time += frame_step;
        printf("rendering frame %u\n", currframe);
        char filename[100];
        sprintf(filename, "%s/step%06d", output_prefix.c_str(), currframe++);
        gpu_system.WriteParticleFile(filename);
    }
    timer.stop();
    return timer.GetTimeSeconds();
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cout << "Usage:\n./test_GPU_milsettle <results_log_file>" << std::endl;
    }

    // up to one million bodies
    double time50k = run_test(100, 100, 100);
    double time500k = run_test(220, 220, 220);
    double time1mil = run_test(280, 280, 280);

    std::cout << "Running settling test!" << std::endl;
    std::cout << "50 thousand bodies took " << time50k << " seconds!" << std::endl;
    std::cout << "500 thousand bodies took " << time500k << " seconds!" << std::endl;
    std::cout << "1 million bodies took " << time1mil << " seconds!" << std::endl;

    std::ofstream ofile(argv[1], std::ofstream::app);
    ofile << "Running settling test!" << std::endl;
    ofile << "50 thousand bodies took " << time50k << " seconds!" << std::endl;
    ofile << "500 thousand bodies took " << time500k << " seconds!" << std::endl;
    ofile << "1 million bodies took " << time1mil << " seconds!" << std::endl;
    return 0;
}

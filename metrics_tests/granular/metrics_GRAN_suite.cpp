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
// Authors: Conlain Kelly, Nic Olsen
// =============================================================================
// Chrono::Granular metrics test suite of various parameters
// =============================================================================

#include <fstream>
#include <iostream>
#include <string>

#include "chrono/core/ChTimer.h"
#include "chrono/utils/ChUtilsSamplers.h"

#include "chrono_granular/physics/ChGranular.h"

#include "chrono_thirdparty/SimpleOpt/SimpleOpt.h"
#include "chrono_thirdparty/filesystem/path.h"

#include "../BaseTest.h"

using namespace chrono;
using namespace chrono::granular;

std::string output_prefix = "./results";
std::string delim = "--------------------------------";

// Default values
float ballRadius = 1.f;
float ballDensity = 2.50f;
float timeEnd = .25f;
float grav_acceleration = -980.f;
float normalStiffness_S2S = 5e7f;
float normalStiffness_S2W = 5e7f;
float normalDampS2S = 10000;
float normalDampS2W = 10000;

float tangentStiffness_S2S = 2e7f;
float tangentStiffness_S2W = 2e7f;
float tangentDampS2S = 0;
float tangentDampS2W = 0;

float timestep = 1e-4;

unsigned int psi_T = 32;
unsigned int psi_L = 16;

GRAN_OUTPUT_MODE write_mode = GRAN_OUTPUT_MODE::BINARY;

// about a half-million bodies
constexpr float3 full_box_size = {280, 280, 280};

int fps = 50;
// assume we run for at least one frame
float frame_step = 1.0f / fps;

// Test class
class GranSuiteTest : public BaseTest {
  public:
    GranSuiteTest(const std::string& testName, const std::string& testProjectName)
        : BaseTest(testName, testProjectName), m_execTime(0) {}

    ~GranSuiteTest() {}

    // Override corresponding functions in BaseTest
    virtual bool execute() override;
    virtual double getExecutionTime() const override { return m_execTime; }

  private:
    double m_execTime;
};

void setupBasicSystem(ChSystemGranularSMC& gran_sys, float3 box_size) {
    gran_sys.set_K_n_SPH2SPH(normalStiffness_S2S);
    gran_sys.set_K_n_SPH2WALL(normalStiffness_S2W);

    gran_sys.set_Gamma_n_SPH2SPH(normalDampS2S);
    gran_sys.set_Gamma_n_SPH2WALL(normalDampS2W);

    gran_sys.setPsiFactors(psi_T, psi_L);

    gran_sys.set_Cohesion_ratio(0);
    gran_sys.set_Adhesion_ratio_S2W(0);
    gran_sys.set_gravitational_acceleration(0.f, 0.f, grav_acceleration);
    gran_sys.setOutputDirectory(output_prefix);
    gran_sys.setOutputMode(write_mode);

    // Fill the bottom half with material
    chrono::utils::HCPSampler<float> sampler(2.4 * ballRadius);  // Add epsilon
    ChVector<float> center(0, 0, -.25 * box_size.z);
    ChVector<float> hdims(box_size.x / 2 - ballRadius, box_size.y / 2 - ballRadius, box_size.z / 4 - ballRadius);
    std::vector<ChVector<float>> body_points = sampler.SampleBox(center, hdims);
    gran_sys.setParticlePositions(body_points);

    gran_sys.set_BD_Fixed(true);
    gran_sys.setVerbose(GRAN_VERBOSITY::QUIET);

    gran_sys.set_fixed_stepSize(timestep);
}

double runGranularSystem(ChSystemGranularSMC& gran_sys, std::string fprefix) {
    float curr_time = 0;
    int currframe = 0;
    ChTimer<double> timer;
    timer.start();
    gran_sys.initialize();

    // Run settling experiments
    while (curr_time < timeEnd) {
        gran_sys.advance_simulation(frame_step);
        curr_time += frame_step;
        printf("rendering frame %u\n", currframe);
        char filename[100];
        sprintf(filename, "%s/%s%06d", output_prefix.c_str(), fprefix.c_str(), currframe++);
        gran_sys.writeFile(filename);
    }
    timer.stop();

    return timer.GetTimeSeconds();
}

// run the basic frictionless system
double runWarmupTest() {
    std::cout << delim << "\n"
              << "Running warmup!\n";
    // run a small system to get the GPU warmed up
    constexpr float3 warmup_box_size = {100, 100, 100};
    ChSystemGranularSMC gran_sys(ballRadius, ballDensity, warmup_box_size);

    setupBasicSystem(gran_sys, warmup_box_size);
    gran_sys.set_friction_mode(FRICTIONLESS);
    gran_sys.set_timeIntegrator(FORWARD_EULER);

    // Run settling experiment and time it
    return runGranularSystem(gran_sys, "nofric");
}

// run the basic frictionless system
double runFrictionlessTest() {
    std::cout << delim << "\n"
              << "Running frictionless euler test!\n";
    ChSystemGranularSMC gran_sys(ballRadius, ballDensity,
                                 make_float3(full_box_size.x, full_box_size.y, full_box_size.z));

    setupBasicSystem(gran_sys, full_box_size);
    gran_sys.set_friction_mode(FRICTIONLESS);
    gran_sys.set_timeIntegrator(FORWARD_EULER);

    // Run settling experiment and time it
    return runGranularSystem(gran_sys, "nofric_euler");
}

// run the basic frictionless system
double runFrictionlessChung() {
    ChSystemGranularSMC gran_sys(ballRadius, ballDensity,
                                 make_float3(full_box_size.x, full_box_size.y, full_box_size.z));
    std::cout << delim << "\n"
              << "Running chung test!\n";
    setupBasicSystem(gran_sys, full_box_size);
    gran_sys.set_friction_mode(FRICTIONLESS);
    gran_sys.set_timeIntegrator(CHUNG);

    // Run settling experiment and time it
    return runGranularSystem(gran_sys, "nofric_chung");
}
// run the basic frictionless system
double runFrictionlessCenteredDiff() {
    ChSystemGranularSMC gran_sys(ballRadius, ballDensity,
                                 make_float3(full_box_size.x, full_box_size.y, full_box_size.z));
    std::cout << delim << "\n"
              << "Running Centered Diff test!\n";
    setupBasicSystem(gran_sys, full_box_size);
    gran_sys.set_friction_mode(FRICTIONLESS);
    gran_sys.set_timeIntegrator(CENTERED_DIFFERENCE);

    // Run settling experiment and time it
    return runGranularSystem(gran_sys, "nofric_CD");
}

// run the basic frictionless system
double runMultistepTest() {
    ChSystemGranularSMC gran_sys(ballRadius, ballDensity,
                                 make_float3(full_box_size.x, full_box_size.y, full_box_size.z));
    std::cout << delim << "\n"
              << "Running multistep friction test!\n";
    setupBasicSystem(gran_sys, full_box_size);
    gran_sys.set_friction_mode(MULTI_STEP);
    gran_sys.set_timeIntegrator(FORWARD_EULER);

    gran_sys.set_K_t_SPH2SPH(tangentStiffness_S2S);
    gran_sys.set_K_t_SPH2WALL(tangentStiffness_S2W);
    gran_sys.set_Gamma_t_SPH2SPH(tangentDampS2S);
    gran_sys.set_Gamma_t_SPH2WALL(tangentDampS2W);

    // Run settling experiment and time it
    return runGranularSystem(gran_sys, "multistep_euler");
}

// run the basic frictionless system
double runSingleStepTest() {
    ChSystemGranularSMC gran_sys(ballRadius, ballDensity,
                                 make_float3(full_box_size.x, full_box_size.y, full_box_size.z));
    std::cout << delim << "\n"
              << "Running single step friction test!\n";
    setupBasicSystem(gran_sys, full_box_size);
    gran_sys.set_friction_mode(SINGLE_STEP);
    gran_sys.set_timeIntegrator(FORWARD_EULER);

    gran_sys.set_K_t_SPH2SPH(tangentStiffness_S2S);
    gran_sys.set_K_t_SPH2WALL(tangentStiffness_S2W);
    gran_sys.set_Gamma_t_SPH2SPH(tangentDampS2S);
    gran_sys.set_Gamma_t_SPH2WALL(tangentDampS2W);

    // Run settling experiment and time it
    return runGranularSystem(gran_sys, "singlestep_euler");
}

// run all tests
bool GranSuiteTest::execute() {
    // run a small frictionless test to make sure everything is warmed up
    runWarmupTest();
    // these are all the metrics tests to run
    double time_nofric = 0;
    double time_singlestep = 0;
    double time_nofricChung = 0;
    double time_nofricCenteredDiff = 0;
    double time_multistep = 0;

    time_nofric = runFrictionlessTest();
    time_singlestep = runSingleStepTest();
    time_nofricChung = runFrictionlessChung();
    time_nofricCenteredDiff = runFrictionlessCenteredDiff();
    time_multistep = runMultistepTest();

    addMetric("time_frictionless_euler", time_nofric);
    addMetric("time_frictionless_chung", time_singlestep);
    addMetric("time_single_step_euler", time_nofricChung);
    addMetric("time_frictionless_CD", time_nofricCenteredDiff);
    addMetric("time_multistep_euler", time_multistep);

    std::cout << "Running metrics suite!" << std::endl;
    std::cout << "Frictionless took " << time_nofric << " seconds!" << std::endl;
    std::cout << "Chung took " << time_nofricChung << " seconds!" << std::endl;
    std::cout << "CenteredDiff took " << time_nofricCenteredDiff << " seconds!" << std::endl;
    std::cout << "Single Step took " << time_singlestep << " seconds!" << std::endl;
    std::cout << "Multistep took " << time_multistep << " seconds!" << std::endl;

    return true;
}

int main(int argc, char* argv[]) {
    filesystem::create_directory(filesystem::path(output_prefix));

    std::string out_dir = "../METRICS";
    if (!filesystem::create_directory(filesystem::path(out_dir))) {
        std::cout << "Error creating directory " << out_dir << std::endl;
        return 1;
    }

    GranSuiteTest test("metrics_GRAN_suite", "Chrono::Granular");
    test.setOutDir(out_dir);
    test.setVerbose(true);
    bool passed = test.run();
    test.print();

    return 0;
}

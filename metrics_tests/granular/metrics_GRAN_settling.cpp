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
// Authors: Nic Olsen
// =============================================================================
// This test compares timing results of three systems:
// (1) Granular-only (2) Granular and a single disabled triangle mesh
// (3) Granular and a single enabled triangle mesh fixed in place
// =============================================================================

#include <cmath>
#include <ctime>
#include <iostream>
#include <string>
#include "chrono/utils/ChUtilsSamplers.h"
#include "chrono_granular/physics/ChGranular.h"
#include "chrono_granular/physics/ChGranularTriMesh.h"
#include "chrono_thirdparty/filesystem/path.h"

#include "chrono/physics/ChGlobal.h"

#include "chrono_granular/utils/ChGranularJsonParser.h"

#include "../BaseTest.h"

using namespace chrono;
using namespace chrono::granular;

using std::cout;
using std::endl;
using std::string;
using std::vector;

float sphere_radius = 1;
float sphere_density = 2.5;

float box_X = 220;
float box_Y = 220;
float box_Z = 220;

float step_size = 5e-5;
float time_end = 1;

float grav_X = 0;
float grav_Y = 0;
float grav_Z = -980;

float normalStiffS2S = 5e7;
float normalStiffS2W = 5e7;
float normalStiffS2M = 5e7;
float normalDampS2S = 10000;
float normalDampS2W = 10000;
float normalDampS2M = 10000;

float tangentStiffS2S = 2e7;
float tangentStiffS2M = 2e7;
float tangentStiffS2W = 2e7;
float tangentDampS2S = 100;
float tangentDampS2M = 100;
float tangentDampS2W = 100;

float cohesion_ratio = 0;
float adhesion_ratio_s2w = 0;
float adhesion_ratio_s2m = 0;
bool verbose = false;

unsigned int psi_T = 32;
unsigned int psi_L = 16;

std::string output_dir = "settling";
GRAN_OUTPUT_MODE write_mode = GRAN_OUTPUT_MODE::BINARY;

enum RUN_MODE { GRAN = 0, GRAN_TRI_DISABLED = 1, GRAN_TRI_ENABLED = 2 };

double fill_top;
double block_mass = 1;

// Test class
class GranSettlingTest : public BaseTest {
  public:
    GranSettlingTest(const std::string& testName, const std::string& testProjectName)
        : BaseTest(testName, testProjectName), m_execTime(0) {}

    ~GranSettlingTest() {}

    // Override corresponding functions in BaseTest
    virtual bool execute() override;
    virtual double getExecutionTime() const override { return m_execTime; }

  private:
    double m_execTime;
};

void ShowUsage() {
    cout << "usage: ./metrics_GRAN_settling <json_file> <optional: test index for single test>" << endl;
}

void SetupGranSystem(ChSystemGranularSMC& m_sys) {
    m_sys.set_K_n_SPH2SPH(normalStiffS2S);
    m_sys.set_K_n_SPH2WALL(normalStiffS2W);
    m_sys.set_K_t_SPH2SPH(tangentStiffS2S);
    m_sys.set_K_t_SPH2WALL(tangentStiffS2W);
    m_sys.set_Gamma_n_SPH2SPH(normalDampS2S);
    m_sys.set_Gamma_n_SPH2WALL(normalDampS2W);
    m_sys.set_Gamma_t_SPH2SPH(tangentDampS2S);
    m_sys.set_Gamma_t_SPH2WALL(tangentDampS2W);

    m_sys.set_Cohesion_ratio(cohesion_ratio);
    m_sys.set_Adhesion_ratio_S2W(adhesion_ratio_s2w);
    m_sys.set_gravitational_acceleration(grav_X, grav_Y, grav_Z);
    m_sys.set_friction_mode(chrono::granular::GRAN_FRICTION_MODE::MULTI_STEP);

    m_sys.setOutputMode(write_mode);
    m_sys.setOutputDirectory(output_dir);

    m_sys.set_timeIntegrator(GRAN_TIME_INTEGRATOR::CENTERED_DIFFERENCE);
    m_sys.set_fixed_stepSize(step_size);
    m_sys.set_BD_Fixed(true);

    // Fill domain with particles
    vector<ChVector<float>> body_points;
    double epsilon = 0.2 * sphere_radius;
    double spacing = 2 * sphere_radius + epsilon;

    utils::HCPSampler<float> sampler(spacing);
    double fill_bottom = -box_Z / 2 + 2 * spacing;
    fill_top = box_Z / 2 - sphere_radius - epsilon;
    ChVector<> hdims(box_X / 2 - sphere_radius - epsilon, box_Y / 2 - sphere_radius - epsilon, 0);
    for (double z = fill_bottom; z < fill_top; z += spacing) {
        ChVector<> center(0, 0, z);
        auto points = sampler.SampleBox(center, hdims);
        body_points.insert(body_points.end(), points.begin(), points.end());
    }

    cout << "Created " << body_points.size() << " spheres" << endl;

    m_sys.setParticlePositions(body_points);
}

void SetupGranTriSystem(ChSystemGranularSMC_trimesh& m_sys) {
    SetupGranSystem(m_sys);

    m_sys.set_K_n_SPH2MESH(normalStiffS2M);
    m_sys.set_K_t_SPH2MESH(tangentStiffS2M);
    m_sys.set_Gamma_n_SPH2MESH(normalDampS2M);
    m_sys.set_Gamma_t_SPH2MESH(tangentDampS2M);
    m_sys.set_Adhesion_ratio_S2M(adhesion_ratio_s2m);

    // Mesh values
    vector<string> mesh_filenames;
    string mesh_filename("../metrics_tests/granular/upward_plane_refined.obj");

    mesh_filenames.push_back(GetChronoDataFile(mesh_filename));

    vector<float3> mesh_scalings;
    float3 scaling = make_float3(box_X / 2, box_Y / 2, 1);
    mesh_scalings.push_back(scaling);

    vector<float> mesh_masses;
    mesh_masses.push_back(block_mass);

    std::vector<bool> mesh_inflated;
    std::vector<float> mesh_inflation_radii;
    mesh_inflated.push_back(false);
    mesh_inflation_radii.push_back(0);

    m_sys.load_meshes(mesh_filenames, mesh_scalings, mesh_masses, mesh_inflated, mesh_inflation_radii);
}

double RunTest(RUN_MODE run_mode) {
    double out_fps = 50;
    float frame_step = 1.0 / out_fps;

    clock_t start = std::clock();
    switch (run_mode) {
        case RUN_MODE::GRAN: {
            cout << "Running Granular system test..." << endl;
            ChSystemGranularSMC m_sys(sphere_radius, sphere_density, make_float3(box_X, box_Y, box_Z));
            SetupGranSystem(m_sys);
            filesystem::create_directory(filesystem::path(output_dir));
            m_sys.initialize();

            unsigned int currframe = 0;
            for (float t = 0; t < time_end; t += frame_step) {
                cout << "Rendering frame " << currframe << endl;
                char filename[100];
                sprintf(filename, "%s/step%06u", output_dir.c_str(), currframe++);
                m_sys.writeFile(string(filename));

                m_sys.advance_simulation(frame_step);
            }
            break;
        }
        case RUN_MODE::GRAN_TRI_DISABLED: {
            cout << "Running Granular system with disabled mesh test..." << endl;
            ChSystemGranularSMC_trimesh m_sys(sphere_radius, sphere_density, make_float3(box_X, box_Y, box_Z));
            SetupGranTriSystem(m_sys);
            m_sys.disableMeshCollision();
            filesystem::create_directory(filesystem::path(output_dir));

            unsigned int nSoupFamilies = m_sys.getNumTriangleFamilies();
            cout << nSoupFamilies << " soup families" << endl;
            double* meshSoupLocOri = new double[7 * nSoupFamilies];
            float* meshVel = new float[6 * nSoupFamilies]();

            m_sys.initialize();

            unsigned int currframe = 0;
            for (float t = 0; t < time_end; t += frame_step) {
                meshSoupLocOri[0] = 0;
                meshSoupLocOri[1] = 0;
                meshSoupLocOri[2] = -box_Z / 2 + 2 * sphere_radius;

                meshSoupLocOri[3] = 1;
                meshSoupLocOri[4] = 0;
                meshSoupLocOri[5] = 0;
                meshSoupLocOri[6] = 0;

                m_sys.meshSoup_applyRigidBodyMotion(meshSoupLocOri, meshVel);
                cout << "Rendering frame " << currframe << endl;
                char filename[100];
                sprintf(filename, "%s/step%06u", output_dir.c_str(), currframe++);
                m_sys.writeFile(string(filename));
                m_sys.write_meshes(string(filename));

                m_sys.advance_simulation(frame_step);
            }
            delete[] meshSoupLocOri;

            break;
        }
        case RUN_MODE::GRAN_TRI_ENABLED: {
            cout << "Running Granular system with enabled mesh test..." << endl;
            ChSystemGranularSMC_trimesh m_sys(sphere_radius, sphere_density, make_float3(box_X, box_Y, box_Z));
            SetupGranTriSystem(m_sys);
            m_sys.enableMeshCollision();
            filesystem::create_directory(filesystem::path(output_dir));

            unsigned int nSoupFamilies = m_sys.getNumTriangleFamilies();
            cout << nSoupFamilies << " soup families" << endl;
            double* meshSoupLocOri = new double[7 * nSoupFamilies];
            float* meshVel = new float[6 * nSoupFamilies]();

            m_sys.initialize();

            unsigned int currframe = 0;
            for (float t = 0; t < time_end; t += frame_step) {
                meshSoupLocOri[0] = 0;
                meshSoupLocOri[1] = 0;
                meshSoupLocOri[2] = -box_Z / 2 + 2 * sphere_radius;

                meshSoupLocOri[3] = 1;
                meshSoupLocOri[4] = 0;
                meshSoupLocOri[5] = 0;
                meshSoupLocOri[6] = 0;

                m_sys.meshSoup_applyRigidBodyMotion(meshSoupLocOri, meshVel);
                cout << "Rendering frame " << currframe << endl;
                char filename[100];
                sprintf(filename, "%s/step%06u", output_dir.c_str(), currframe++);
                m_sys.writeFile(string(filename));
                m_sys.write_meshes(string(filename));

                m_sys.advance_simulation(frame_step);
            }
            delete[] meshSoupLocOri;

            break;
        }
    }
    clock_t end = std::clock();
    double total_time = ((double)(end - start)) / CLOCKS_PER_SEC;

    return total_time;
}

bool GranSettlingTest::execute() {
    double time_gran = -1;  // RunTest(RUN_MODE::GRAN);
    double time_tri_disabled = RunTest(RUN_MODE::GRAN_TRI_DISABLED);
    double time_tri_enabled = RunTest(RUN_MODE::GRAN_TRI_ENABLED);

    cout << "================== Results ==================" << endl;
    cout << "Granular system: " << time_gran << " seconds" << endl;
    cout << "Granular system with disabled triangles: " << time_tri_disabled << " seconds" << endl;
    cout << "Granular system with enabled triangles: " << time_tri_enabled << " seconds" << endl;

    addMetric("time_no_mesh", time_gran);
    addMetric("time_mesh_disabled", time_tri_disabled);
    addMetric("time_mesh_enabled", time_tri_enabled);

    return true;
}

int main(int argc, char* argv[]) {
    // Set the path to the Chrono data folder
    SetChronoDataPath(CHRONO_DATA_DIR);

    std::string out_dir = "../METRICS";
    if (!filesystem::create_directory(filesystem::path(out_dir))) {
        std::cout << "Error creating directory " << out_dir << std::endl;
        return 1;
    }

    GranSettlingTest test("metrics_GRAN_settling", "Chrono::Granular");
    test.setOutDir(out_dir);
    test.setVerbose(true);
    bool passed = test.run();
    test.print();
    return 0;
}

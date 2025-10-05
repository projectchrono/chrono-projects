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
// Authors: Nic Olsen, Conlain Kelly, Ruochun Zhang
// =============================================================================
// Set of simple tests for validating low-level behavior of a Chrono::DEM
// system.
// =============================================================================

#include <cmath>
#include <iostream>
#include <string>

#include "chrono/core/ChDataPath.h"
#include "chrono/utils/ChUtilsSamplers.h"
#include "chrono_dem/physics/ChSystemDem.h"
#include "chrono_dem/utils/ChDemJsonParser.h"
#include "chrono_thirdparty/filesystem/path.h"

#include "../utils.h"

using namespace chrono;
using namespace chrono::dem;

std::string output_dir = "../test_results";
std::string delim = "--------------------------------";

// Default values
constexpr float sphere_radius = 1.f;
constexpr float sphere_density = 2.50f;
constexpr float grav_Z = -980.f;
constexpr float normalStiffness_S2S = 1e8;
constexpr float normalStiffness_S2W = 1e8;
constexpr float normalStiffness_S2M = 1e8;
constexpr float normalDampS2S = 10000;
constexpr float normalDampS2W = 10000;
constexpr float normalDampS2M = 10000;

constexpr float tangentStiffness_S2S = 3e7;
constexpr float tangentStiffness_S2W = 3e7;
constexpr float tangentStiffness_S2M = 3e7;
constexpr float tangentDampS2S = 500;
constexpr float tangentDampS2W = 500;
constexpr float tangentDampS2M = 500;

constexpr float static_friction_coeff = 0.5f;

constexpr float cohes = 0;

constexpr float timestep = 2e-5f;

constexpr unsigned int psi_T = 16;
constexpr unsigned int psi_L = 16;

float box_X = 400.f;
float box_Y = 100.f;
float box_Z = 50.f;
float timeEnd = 5.f;

double fill_top;
double step_mass = 1;
double step_height = -1;

constexpr int fps = 100;
constexpr float frame_step = 1.f / fps;

CHDEM_OUTPUT_MODE write_mode = CHDEM_OUTPUT_MODE::CSV;

// assume we run for at least one frame
float curr_time = 0;
int currframe = 0;

// Bowling ball starts on incline to accelerate
enum TEST_TYPE { ROTF = 0, PYRAMID = 1, ROTF_MESH = 2, PYRAMID_MESH = 3, MESH_STEP = 4, MESH_FORCE = 5 };

// Set common set of parameters for all tests
void setCommonParameters(ChSystemDem& dem_sys) {
    dem_sys.SetPsiFactors(psi_T, psi_L);
    dem_sys.SetKn_SPH2SPH(normalStiffness_S2S);
    dem_sys.SetKn_SPH2WALL(normalStiffness_S2W);
    dem_sys.SetGn_SPH2SPH(normalDampS2S);
    dem_sys.SetGn_SPH2WALL(normalDampS2W);

    dem_sys.SetKt_SPH2SPH(tangentStiffness_S2S);
    dem_sys.SetKt_SPH2WALL(tangentStiffness_S2W);
    dem_sys.SetGt_SPH2SPH(tangentDampS2S);
    dem_sys.SetGt_SPH2WALL(tangentDampS2W);

    dem_sys.SetCohesionRatio(cohes);
    dem_sys.SetAdhesionRatio_SPH2WALL(cohes);
    dem_sys.SetGravitationalAcceleration(ChVector3d(0, 0, grav_Z));
    dem_sys.SetParticleOutputMode(write_mode);
    dem_sys.SetStaticFrictionCoeff_SPH2SPH(static_friction_coeff);
    dem_sys.SetStaticFrictionCoeff_SPH2WALL(static_friction_coeff);

    dem_sys.SetRollingCoeff_SPH2SPH(static_friction_coeff / 2.);
    dem_sys.SetRollingCoeff_SPH2WALL(static_friction_coeff / 2.);

    dem_sys.SetTimeIntegrator(CHDEM_TIME_INTEGRATOR::CENTERED_DIFFERENCE);
    dem_sys.SetFixedStepSize(timestep);
    filesystem::create_directory(filesystem::path(output_dir));
    dem_sys.SetBDFixed(true);
}

void setCommonMeshParameters(ChSystemDemMesh& dem_sys) {
    dem_sys.SetKn_SPH2MESH(normalStiffness_S2M);
    dem_sys.SetGn_SPH2MESH(normalDampS2M);
    dem_sys.SetKt_SPH2MESH(tangentStiffness_S2M);
    dem_sys.SetGt_SPH2MESH(tangentDampS2M);
}

void writeGranFile(ChSystemDem& dem_sys) {
    printf("rendering frame %u\n", currframe);
    char filename[100];
    sprintf(filename, "%s/step%06d", output_dir.c_str(), ++currframe);
    dem_sys.WriteParticleFile(std::string(filename));
}

void advanceGranSim(ChSystemDem& dem_sys) {
    dem_sys.AdvanceSimulation(frame_step);
    curr_time += frame_step;
}

void run_ROTF() {
    ChSystemDemMesh dem_sys(sphere_radius, sphere_density, ChVector3f(box_X, box_Y, box_Z));
    setCommonParameters(dem_sys);

    float ramp_angle = (float)CH_PI / 4;
    // ramp normal is 45 degrees about y
    float nx = std::cos(ramp_angle);
    float nz = std::sin(ramp_angle);

    ChVector3d plane_normal(nx, 0.f, nz);
    printf("Plane normal: (%f, %f, %f)\n", plane_normal.x(), plane_normal.y(), plane_normal.z());
    // place so that plane intersects wall near z = 0
    ChVector3d plane_pos(-box_X / 2.f, 0.f, 0.);

    std::vector<ChVector3f> points;
    // start at far-x wall, halfway up
    ChVector3f sphere_pos(-box_X / 2.f + 2.f * sphere_radius, 0, 2 * sphere_radius);
    points.push_back(sphere_pos);
    dem_sys.SetParticles(points);

    printf("Plane pos: (%f, %f, %f)\n", plane_pos.x(), plane_pos.y(), plane_pos.z());

    size_t slope_plane_id = dem_sys.CreateBCPlane(plane_pos, plane_normal, true);
    // add bottom plane to capture bottom forces
    ChVector3d bot_plane_pos(0, 0, -box_Z / 2 + 2 * sphere_radius);
    ChVector3f bot_plane_normal(0, 0, 1);
    size_t bottom_plane_id = dem_sys.CreateBCPlane(bot_plane_pos, bot_plane_normal, true);
    // Finalize settings and Initialize for runtime

    dem_sys.SetFrictionMode(CHDEM_FRICTION_MODE::MULTI_STEP);
    dem_sys.SetRollingMode(CHDEM_ROLLING_MODE::NO_RESISTANCE);
    dem_sys.Initialize();

    ChVector3f reaction_forces;

    // total distance traveled parallel to slope
    float total_dist = (1 / std::cos(ramp_angle)) * box_Z / 2;
    float estimated_time_to_bot = std::sqrt(2 * total_dist / std::abs(grav_Z * std::cos(ramp_angle)));
    printf("total dist is %f, estimated time is %f\n", total_dist, estimated_time_to_bot);

    // Run settling experiments
    while (curr_time < timeEnd) {
        bool success = dem_sys.GetBCReactionForces(slope_plane_id, reaction_forces);
        if (!success) {
            printf("ERROR! Get contact forces for plane failed\n");
        } else {
            printf("curr time is %f, slope plane force is (%f, %f, %f) Newtons\n", curr_time, reaction_forces.x(),
                   reaction_forces.x(), reaction_forces.z());
        }

        success = dem_sys.GetBCReactionForces(bottom_plane_id, reaction_forces);
        if (!success) {
            printf("ERROR! Get contact forces for plane failed\n");
        } else {
            printf("curr time is %f, bottom plane force is (%f, %f, %f) Newtons\n", curr_time, reaction_forces.x(),
                   reaction_forces.y(), reaction_forces.z());
        }
        writeGranFile(dem_sys);
        advanceGranSim(dem_sys);
    }
}

void run_ROTF_MESH() {
    // overwrite olds system
    ChSystemDemMesh dem_sys(sphere_radius, sphere_density, ChVector3f(box_X, box_Y, box_Z));
    setCommonParameters(dem_sys);

    setCommonMeshParameters(dem_sys);

    // place so that plane intersects wall near z = 0
    ChVector3d plane_pos(-box_X / 2.f, 0.f, 0.);

    std::vector<ChVector3f> points;
    // start at far-x wall, halfway up
    ChVector3f sphere_pos(-box_X / 2.f + 2.f * sphere_radius, 0, 2 * sphere_radius);
    points.push_back(sphere_pos);
    dem_sys.SetParticles(points);

    printf("Plane pos: (%f, %f, %f)\n", plane_pos.x(), plane_pos.y(), plane_pos.z());

    // add bottom plane to capture bottom forces
    ChVector3d bot_plane_pos(0, 0, -box_Z / 2 + 2 * sphere_radius);
    ChVector3f bot_plane_normal(0, 0, 1);

    std::vector<std::string> mesh_filenames;
    std::vector<ChMatrix33<float>> mesh_rotscales;
    std::vector<ChVector3f> mesh_translations;
    std::vector<float> mesh_masses;

    ChMatrix33<float> mesh_scaling(ChVector3f(100, 100, 100));

    // make two plane meshes, one for ramp and one for bottom
    mesh_filenames.push_back(GetProjectsDataFile("dem/meshes/testsuite/square_plane_fine.obj"));
    mesh_rotscales.push_back(mesh_scaling);
    mesh_translations.push_back(ChVector3f(0, 0, 0));
    mesh_masses.push_back(10.f);

    mesh_filenames.push_back(GetProjectsDataFile("dem/meshes/testsuite/square_plane_fine.obj"));
    mesh_rotscales.push_back(mesh_scaling);
    mesh_translations.push_back(ChVector3f(0, 0, 0));
    mesh_masses.push_back(10.f);

    dem_sys.AddMeshes(mesh_filenames, mesh_translations, mesh_rotscales, mesh_masses);

    // Finalize settings and Initialize for runtime
    dem_sys.SetFrictionMode(CHDEM_FRICTION_MODE::MULTI_STEP);
    dem_sys.SetRollingMode(CHDEM_ROLLING_MODE::NO_RESISTANCE);
    dem_sys.Initialize();

    unsigned int nSoupFamilies = dem_sys.GetNumMeshes();
    printf("%u soup families\n", nSoupFamilies);

    double* meshSoupLocOri = new double[7 * nSoupFamilies];

    // bottom plane faces upwards
    auto bot_quat = QuatFromAngleY(0);

    // this is a quaternion
    auto rot_quat = QuatFromAngleY(CH_PI / 4);

    // Run settling experiments
    while (curr_time < timeEnd) {
        dem_sys.ApplyMeshMotion(0, bot_plane_pos, bot_quat, ChVector3d(0, 0, 0), ChVector3d(0, 0, 0));
        dem_sys.ApplyMeshMotion(1, plane_pos, rot_quat, ChVector3d(0, 0, 0), ChVector3d(0, 0, 0));
        writeGranFile(dem_sys);
        advanceGranSim(dem_sys);
    }
}

void run_PYRAMID() {
    ChSystemDemMesh dem_sys(sphere_radius, sphere_density, ChVector3f(box_X, box_Y, box_Z));

    setCommonParameters(dem_sys);

    timeEnd = 1;
    // slightly inflated diameter to ensure no penetration
    float diam_delta = 2.01f;
    // add plane just below origin
    ChVector3d bot_plane_pos(0, 0, -1.02 * sphere_radius);
    ChVector3d bot_plane_normal(0, 0, 1);
    size_t bottom_plane_id = dem_sys.CreateBCPlane(bot_plane_pos, bot_plane_normal, true);

    dem_sys.SetFrictionMode(CHDEM_FRICTION_MODE::MULTI_STEP);
    dem_sys.SetRollingMode(CHDEM_ROLLING_MODE::NO_RESISTANCE);
    ChVector3d reaction_forces;

    // just above origin
    ChVector3d base_sphere_1(0, 0, 0);
    // down the x a little
    ChVector3d base_sphere_2(diam_delta * sphere_radius, 0, 0);
    // top of the triangle
    ChVector3d base_sphere_3(diam_delta * sphere_radius * std::cos(CH_PI / 3),
                             diam_delta * sphere_radius * std::sin(CH_PI / 3), 0);
    // top of pyramid in middle (average x, y)
    ChVector3d top_sphere((base_sphere_1.x() + base_sphere_2.x() + base_sphere_3.x()) / 3.,
                          (base_sphere_1.y() + base_sphere_2.y() + base_sphere_3.y()) / 3.,
                          2.0 * sphere_radius * std::sin(CH_PI / 3));

    std::vector<ChVector3f> points;

    points.push_back(base_sphere_1);
    points.push_back(base_sphere_2);
    points.push_back(base_sphere_3);
    points.push_back(top_sphere);
    dem_sys.SetParticles(points);

    dem_sys.Initialize();

    while (curr_time < timeEnd) {
        writeGranFile(dem_sys);
        advanceGranSim(dem_sys);
    }
}

void run_PYRAMID_MESH() {
    ChSystemDemMesh dem_sys(sphere_radius, sphere_density, ChVector3f(box_X, box_Y, box_Z));
    setCommonParameters(dem_sys);

    setCommonMeshParameters(dem_sys);

    timeEnd = 1;
    // slightly inflated diameter to ensure no penetration
    float diam_delta = 2.01f;
    // add plane just below origin
    ChVector3d bot_plane_pos(0, 0, -1.02 * sphere_radius);
    ChVector3d bot_plane_normal(0, 0, 1);

    // Add mesh
    dem_sys.AddMesh(GetProjectsDataFile("dem/meshes/testsuite/tiny_triangle.obj"), ChVector3f(0),
                    ChMatrix33<float>(ChVector3f(1, 1, 1)), 10.0f);

    dem_sys.SetFrictionMode(CHDEM_FRICTION_MODE::MULTI_STEP);
    dem_sys.SetRollingMode(CHDEM_ROLLING_MODE::NO_RESISTANCE);
    ChVector3d reaction_forces;

    // just above origin
    ChVector3d base_sphere_1(0, 0, 0);
    // down the x a little
    ChVector3d base_sphere_2(diam_delta * sphere_radius, 0, 0);
    // top of the triangle
    ChVector3d base_sphere_3(diam_delta * sphere_radius * std::cos(CH_PI / 3),
                             diam_delta * sphere_radius * std::sin(CH_PI / 3), 0);
    // top of pyramid in middle (average x, y)
    ChVector3d top_sphere((base_sphere_1.x() + base_sphere_2.x() + base_sphere_3.x()) / 3.,
                          (base_sphere_1.y() + base_sphere_2.y() + base_sphere_3.y()) / 3.,
                          2.0 * sphere_radius * std::sin(CH_PI / 3));

    std::vector<ChVector3f> points;

    points.push_back(base_sphere_1);
    points.push_back(base_sphere_2);
    points.push_back(base_sphere_3);
    points.push_back(top_sphere);
    dem_sys.SetParticles(points);

    dem_sys.Initialize();

    unsigned int nSoupFamilies = dem_sys.GetNumMeshes();
    printf("%u soup families\n", nSoupFamilies);

    double* meshSoupLocOri = new double[7 * nSoupFamilies];

    // bottom plane faces upwards
    auto quat = QuatFromAngleY(0);

    while (curr_time < timeEnd) {
        dem_sys.ApplyMeshMotion(0, bot_plane_pos, quat, ChVector3d(0, 0, 0), ChVector3d(0, 0, 0));
        writeGranFile(dem_sys);
        char filename[100];

        sprintf(filename, "%s/step%06d_meshes", output_dir.c_str(), currframe);

        dem_sys.WriteMeshes(std::string(filename));

        advanceGranSim(dem_sys);
    }
}

void run_MESH_STEP() {
    ChSystemDemMesh dem_sys(sphere_radius, sphere_density, ChVector3f(box_X, box_Y, box_Z));
    setCommonParameters(dem_sys);
    setCommonMeshParameters(dem_sys);

    dem_sys.AddMesh(GetProjectsDataFile("dem/meshes/testsuite/step.obj"), ChVector3f(0),
                    ChMatrix33<float>(ChVector3f(box_X / 2, box_Y / 2, (float)step_height)), (float)step_mass);

    // Fill domain with particles
    std::vector<ChVector3f> body_points;
    double epsilon = 0.2 * sphere_radius;
    double spacing = 2 * sphere_radius + epsilon;

    utils::ChPDSampler<float> sampler((float)spacing);
    double fill_bottom = -box_Z / 2 + step_height + 2 * spacing;
    fill_top = box_Z / 2 - sphere_radius - epsilon;
    ChVector3d hdims(box_X / 2 - sphere_radius - epsilon, box_Y / 2 - sphere_radius - epsilon, 0);
    for (double z = fill_bottom; z < fill_top; z += spacing) {
        ChVector3d center(0, 0, z);
        auto points = sampler.SampleBox(center, hdims);
        body_points.insert(body_points.end(), points.begin(), points.end());
    }

    std::cout << "Created " << body_points.size() << " spheres" << std::endl;

    dem_sys.SetParticles(body_points);

    unsigned int nSoupFamilies = dem_sys.GetNumMeshes();
    std::cout << nSoupFamilies << " soup families" << std::endl;

    ChVector3d meshSoupLoc(0, 0, -box_Z / 2 + 2 * sphere_radius);
    auto quat = QuatFromAngleY(0);

    dem_sys.Initialize();

    for (float t = 0; t < timeEnd; t += frame_step) {
        dem_sys.ApplyMeshMotion(0, meshSoupLoc, quat, ChVector3d(0, 0, 0), ChVector3d(0, 0, 0));
        std::cout << "Rendering frame " << currframe << std::endl;
        char filename[100];
        sprintf(filename, "%s/step%06u", output_dir.c_str(), currframe);
        writeGranFile(dem_sys);
        dem_sys.WriteMeshes(std::string(filename));
        advanceGranSim(dem_sys);
    }
}

void run_MESH_FORCE() {
    // TODO Adapt sizing
    std::cout << "MESH_FORCE not implemented" << std::endl;
    return;

    ChSystemDemMesh dem_sys(sphere_radius, sphere_density, ChVector3f(box_X, box_Y, box_Z));
    setCommonParameters(dem_sys);
    setCommonMeshParameters(dem_sys);

    utils::ChHCPSampler<float> sampler(2.1f * sphere_radius);
    auto pos = sampler.SampleBox(ChVector3d(0, 0, 26), ChVector3d(38, 38, 10));

    auto n_spheres = pos.size();
    std::cout << "Created " << n_spheres << " spheres" << std::endl;
    double sphere_mass = sphere_density * 4.0 * CH_PI * sphere_radius * sphere_radius * sphere_radius / 3.0;

    double total_mass = sphere_mass * n_spheres;
    double sphere_weight = sphere_mass * std::abs(grav_Z);
    double total_weight = total_mass * std::abs(grav_Z);

    dem_sys.SetParticles(pos);

    // Add mesh
    dem_sys.AddMesh(GetProjectsDataFile("dem/meshes/testsuite/square_box.obj"), ChVector3f(0),
                    ChMatrix33<float>(ChVector3f(40, 40, 40)), 1.0f);

    unsigned int nSoupFamilies = dem_sys.GetNumMeshes();
    std::cout << nSoupFamilies << " soup families" << std::endl;

    // Triangle remains at the origin
    ChVector3d meshLoc(0, 0, 0);
    auto quat = QuatFromAngleY(0);

    dem_sys.Initialize();

    // Run a loop that is typical of co-simulation. For instance, the wheeled is moved a bit, which moves the
    // particles. Conversely, the particles impress a force and torque upon the mesh soup
    for (float t = 0; t < timeEnd; t += frame_step) {
        dem_sys.ApplyMeshMotion(0, meshLoc, quat, ChVector3d(0, 0, 0), ChVector3d(0, 0, 0));
        std::cout << "Rendering frame " << currframe << std::endl;
        char filename[100];
        sprintf(filename, "%s/step%06u", output_dir.c_str(), currframe);
        writeGranFile(dem_sys);
        dem_sys.WriteMeshes(std::string(filename));
        // dem_sys.checkSDCounts(std::string(filename) + "SU", true, false);
        ChVector3d force;
        ChVector3d torque;
        dem_sys.CollectMeshContactForces(0, force, torque);
        std::cout << "force_z: " << force.z() << "; total weight: " << total_weight << "; sphere weight "
                  << sphere_weight << std::endl;
        std::cout << "torque: " << torque.x() << ", " << torque.y() << ", " << torque.z() << std::endl;

        advanceGranSim(dem_sys);
    }
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cout << "Usage:\n./test_DEM_testsuite <test_type>" << std::endl;
        std::cout << "  test_type : 0 - ROTF, 1 - PYRAMID, 2 - ROTF_MESH, 3 - PYRAMID_MESH, 4 - MESH_STEP 5:MESH_FORCE"
                  << std::endl;
        return 1;
    }

    TEST_TYPE curr_test = static_cast<TEST_TYPE>(std::atoi(argv[1]));

    std::cout << "frame step is " << frame_step << std::endl;

    switch (curr_test) {
        case ROTF: {
            run_ROTF();
            break;
        }
        case ROTF_MESH: {
            run_ROTF_MESH();
            break;
        }
        case PYRAMID: {
            run_PYRAMID();
            break;
        }
        case PYRAMID_MESH: {
            run_PYRAMID_MESH();
            break;
        }
        case MESH_STEP: {
            run_MESH_STEP();
            break;
        }
        case MESH_FORCE: {
            run_MESH_FORCE();
            break;
        }
        default: {
            std::cout << "Invalid test" << std::endl;
            return 1;
        }
    }
    return 0;
}

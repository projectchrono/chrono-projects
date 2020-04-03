
// Chrono utility header files
#include "chrono/utils/ChUtilsGeometry.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsGenerators.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_parallel/solver/ChSystemDescriptorParallel.h"
#include "chrono_parallel/collision/ChNarrowphaseRUtils.h"

using namespace chrono;
using namespace chrono::collision;
using namespace chrono::utils;

// Container dimensions
const real conversion = .3048;  // meters per foot
								// note that when specifying dimensions for the geometry half lengths are used.
real dim_a = 20 * conversion * 0.5;    // length of end platforms default : 28.5
real dim_b = 5 * conversion * 0.5;     // length of slope at top
real dim_c = 12 * conversion * 0.5;    // length of submerged slope default: 48.5
real dim_d = 15 * conversion * 0.5;    // length of bottom default: 100
real dim_e = 8 * conversion * 0.5;     // full depth of trench
real dim_w = 14 * conversion * 0.5;    // width of trench default: 20
real dim_t = 5.0 / 12.0 * conversion;  // wall thickness default : 10
bool add_top = true;

// Initial vehicle position and orientation
ChVector<> initLoc( -(dim_d + (dim_b + dim_c) * 2 + dim_a  * .9) + 1.9, 0, dim_e * 2 + 1.0);
ChQuaternion<> initRot(1, 0, 0, 0);
double dist_end = (dim_d + (dim_b + dim_c) * 2 + dim_a * 1.6);  // When to apply brakes

float container_friction = 0.8f;

void InitializeObject(std::shared_ptr<ChBody> body,
                      double mass,
                      const ChVector<>& pos,
                      const ChQuaternion<>& rot,
                      bool collide,
                      bool fixed,
                      int collision_family,
                      int do_not_collide_with) {
    body->SetMass(mass);
    body->SetPos(pos);
    body->SetRot(rot);
    body->SetCollide(collide);
    body->SetBodyFixed(fixed);
    body->GetCollisionModel()->ClearModel();
    body->GetCollisionModel()->SetFamily(collision_family);
    body->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(do_not_collide_with);
}

void FinalizeObject(std::shared_ptr<ChBody> body, ChSystem* system) {
    assert(body->GetContactMethod() == system->GetContactMethod());

    body->GetCollisionModel()->BuildModel();
    system->AddBody(body);
}

void CreateContainer(ChSystemParallelNSC* system) {
	// the extra .03 makes the slope meet up with the end platform
	real dim_slope = sqrt(dim_e * dim_e + (dim_b + dim_c) * (dim_b + dim_c));  // + .03;
	real angle = atan(dim_e / (dim_b + dim_c));
	real width = dim_w + dim_t * 2.0;

    std::shared_ptr<ChBody> bottom_plate = std::make_shared<ChBody>(std::make_shared<ChCollisionModelParallel>());
    std::shared_ptr<ChBody> side_plate_1 = std::make_shared<ChBody>(std::make_shared<ChCollisionModelParallel>());
    std::shared_ptr<ChBody> side_plate_2 = std::make_shared<ChBody>(std::make_shared<ChCollisionModelParallel>());
    std::shared_ptr<ChBody> end_plate_1 = std::make_shared<ChBody>(std::make_shared<ChCollisionModelParallel>());
    std::shared_ptr<ChBody> end_plate_2 = std::make_shared<ChBody>(std::make_shared<ChCollisionModelParallel>());
    std::shared_ptr<ChBody> end_slope_1 = std::make_shared<ChBody>(std::make_shared<ChCollisionModelParallel>());
    std::shared_ptr<ChBody> end_slope_2 = std::make_shared<ChBody>(std::make_shared<ChCollisionModelParallel>());

    auto material = std::make_shared<ChMaterialSurfaceNSC>();
	material->SetFriction(container_friction);
	material->SetCompliance(1e-9f);
	material->SetCohesion(0.0f);

	Vector c_pos = Vector(0, 0, 0);
	InitializeObject(bottom_plate, 1, Vector(0, 0, 0) + c_pos, Quaternion(1, 0, 0, 0), true, true, 2, 2);
	InitializeObject(side_plate_1, 1, Vector(0, 0, 0) + c_pos, Quaternion(1, 0, 0, 0), true, true, 2, 2);
	InitializeObject(side_plate_2, 1, Vector(0, 0, 0) + c_pos, Quaternion(1, 0, 0, 0), true, true, 2, 2);
	InitializeObject(end_plate_1, 1, Vector(0, 0, 0) + c_pos, Quaternion(1, 0, 0, 0), true, true, 2, 2);
	InitializeObject(end_plate_2, 1, Vector(0, 0, 0) + c_pos, Quaternion(1, 0, 0, 0), true, true, 2, 2);
	InitializeObject(end_slope_1, 1, Vector(0, 0, 0) + c_pos, Quaternion(1, 0, 0, 0), true, true, 2, 2);
	InitializeObject(end_slope_2, 1, Vector(0, 0, 0) + c_pos, Quaternion(1, 0, 0, 0), true, true, 2, 2);

	// Bottom plate
	AddBoxGeometry(bottom_plate.get(), material, Vector(dim_d * 1.1, width, dim_t), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	// Bottom Plate Edge
	AddCylinderGeometry(bottom_plate.get(), material, dim_t * 1.05, width, Vector(-dim_d, 0, 0), Quaternion(1, 0, 0, 0));
	AddCylinderGeometry(bottom_plate.get(), material, dim_t * 1.05, width, Vector(dim_d, 0, 0), Quaternion(1, 0, 0, 0));

	// Side walls
	AddBoxGeometry(side_plate_1.get(), material, Vector(dim_d + (dim_b + dim_c + dim_a) * 2, dim_t, dim_e * 2.5 + dim_t * 2),
		Vector(0, +(dim_w + dim_t), dim_e * 2.2), Quaternion(1, 0, 0, 0));
	AddBoxGeometry(side_plate_2.get(), material, Vector(dim_d + (dim_b + dim_c + dim_a) * 2, dim_t, dim_e * 2.5 + dim_t * 2),
		Vector(0, -(dim_w + dim_t), dim_e * 2.2), Quaternion(1, 0, 0, 0));
	// End Platforms
	AddBoxGeometry(end_plate_1.get(), material, Vector(dim_a, width, dim_t),
		Vector(+(dim_d + dim_c * 2 + dim_b * 2 + dim_a), 0, dim_e * 2.0), Quaternion(1, 0, 0, 0));
	AddBoxGeometry(end_plate_2.get(), material, Vector(dim_a, width, dim_t),
		Vector(-(dim_d + dim_c * 2 + dim_b * 2 + dim_a), 0, dim_e * 2.0), Quaternion(1, 0, 0, 0));
	// Slopes
	AddBoxGeometry(end_slope_1.get(), material, Vector(dim_slope, dim_w + dim_t * 2.0, dim_t),
		Vector(+(dim_d + (dim_c + dim_b) + sin(angle) * dim_t * 0.5), 0, dim_e),
		Q_from_AngAxis(-angle, VECT_Y));
	AddBoxGeometry(end_slope_2.get(), material, Vector(dim_slope, dim_w + dim_t * 2.0, dim_t),
		Vector(-(dim_d + (dim_c + dim_b) + sin(angle) * dim_t * 0.5), 0, dim_e),
		Q_from_AngAxis(+angle, VECT_Y));

	// Slope-platform edge
	AddCylinderGeometry(end_plate_1.get(), material, dim_t * 1.02, width,
		Vector(+(dim_d + dim_c * 2 + dim_b * 2), 0, dim_e * 2.0), Quaternion(1, 0, 0, 0));
	AddCylinderGeometry(end_plate_2.get(), material, dim_t * 1.02, width,
		Vector(-(dim_d + dim_c * 2 + dim_b * 2), 0, dim_e * 2.0), Quaternion(1, 0, 0, 0));
	// Top
	if (add_top) {
		AddBoxGeometry(bottom_plate.get(), material, Vector(dim_d + (dim_b + dim_c + dim_a) * 2, width, dim_t),
			Vector(0, 0, dim_e * 4.7 + dim_t * 2.0), Quaternion(1, 0, 0, 0));
	}
	// End Caps
	AddBoxGeometry(end_plate_1.get(), material, Vector(dim_t * .25, width, dim_e * 1.2 + dim_t * 2.0),
		Vector(+(dim_d + dim_c * 2 + dim_b * 2 + dim_a * 2), 0, dim_e * 3.1 + dim_t),
		Quaternion(1, 0, 0, 0));
	AddBoxGeometry(end_plate_2.get(), material, Vector(dim_t * .25, width, dim_e * 1.2 + dim_t * 2.0),
		Vector(-(dim_d + dim_c * 2 + dim_b * 2 + dim_a * 2), 0, dim_e * 3.1 + dim_t),
		Quaternion(1, 0, 0, 0));

	FinalizeObject(bottom_plate, (ChSystemParallel*)system);
	FinalizeObject(side_plate_1, (ChSystemParallel*)system);
	FinalizeObject(side_plate_2, (ChSystemParallel*)system);
	FinalizeObject(end_plate_1, (ChSystemParallel*)system);
	FinalizeObject(end_plate_2, (ChSystemParallel*)system);
	FinalizeObject(end_slope_1, (ChSystemParallel*)system);
	FinalizeObject(end_slope_2, (ChSystemParallel*)system);

	Vector side_dim(dim_d + (dim_b + dim_c + dim_a) * 2, dim_t, dim_e * 2 + dim_t * 2);
	Vector top_dim(dim_d + (dim_b + dim_c + dim_a) * 2, width, dim_t);

    printf("Containment dims: [%f %f %f] [%f %f %f]\n", -side_dim.x(), -width, -side_dim.z(), side_dim.x(), width,
           side_dim.z());
}

#include "chrono_parallel/physics/Ch3DOFContainer.h"

real stiffness = 1.4e5;

real youngs_modulus = stiffness;
// real youngs_modulus_fluid = 1.4e6/2.0;
real poissons_ratio = 0.2;

real theta_c = 2.5e-2;
real theta_s = 7.5e-3;
real lame_lambda = youngs_modulus * poissons_ratio / ((1. + poissons_ratio) * (1. - 2. * poissons_ratio));
real lame_mu = youngs_modulus / (2. * (1. + poissons_ratio));
real alpha_flip = .95;
real hardening_coefficient = 0.0;
int mpm_iterations = 50;
real dist = 0;
real alpha = 0.001;
real contact_compliance = 1e-9;
real compliance = 1e-12;
real mass;
double fluid_r = 0.016;


double rho = 207;
double contact_mu = 0;

double chassis_mu = 0.00;

Ch3DOFContainer* dof_container;


static bool PointInTriangle(real2 p, real2 p0, real2 p1, real2 p2) {
	real s = p0.y * p2.x - p0.x * p2.y + (p2.y - p0.y) * p.x + (p0.x - p2.x) * p.y;
	real t = p0.x * p1.y - p0.y * p1.x + (p0.y - p1.y) * p.x + (p1.x - p0.x) * p.y;

	if ((s < 0) != (t < 0)) {
		return false;
	}
	real A = -p1.y * p2.x + p0.y * (p2.x - p1.x) + p0.x * (p1.y - p2.y) + p1.x * p2.y;
	if (A < 0.0) {
		s = -s;
		t = -t;
		A = -A;
	}
	return s > 0 && t > 0 && (s + t) < A;
}

real TriArea(real2 p0, real2 p1, real2 p2) {
	double dArea = ((p1.x - p0.x) * (p2.y - p0.y) - (p2.x - p0.x) * (p1.y - p0.y)) / 2.0;
	return (dArea > 0.0) ? dArea : -dArea;
}

void CreateFluid(ChSystemParallelNSC* system) {
	youngs_modulus = stiffness;

	dof_container->kernel_radius = fluid_r * 2;
	dof_container->collision_envelope = 0;
	dof_container->contact_recovery_speed = 50;
	dof_container->max_velocity = 20;
	dof_container->contact_compliance = 0;
	dof_container->contact_cohesion = 0;
	dof_container->contact_mu = contact_mu;

    if (ChFluidContainer* fluid_container = dynamic_cast<ChFluidContainer*>(dof_container)) {
        fluid_container->alpha = alpha;
        fluid_container->mass = 1;
        fluid_container->contact_compliance = contact_compliance;
        fluid_container->epsilon = 1e-8;

    } else if (ChParticleContainer* fluid_container = dynamic_cast<ChParticleContainer*>(dof_container)) {
        fluid_container->alpha = alpha;
        fluid_container->mass = 1;
    }

    if (ChFluidContainer* fluid_container = dynamic_cast<ChFluidContainer*>(dof_container)) {
        fluid_container->rho = rho;
        fluid_container->tau = 1e-3 * 2;

        fluid_container->viscosity = 1;  // Viscosity of water.
        fluid_container->enable_viscosity = false;
        fluid_container->artificial_pressure = false;
        fluid_container->artificial_pressure_k = .01;
        fluid_container->artificial_pressure_dq = .2 * fluid_container->kernel_radius;
        fluid_container->artificial_pressure_n = 4;

        dist = fluid_container->kernel_radius * .9;

    } else if (ChParticleContainer* fluid_container = dynamic_cast<ChParticleContainer*>(dof_container)) {
        fluid_container->mu = 0;
        fluid_container->cohesion = 0;
        fluid_container->kernel_radius *= .9;
        fluid_container->mu = 0;
        fluid_container->compliance = compliance;
        dist = fluid_container->kernel_radius;
    }

    real offset_z = dof_container->kernel_radius;
    real height = dim_e * .75;

	utils::HCPSampler<> sampler(dist);
	utils::Generator::PointVector points = sampler.SampleBox(
		ChVector<>(0, 0, dim_t + offset_z + height * .5),
		ChVector<>((dim_d + dim_c + dim_b) * 2, dim_w - dof_container->kernel_radius * 2, +height * .5 + offset_z));

	real slope = (dim_e) / (dim_c + dim_b);
	real x_pos = (height + offset_z) / slope;

	std::vector<real3> pos_fluid;
	std::vector<real3> vel_fluid;
	real2 A = real2(-(dim_d + x_pos), height + dim_t + offset_z + offset_z);
	real2 B = real2((dim_d + x_pos), height + dim_t + offset_z + offset_z);
	real2 C = real2((dim_d), dim_t + offset_z);
	real2 D = real2((-dim_d), dim_t + offset_z);

	for (int i = 0; i < points.size(); i++) {
		real2 p = real2(points[i].x(), points[i].z());
		if (PointInTriangle(p, A, B, D) || PointInTriangle(p, B, C, D)) {
			pos_fluid.push_back(real3(points[i].x(), points[i].y(), points[i].z()));
			vel_fluid.push_back(real3(0));
		}
	}
    real vol = dist * dist * dist * .8;
    mass = rho * vol;
    std::cout << "fluid_mass: " << mass << " " << pos_fluid.size() << std::endl;
    if (ChFluidContainer* fluid_container = dynamic_cast<ChFluidContainer*>(dof_container)) {
        fluid_container->mass = mass;
        fluid_container->AddBodies(pos_fluid, vel_fluid);
    } else if (ChParticleContainer* fluid_container = dynamic_cast<ChParticleContainer*>(dof_container)) {
        fluid_container->mass = mass;
        fluid_container->AddBodies(pos_fluid, vel_fluid);
    }
}

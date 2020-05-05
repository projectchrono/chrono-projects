#pragma once

#include <iostream>
#include <sstream>
#include <cstring>
#include <thread>
#include <zlib.h>

#include "chrono/assets/ChVisualization.h"
#include "chrono/assets/ChSphereShape.h"
#include "chrono/assets/ChEllipsoidShape.h"
#include "chrono/assets/ChBoxShape.h"
#include "chrono/assets/ChConeShape.h"
#include "chrono/assets/ChCylinderShape.h"

#include "chrono/utils/ChUtilsGeometry.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsInputOutput.h"
#include "chrono/utils/ChUtilsGenerators.h"

#include "chrono_parallel/physics/ChSystemParallel.h"

using namespace chrono;
using namespace chrono::collision;
using std::cout;
using std::endl;

class CSVGen {
public:
	CSVGen() {
		delim = ",";
		binary = false;
		gz_file = 0;
	}
	~CSVGen() {}

	void OpenFile(std::string filename, bool bin = false) {
		if (bin) {
			binary = true;
			gz_file = gzopen(filename.c_str(), "wb");
		}
		else {
			ofile.open(filename.c_str(), std::ios::out);
		}
	}

	void CloseFile() {
		if (binary) {
			unsigned long int file_size = sizeof(char) * ss.str().size();

			gzwrite(gz_file, (void*)&file_size, sizeof(file_size));  // writing size of file
			gzwrite(gz_file, (void*)(ss.str().data()), file_size);
			gzclose(gz_file);
		}
		else {
			ofile << ss.str();
			ofile.close();
		}
	}
	template <class T>
	void operator<<(const T& token) {
		WriteToken(token);
	}
	void WriteToken(real token) { ss << token << delim; }
	void WriteToken(real2 token) { ss << token.x << delim << token.y << delim; }
	void WriteToken(real3 token) { ss << token.x << delim << token.y << delim << token.z << delim; }
	void WriteToken(real4 token) { ss << token.x << delim << token.y << delim << token.z << delim << token.w << delim; }
	void WriteToken(quaternion token) {
		ss << token.w << delim << token.x << delim << token.y << delim << token.z << delim;
	}
	void WriteToken(ChVector<> token) { ss << token.x() << delim << token.y() << delim << token.z() << delim; }
	void WriteToken(ChQuaternion<> token) {
		ss << token.e0() << delim << token.e1() << delim << token.e2() << delim << token.e3() << delim;
	}
	void WriteToken(std::string token) { ss << token << delim; }
	void endline() { ss << std::endl; }

	std::string delim;
	std::ofstream ofile;
	std::stringstream ss;
	gzFile gz_file;
	bool binary;
};


class BinaryGen {
public:
	BinaryGen() {}
	void OpenFile(std::string filename) { bin_file.open(filename.c_str(), std::ios::out | std::ofstream::binary); }
	void CloseFile() { bin_file.close(); }
	template <class T>
	void operator<<(const std::vector<T>& data) {
		size_t size = data.size();
		// write size of data type and then the data;
		bin_file.write(reinterpret_cast<const char*>(&size), sizeof(size));
		bin_file.write(reinterpret_cast<const char*>(&data[0]), size * sizeof(T));
	}
	template <typename T>
	void Write(const std::vector<T>& data) {
		size_t size = data.size();
		bin_file.write(reinterpret_cast<const char*>(&size), sizeof(size));
		bin_file.write(reinterpret_cast<const char*>(&data[0]), size * sizeof(T));
	}
	std::ofstream bin_file;
};

static std::vector<std::thread> writethreads;

void static WriteLocalData(const std::string&& filename,
	const std::vector<real3>&& pos_particle,
	const std::vector<real3>&& vel_particle,
	const int&& num_particles) {
	BinaryGen bin_output;
	bin_output.OpenFile(filename.c_str());
	bin_output.Write(pos_particle);
	bin_output.Write(vel_particle);
	bin_output.CloseFile();
}

void static DumpFluidData(chrono::ChSystemParallelNSC* system, std::string filename, bool binary = true) {
	std::vector<real3> pos_particle;
	std::vector<real3> vel_particle;
	int num_particles = system->data_manager->num_fluid_bodies;
	if (num_particles <= 0) {
		std::cout << "No fluid to write!!\n";
		return;
	}

	pos_particle = system->data_manager->host_data.pos_3dof;
	vel_particle = system->data_manager->host_data.vel_3dof;

	writethreads.push_back(std::thread(WriteLocalData, std::move(filename), std::move(pos_particle),
		std::move(vel_particle), std::move(num_particles)));
	if (writethreads.size() > 8) {
		for (std::thread& t : writethreads) {
			t.join();
		}
		writethreads.clear();
	}
}


void static DumpAllObjectsWithGeometryPovray(ChSystem* mSys, std::string filename, bool binary = false) {
	CSVGen csv_output;
    csv_output.OpenFile(filename.c_str(), binary);

    for (int i = 0; i < mSys->Get_bodylist().size(); i++) {
        auto abody = mSys->Get_bodylist().at(i);  // Get body
        const Vector pos = abody->GetFrame_REF_to_abs().GetPos();
        const Vector vel = abody->GetPos_dt();  // Get velocity
        Quaternion rot = abody->GetFrame_REF_to_abs().GetRot();
        Vector pos_final, rad_final;
        ChCollisionShape::Type type = ChCollisionShape::Type::SPHERE;
        // Get each asset
        for (int j = 0; j < abody->GetAssets().size(); j++) {
            auto asset = abody->GetAssets().at(j);
            if (!std::dynamic_pointer_cast<ChVisualization>(asset)) {
                continue;
            }
            ChVisualization* visual_asset = ((ChVisualization*)(asset.get()));
            Vector center = visual_asset->Pos;
            center = rot.Rotate(center);
            pos_final = pos + center;
            Quaternion lrot = visual_asset->Rot.Get_A_quaternion();
            lrot = rot * lrot;
            lrot.Normalize();
			if (std::dynamic_pointer_cast<ChSphereShape>(asset)) {
				ChSphereShape* sphere_shape = ((ChSphereShape*)(asset.get()));
				real radius = sphere_shape->GetSphereGeometry().rad;
				rad_final.x() = radius;
				rad_final.y() = radius;
				rad_final.z() = radius;
                type = ChCollisionShape::Type::SPHERE;
			}

			else if (std::dynamic_pointer_cast<ChEllipsoidShape>(asset)) {
				ChEllipsoidShape* ellipsoid_shape = ((ChEllipsoidShape*)(asset.get()));
				rad_final = ellipsoid_shape->GetEllipsoidGeometry().rad;
                type = ChCollisionShape::Type::ELLIPSOID;
			}
			else if (std::dynamic_pointer_cast<ChBoxShape>(asset)) {
				ChBoxShape* box_shape = ((ChBoxShape*)(asset.get()));
				rad_final = box_shape->GetBoxGeometry().Size;
                type = ChCollisionShape::Type::BOX;
			}
			else if (std::dynamic_pointer_cast<ChCylinderShape>(asset)) {
				ChCylinderShape* cylinder_shape = ((ChCylinderShape*)(asset.get()));
				double rad = cylinder_shape->GetCylinderGeometry().rad;
				double height = cylinder_shape->GetCylinderGeometry().p1.y() - cylinder_shape->GetCylinderGeometry().p2.y();
				rad_final.x() = rad;
				rad_final.y() = height;
				rad_final.z() = rad;
                type = ChCollisionShape::Type::CYLINDER;
			}
			else if (std::dynamic_pointer_cast<ChConeShape>(asset)) {
				ChConeShape* cone_shape = ((ChConeShape*)(asset.get()));
				rad_final.x() = cone_shape->GetConeGeometry().rad.x();
				rad_final.y() = cone_shape->GetConeGeometry().rad.y();
				rad_final.z() = cone_shape->GetConeGeometry().rad.z();
                type = ChCollisionShape::Type::CONE;
			}

			csv_output << pos_final;
			csv_output << lrot;
			csv_output << vel;

			if (std::dynamic_pointer_cast<ChSphereShape>(asset)) {
				csv_output << type;
				csv_output << rad_final.x();
				csv_output.endline();
			}
			else if (std::dynamic_pointer_cast<ChEllipsoidShape>(asset)) {
				csv_output << type;
				csv_output << real3(rad_final.x(), rad_final.y(), rad_final.z());
				csv_output.endline();
			}
			else if (std::dynamic_pointer_cast<ChBoxShape>(asset)) {
				csv_output << type;
				csv_output << real3(rad_final.x(), rad_final.y(), rad_final.z());
				csv_output.endline();
			}
			else if (std::dynamic_pointer_cast<ChCylinderShape>(asset)) {
				csv_output << type;
				csv_output << real2(rad_final.x(), rad_final.y());
				csv_output.endline();
			}
			else if (std::dynamic_pointer_cast<ChConeShape>(asset)) {
				csv_output << type;
				csv_output << real2(rad_final.x(), rad_final.y());
				csv_output.endline();
			}
			else {
				csv_output << -1;
				csv_output.endline();
			}
		}
	}

	csv_output.CloseFile();
}



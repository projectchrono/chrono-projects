#pragma once

#include <cstdint>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <thread>

#include "chrono/assets/ChVisualShapeSphere.h"
#include "chrono/assets/ChVisualShapeEllipsoid.h"
#include "chrono/assets/ChVisualShapeBox.h"
#include "chrono/assets/ChVisualShapeCone.h"
#include "chrono/assets/ChVisualShapeCylinder.h"
#include "chrono/collision/ChCollisionShape.h"

#include "chrono/utils/ChUtilsGeometry.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsGenerators.h"

#include "chrono_multicore/physics/ChSystemMulticore.h"

// Use the zlib compressor from stb_image_write (bundled with Chrono) to write gzip files
#define STB_IMAGE_WRITE_STATIC
#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STBI_WRITE_NO_STDIO
#include "chrono_thirdparty/stb/stb_image_write.h"

using namespace chrono;
using std::cout;
using std::endl;

// CRC-32 checksum, as required in the gzip trailer.
static uint32_t Crc32(const std::string& data) {
    uint32_t crc = 0xffffffffu;
    for (unsigned char c : data) {
        crc ^= c;
        for (int k = 0; k < 8; k++)
            crc = (crc >> 1) ^ (0xedb88320u & (0u - (crc & 1u)));
    }
    return ~crc;
}

// Write the given data to a gzip file (RFC 1952).
// The deflate stream from stbi_zlib_compress is wrapped in a zlib header and Adler-32 trailer, which are replaced with
// the gzip header and trailer.
static bool WriteGzipFile(const std::string& filename, const std::string& data) {
    int zlen = 0;
    unsigned char* zdata = stbi_zlib_compress((unsigned char*)data.data(), (int)data.size(), &zlen, 8);
    if (!zdata)
        return false;

    std::ofstream ofile(filename, std::ios::out | std::ios::binary);
    const unsigned char header[10] = {0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 255};  // deflate, no flags, unknown OS
    ofile.write((const char*)header, sizeof(header));
    ofile.write((const char*)zdata + 2, zlen - 6);
    free(zdata);

    uint32_t crc = Crc32(data);
    uint32_t size = (uint32_t)data.size();
    unsigned char trailer[8];
    for (int i = 0; i < 4; i++) {
        trailer[i] = (unsigned char)(crc >> (8 * i));
        trailer[4 + i] = (unsigned char)(size >> (8 * i));
    }
    ofile.write((const char*)trailer, sizeof(trailer));

    return ofile.good();
}

class CSVGen {
public:
	CSVGen() {
		delim = ",";
        binary = false;
    }
	~CSVGen() {}

    // Open the output file. In binary mode, the output is a gzip file containing the size of the text (uint64_t),
    // followed by the text.
    void OpenFile(std::string filename, bool bin = false) {
		if (bin) {
			binary = true;
            gz_filename = filename;
        }
		else {
			ofile.open(filename.c_str(), std::ios::out);
		}
	}

	void CloseFile() {
		if (binary) {
            std::string text = ss.str();
            uint64_t file_size = text.size();

            std::string data((const char*)&file_size, sizeof(file_size));  // size of file
            data += text;
            if (!WriteGzipFile(gz_filename, data))
                std::cerr << "Error writing " << gz_filename << std::endl;
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
    void WriteToken(int token) { ss << token << delim; }
    void WriteToken(float token) { ss << token << delim; }
    void WriteToken(double token) { ss << token << delim; }
    void WriteToken(real2 token) { ss << token.x << delim << token.y << delim; }
	void WriteToken(real3 token) { ss << token.x << delim << token.y << delim << token.z << delim; }
	void WriteToken(real4 token) { ss << token.x << delim << token.y << delim << token.z << delim << token.w << delim; }
	void WriteToken(quaternion token) {
		ss << token.w << delim << token.x << delim << token.y << delim << token.z << delim;
	}
    void WriteToken(ChVector3d token) { ss << token.x() << delim << token.y() << delim << token.z() << delim; }
    void WriteToken(ChQuaterniond token) { ss << token.e0() << delim << token.e1() << delim << token.e2() << delim << token.e3() << delim; }
    void WriteToken(std::string token) { ss << token << delim; }
	void endline() { ss << std::endl; }

	std::string delim;
	std::ofstream ofile;
	std::stringstream ss;
    std::string gz_filename;
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

void static DumpFluidData(chrono::ChSystemMulticoreNSC* system, std::string filename, bool binary = true) {
	std::vector<real3> pos_particle;
	std::vector<real3> vel_particle;
    int num_particles = system->data_manager->num_particles;
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

// Wait for all pending fluid data output to be written.
void static FinishFluidDataOutput() {
    for (std::thread& t : writethreads)
        t.join();
    writethreads.clear();
}

// Write one line for each visual shape in the system: absolute position and orientation of the shape, body velocity,
// shape type (a ChCollisionShape::Type value, -1 for unsupported shapes), and shape dimensions: radius for spheres,
// semiaxes for ellipsoids, half-lengths for boxes, and radius and height for cylinders and cones. Cylinders and cones
// are reported with their axis along the Y axis of the shape frame.
void static DumpAllObjectsWithGeometryPovray(ChSystem* mSys, std::string filename, bool binary = false) {
    CSVGen csv_output;
    csv_output.OpenFile(filename.c_str(), binary);

    for (const auto& body : mSys->GetBodies()) {
        auto vis_model = body->GetVisualModel();
        if (!vis_model)
            continue;

        ChFramed body_frame = body->GetVisualModelFrame();
        ChVector3d vel = body->GetPosDt();

        for (const auto& instance : vis_model->GetShapeInstances()) {
            ChFramed frame = body_frame.TransformLocalToParent(instance.frame);
            ChQuaterniond rot = frame.GetRot();

            int type = -1;
            int num_dims = 0;
            ChVector3d dims;
            if (auto sphere = std::dynamic_pointer_cast<ChVisualShapeSphere>(instance.shape)) {
                type = ChCollisionShape::Type::SPHERE;
                num_dims = 1;
                dims = ChVector3d(sphere->GetRadius());
            } else if (auto ellipsoid = std::dynamic_pointer_cast<ChVisualShapeEllipsoid>(instance.shape)) {
                type = ChCollisionShape::Type::ELLIPSOID;
                num_dims = 3;
                dims = ellipsoid->GetSemiaxes();
            } else if (auto box = std::dynamic_pointer_cast<ChVisualShapeBox>(instance.shape)) {
                type = ChCollisionShape::Type::BOX;
                num_dims = 3;
                dims = box->GetHalflengths();
            } else if (auto cylinder = std::dynamic_pointer_cast<ChVisualShapeCylinder>(instance.shape)) {
                type = ChCollisionShape::Type::CYLINDER;
                num_dims = 2;
                dims = ChVector3d(cylinder->GetRadius(), cylinder->GetHeight(), 0);
                rot = rot * Q_ROTATE_Y_TO_Z;
            } else if (auto cone = std::dynamic_pointer_cast<ChVisualShapeCone>(instance.shape)) {
                type = ChCollisionShape::Type::CONE;
                num_dims = 2;
                dims = ChVector3d(cone->GetRadius(), cone->GetHeight(), 0);
                rot = rot * Q_ROTATE_Y_TO_Z;
            }

            csv_output << frame.GetPos();
            csv_output << rot;
            csv_output << vel;
            csv_output << type;
            for (int i = 0; i < num_dims; i++)
                csv_output << dims[i];
            csv_output.endline();
        }
    }

    csv_output.CloseFile();
}



#include "chrono/utils/ChUtilsGeometry.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsInputOutput.h"
#include "chrono/utils/ChUtilsGenerators.h"

void static LoadConvexHull(const std::string& file_name,
	geometry::ChTriangleMeshConnected& convex_mesh,
	std::vector<std::vector<ChVector<double> > >& convex_hulls,
	const ChVector<>& pos,
	const ChQuaternion<>& rot) {
	convex_mesh.LoadWavefrontMesh(file_name, true, false);

	std::vector<tinyobj::shape_t> shapes;
	std::string err = tinyobj::LoadObj(shapes, file_name.c_str());

	convex_hulls.resize(shapes.size());
	std::cout << "NUM HULLS: " << shapes.size() << std::endl;
	for (int i = 0; i < shapes.size(); i++) {
		convex_hulls[i].resize(shapes[i].mesh.input_pos.size() / 3);
		std::cout << "HULL: " << i << " "
			<< shapes[i].mesh.input_pos.size() / 3 << std::endl;

		for (int j = 0; j < shapes[i].mesh.input_pos.size() / 3; j++) {
			ChVector<double> pos(shapes[i].mesh.input_pos[j * 3 + 0], shapes[i].mesh.input_pos[j * 3 + 1],
				shapes[i].mesh.input_pos[j * 3 + 2]);

			//std::cout << pos.x << " " << pos.y << " " << pos.z << "\n";

			convex_hulls[i][j] = pos;
		}
	}
}
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdio>
#include <cstdlib>

#include <math.h>
#include <vector>
#include <functional>
#include <array>
#include <algorithm>

#include "../math/Stride.h"
#include "settings.h"
#include "material_info.hpp"
#include "solver2d.hpp"
#include "solver3d.hpp"
#include "../math/Spaces.h"
#include "../mesh/MeshIO.h"

using Space = Space3;
using MaterialInfoIndex = int;


template <typename Space>
typename Space::IndexType idx(const typename Space::IndexVector& coord, const typename Space::IndexVector& stride)
{
	SPACE_TYPEDEFS;

	IndexType resultIdx = 0;
	IndexType mult = 1;

	for (auto dimIndex = 0; dimIndex < Space::Dimension; ++dimIndex)
	{
		resultIdx += coord.Get(dimIndex) * mult;
		mult *= stride.Get(dimIndex);
	}
	return resultIdx;
}

std::vector<MaterialInfoIndex> loadMesh(
	const std::string& fileName,
	const Space::IndexVector& stride,
	Space::AABB& meshAABB)
{
	SPACE_TYPEDEFS;

	MeshIO<Space> meshIO;
	meshIO.Load(AddExtensionToFileName(fileName, ".mesh"), IO::Ascii);
	meshAABB = getAABB(meshIO);

	std::string paramsFileName = AddExtensionToFileName(fileName, ".params");
	FILE* paramsFile = nullptr;
	fopen_s(&paramsFile, paramsFileName.c_str(), "rb");
	assert(paramsFile);

	std::vector<MaterialInfoIndex> cellMaterialInfos(meshIO.GetCellsCount());

	if (paramsFile)
	{
		for (auto cellIndex = 0; cellIndex < meshIO.GetCellsCount(); ++cellIndex)
		{
			char currCellSubmeshIndex = char(-1);
			IndexType bytesRead = fread(&currCellSubmeshIndex, sizeof(char), 1, paramsFile);
			assert(bytesRead == 1);

			cellMaterialInfos[cellIndex] = currCellSubmeshIndex;
		}
		fclose(paramsFile);
	}

	std::vector<MaterialInfoIndex> materialInfoIndices(stride.GetVolume(), 0 /* todo: default value from config */);
	
	for (IndexType cellIndex = 0; cellIndex < meshIO.GetCellsCount(); ++cellIndex)
	{
		Vector points[Space::NodesPerCell];
		for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
		{
			auto nodeIndex = meshIO.indices[cellIndex * Space::NodesPerCell + nodeNumber];
			points[nodeNumber] = meshIO.vertices[nodeIndex];
		}

		std::array<Scalar, Space::Dimension> mins;
		std::array<Scalar, Space::Dimension> maxs;

		for (auto dimIndex = 0; dimIndex < Space::Dimension; ++dimIndex)
		{
			maxs[dimIndex] = mins[dimIndex] = points[0].Get(dimIndex);
		}

		for (auto nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
		{
			for (auto dimIndex = 0; dimIndex < Space::Dimension; ++dimIndex)
			{
				if (points[nodeNumber].Get(dimIndex) < mins[dimIndex]) mins[dimIndex] = points[nodeNumber].Get(dimIndex);
				if (points[nodeNumber].Get(dimIndex) > maxs[dimIndex]) maxs[dimIndex] = points[nodeNumber].Get(dimIndex);
			}
		}

		std::array<int, Space::Dimension> imins;
		std::array<int, Space::Dimension> imaxs;

		for (auto dimIndex = 0; dimIndex < Space::Dimension; ++dimIndex)
		{
			auto minRefValue = (mins[dimIndex] - meshAABB.boxPoint1.Get(dimIndex)) / (meshAABB.boxPoint2.Get(dimIndex) - meshAABB.boxPoint1.Get(dimIndex));
			imins[dimIndex] = int(minRefValue * (stride.Get(dimIndex) - 1));

			auto maxRefValue = (maxs[dimIndex] - meshAABB.boxPoint1.Get(dimIndex)) / (meshAABB.boxPoint2.Get(dimIndex) - meshAABB.boxPoint1.Get(dimIndex));
			imaxs[dimIndex] = int(maxRefValue * (stride.Get(dimIndex) - 1)) + 1;

			adjustValueToRange<int>(imins[dimIndex], 0, stride.Get(dimIndex) - 1);
			adjustValueToRange<int>(imaxs[dimIndex], 0, stride.Get(dimIndex) - 1);
		}

		Vector pointToCheck;
		IndexVector iv;

		std::function<void(int)> checkPoints = [&](int dimIndex)
		{
			if (dimIndex == Space::Dimension)
			{
				if (PointInCell<Scalar>(points, pointToCheck))
				{
					materialInfoIndices[idx<Space>(iv, stride)] = cellMaterialInfos[cellIndex];
				}
			} else
			{ 
				for (int coord = imins[dimIndex]; coord <= imaxs[dimIndex]; ++coord)
				{
					pointToCheck[dimIndex] = meshAABB.boxPoint1.Get(dimIndex) +
						(meshAABB.boxPoint2.Get(dimIndex) - meshAABB.boxPoint1.Get(dimIndex)) * coord / Scalar(stride.Get(dimIndex) - 1);
					iv[dimIndex] = coord;
					checkPoints(dimIndex + 1);
				}
			}
		};

		checkPoints(0);
	}

	return materialInfoIndices;
}


int main()
{
	SPACE_TYPEDEFS;

	const IndexVector stride = { 61, 61, 61 };
	AABB meshAABB;

	auto materialIndices = loadMesh("./meshes/test", stride, meshAABB);

	std::vector<material_info<Scalar>> mis;
	mis.push_back({ 267, 267, 1000.0, 920.0, 0.591, 2.22, 334000.0, 4200.0, 1100.0 });
	mis.push_back({ 270, 270, 1000.0, 920.0, 0.591, 2.22, 334000.0, 4200.0, 2100.0 });
	//mis.push_back(material_info(273, 273, 1000.0, 920.0, 0.591, 2.22, 334000.0, 4200.0, 2100.0));
	//mis.push_back(material_info(50, 50, 1.3, 1.3, 0.0243, 0.0243, 1e9, 1005.0, 1005.0));
	//mis.push_back(material_info(5000, 5000, 2837.0, 2837.0, 1.4, 1.4, 1e9, 1480, 1480));
	

	std::function<material_info<Scalar>(const typename Space::Vector&)> materialParamsGetter = [&mis, &meshAABB, &stride, &materialIndices](const typename Space::Vector& point)
	{ 
		Space::iVector iv;

		for (auto dimIndex = 0; dimIndex < Space::Dimension; ++dimIndex)
		{
			auto refValue = (point.Get(dimIndex) - meshAABB.boxPoint1.Get(dimIndex)) / (meshAABB.boxPoint2.Get(dimIndex) - meshAABB.boxPoint1.Get(dimIndex));
			iv[dimIndex] = int(refValue * (stride.Get(dimIndex) - 1));
			adjustValueToRange<int>(iv[dimIndex], 0, stride.Get(dimIndex) - 1);
		}

		auto index = idx<Space>(iv, stride);
		return mis.at(materialIndices.at(index));
	};

	const auto meshSize = meshAABB.Size();
	auto mesh = Mesh<Space>(stride, meshAABB, materialParamsGetter);

	// set initial temperature
	std::vector<Scalar> td(stride.GetVolume(), 0); // todo: default value from config 
					
	for (const auto& node : AllNodes<Space>(stride))
	{
		const Vector step = meshSize.ComponentDivide(stride - Vector::one());
		const Vector point = meshAABB.boxPoint1 + step.ComponentMultiply(node);
		//Scalar x = meshAABB.boxPoint1.x + meshSize.x * n / (stride.x - 1);
		//Scalar y = meshAABB.boxPoint1.y + meshSize.y * i / (stride.y - 1);
		//Scalar z = meshAABB.boxPoint1.z + meshSize.z * l / (stride.z - 1);
//				if (x > 2 && x < 4 && y > 2 && y < 4 && z > 2 && z < 4)
		if (materialIndices[idx<Space>(node, stride)] == 0)
		{
			td[idx<Space>(node, stride)] = 272.0;
		}
		else
		{
			td[idx<Space>(node, stride)] = 303.0;
		}
	}

	auto sol = solver3d<Space>(mesh, td.data());
	for (int i = 0; i < 3000; ++i)
	{
		if (i%10 == 0)
		{
			std::stringstream ss;
			ss << i;
			std::cout << "step: " << i << std::endl;
			sol.save_to_vtk("out/result_" + ss.str() + ".vtk");
		}
		sol.step(1000);	
	}

	return 0;
}

/*
int main()
{
	SPACE_TYPEDEFS;

	const IndexVector stride = { 241, 121 };
	AABB meshAABB;

	std::vector<material_info<Scalar>> mis;
	mis.push_back({ 267, 267, 1000.0, 920.0, 0.591, 2.22, 334000.0, 4200.0, 1100.0 });
	mis.push_back({ 270, 270, 1000.0, 920.0, 0.591, 2.22, 334000.0, 4200.0, 2100.0 });
	mis.push_back({ 273, 273, 1000.0, 920.0, 0.591, 2.22, 334000.0, 4200.0, 2100.0 });
	mis.push_back({ 50, 50, 1.3, 1.3, 0.0243, 0.0243, 1e9, 1005.0, 1005.0 });
	mis.push_back({ 5000, 5000, 2837.0, 2837.0, 1.4, 1.4, 1e9, 1480, 1480 });
	
	auto materialIndices = loadMesh("./meshes/test2", stride, meshAABB);

	std::function<material_info<Scalar>(const Vector&)> materialParamsGetter = 
		[&mis, &meshAABB, &stride, &materialIndices](const Vector& point)
	{
		iVector iv;

		for (auto dimIndex = 0; dimIndex < Space::Dimension; ++dimIndex)
		{
			auto refValue = (point.Get(dimIndex) - meshAABB.boxPoint1.Get(dimIndex)) / (meshAABB.boxPoint2.Get(dimIndex) - meshAABB.boxPoint1.Get(dimIndex));
			iv[dimIndex] = int(refValue * (stride.Get(dimIndex) - 1));
			adjustValueToRange<int>(iv[dimIndex], 0, stride.Get(dimIndex) - 1);
		}

		auto index = idx<Space>(iv, stride);
		return mis.at(materialIndices.at(index));
	};

	const auto meshSize = meshAABB.Size();
	auto mesh = Mesh<Space>(stride, meshAABB, materialParamsGetter);

	std::vector<Scalar> td(stride.GetVolume());

	for (IndexType n = 0; n < stride.x; ++n)
		for (IndexType i = 0; i < stride.y; ++i)
		{
			//double x = meshAABB.boxPoint1.x + meshSize.x * n / (stride.x - 1);
			//double y = meshAABB.boxPoint1.y + meshSize.y * n / (stride.y - 1);
			if (materialIndices[idx<Space>({ n, i }, stride)] == 0)
			{
				td[idx<Space>({ n, i }, stride)] = 263.0;
			} else
			{
				td[idx<Space>({ n, i }, stride)] = 273.0;
			}
		}

	auto sol = solver2d<Space>(mesh, td.data());
	for (int i = 0; i < 300; ++i)
	{
		if (i%10 == 0)
		{
			std::stringstream ss;
			ss << i;
			std::cout << "step: " << i << std::endl;
			sol.save_to_vtk("out/result_" + ss.str() + ".vtk");
		}
		sol.step(300);	
	}
	return 0;
}
*/
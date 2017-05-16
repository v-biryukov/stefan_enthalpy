#include <algorithm>

template <typename Space>
void MeshIO<Space>::Load(const std::string& fileName, IO::FileType fileType)
{
  std::fstream file;

  switch (fileType)
  {
    case IO::Ascii: file.open(fileName.c_str(), std::fstream::in); break;
    case IO::Binary: file.open(fileName.c_str(), std::fstream::in | std::fstream::binary); break;
  }
  
  if (file.fail())
  {
    std::cerr << "Can`t open file " << fileName << " for loading" << std::endl;
    throw;
  }

  Load(file, fileType);
  file.close();
}

template <typename Space>
typename Space::IndexType MeshIO<Space>::GetCellsCount() const
{
  return indices.size() / Space::NodesPerCell;
}

template <typename Space>
typename Space::Vector MeshIO<Space>::GetCellCenter(IndexType cellIndex) const
{
  Vector center = Vector::zero();
  for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
  {
    center += vertices[indices[cellIndex * Space::NodesPerCell + nodeNumber]];
  }
  center /= Scalar(Space::NodesPerCell);
  return center;
}

template <typename Space>
void MeshIO<Space>::Load(std::fstream& file, IO::FileType fileType)
{
  switch (fileType)
  {
    case IO::Ascii:
    {
      // vertices
      IndexType verticesCount;
      file >> verticesCount;
      vertices.resize(verticesCount);

      for (IndexType vertexIndex = 0; vertexIndex < vertices.size(); ++vertexIndex)
      {
        for (IndexType componentNumber = 0; componentNumber < Space::Dimension; ++componentNumber)
        {
          file >> vertices[vertexIndex][componentNumber];
        }
      }

      // indices
      IndexType indicesCount;
      file >> indicesCount;
      indices.resize(indicesCount);
      for (IndexType index = 0; index < indices.size(); ++index)
      {
        file >> indices[index];
      }

      // detectors
      IndexType detectorsCount;
      file >> detectorsCount;
      detectorsPositions.resize(detectorsCount);
      for (IndexType detectorIndex = 0; detectorIndex < detectorsPositions.size(); ++detectorIndex)
      {
        for (IndexType componentNumber = 0; componentNumber < Space::Dimension; ++componentNumber)
        {
          file >> detectorsPositions[detectorIndex][componentNumber];
        }
      }
    } break;

    case IO::Binary:
    {
      // vertices
      IndexType verticesCount;
      IO::Read(file, verticesCount);
      vertices.resize(verticesCount);
      IO::Read(file, vertices.data(), verticesCount);

      // indices
      IndexType indicesCount;
      IO::Read(file, indicesCount);
      indices.resize(indicesCount);
      IO::Read(file, indices.data(), indicesCount);

      // detectors
      IndexType detectorsCount;
      IO::Read(file, detectorsCount);
      detectorsPositions.resize(detectorsCount);
      IO::Read(file, detectorsPositions.data(), detectorsCount);
    } break;
  }
}

template <typename Space>
typename Space::AABB getAABB(const MeshIO<Space>& meshIO)
{
	SPACE_TYPEDEFS
	Space::AABB aabb;
	for (auto cellIndex = 0; cellIndex < meshIO.GetCellsCount(); ++cellIndex)
	{
		for (auto nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
		{
			auto nodeIndex = meshIO.indices[cellIndex * Space::NodesPerCell + nodeNumber];
			aabb.Expand(meshIO.vertices[nodeIndex]);
		}
	}
	return aabb;
}
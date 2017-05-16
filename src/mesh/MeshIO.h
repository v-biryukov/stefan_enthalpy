#pragma once

#include "MeshIOBase.h"
#include "../utils/Utils.h"

template <typename Space>
struct MeshIO
{
  SPACE_TYPEDEFS
  virtual ~MeshIO() = default;
  virtual void Load(const std::string& fileName, IO::FileType fileType = IO::Binary);

  // file should be opened
  void Load(std::fstream& file, IO::FileType fileType);

  std::vector<Vector>    vertices;
  std::vector<IndexType> indices;

  std::vector<Vector> detectorsPositions;

  IndexType GetCellsCount() const;
  Vector GetCellCenter(IndexType cellIndex) const;
};

template <typename Space>
typename Space::AABB getAABB(const MeshIO<Space>& meshIO);

#include "MeshIO.inl"
#include "MeshIO2.inl"
#include "MeshIO3.inl"

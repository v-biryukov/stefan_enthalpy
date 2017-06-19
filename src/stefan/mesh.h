#include <array>
#include <vector>
#include <numeric>
#include "material_info.h"

template <int Dims>
class Mesh
{
	std::array<int, Dims> nums_;
	std::array<double, Dims> lengths_;
	std::array<double, Dims> steps_;
	int number_of_nodes_;

	std::vector<MaterialInfo> material_infos_;
	std::vector<int> material_indexes_;

public:

	Mesh() = default;
	Mesh( const Mesh & ) = default;
	Mesh(std::array<int, Dims> nums, std::array<double, Dims> lengths,
		 const std::vector<MaterialInfo>& material_infos,
		 const std::vector<int>& material_indexes_)
		: nums_(nums), lengths_(lengths), steps_(steps_),
		material_infos_(material_infos), material_indexes_(material_indexes_)
	{
		for (int i = 0; i < Dims; ++i)
			steps_[i] = lengths_[i] / (nums_[i] - 1);

		number_of_nodes_ = std::accumulate(nums_.begin(), nums_.end(), 1, std::multiplies<int>());
	}

	int GetGlobalId(std::array<int, Dims> ids) const;

	const MaterialInfo& GetMaterialInfo(std::array<int, Dims> ids) const
	{
		return material_infos_[material_indexes_[GetGlobalId(ids)]];
	}

	const MaterialInfo& GetMaterialInfo(int global_id) const
	{
		return material_infos_[material_indexes_[global_id]];
	}

	const int GetMaterialIndex(std::array<int, Dims> ids) const
	{
		return material_indexes_[GetGlobalId(ids)];
	}

	const int GetMaterialIndex(int global_id) const
	{
		return material_indexes_[global_id];
	}

	inline std::array<double, Dims> GetLengths() const {return lengths_;}
	inline std::array<double, Dims> GetSteps() const {return steps_;}
	inline std::array<int, Dims> GetNums() const {return nums_;}
	inline int GetNumberOfNodes() const {return number_of_nodes_;}

	//const std::vector<int> GetMaterialIndexes() const {return material_indexes_;}
	//const std::vector<MaterialInfo> GetMaterialInfos() const {return material_infos_;}
};

template<>
int Mesh<2>::GetGlobalId(std::array<int, 2> ids) const
{
	return ids[0] + ids[1]*nums_[0];
}

template<>
int Mesh<3>::GetGlobalId(std::array<int, 3> ids) const
{
	return ids[0] + ids[1]*nums_[0] + ids[2]*nums_[0]*nums_[1];
}

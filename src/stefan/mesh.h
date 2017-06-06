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
	Mesh(std::array<int, Dims> nums, std::array<double, Dims> lengths,
		 std::vector<MaterialInfo> material_infos,
		 std::vector<int> material_indexes_) 
		: nums_(nums), lengths_(lengths), steps_(steps_), 
		material_infos_(material_infos), material_indexes_(material_indexes_)
	{
		for (int i = 0; i < Dims; ++i)
			steps_[i] = lengths_[i] / (nums_[i] - 1);
		
		number_of_nodes_ = std::accumulate(nums_.begin(), nums_.end(), 1, std::multiplies<int>());
	}

	int GetGlobalId(std::array<int, Dims> ids);

	MaterialInfo & GetMaterialInfo(std::array<int, Dims> ids)
	{
		return material_infos_[material_indexes_[GetGlobalId(ids)]];
	}

	MaterialInfo & GetMaterialInfo(int global_id)
	{
		return material_infos_[material_indexes_[global_id]];
	}

	inline std::array<double, Dims> GetLengths() {return lengths_;}
	inline std::array<double, Dims> GetSteps() {return steps_;}
	inline std::array<int, Dims> GetNums() {return nums_;}
	inline int GetNumberOfNodes() {return number_of_nodes_;}
};

template<>
int Mesh<2>::GetGlobalId(std::array<int, 2> ids)
{
	return ids[0] + ids[1]*nums_[0];
}

template<>
int Mesh<3>::GetGlobalId(std::array<int, 3> ids)
{
	return ids[0] + ids[1]*nums_[0] + ids[2]*nums_[0]*nums_[1];
}

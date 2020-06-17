#pragma once
#include <array>
#include <vector>
#include <numeric>
#include "material_info.h"


bool EndsWith(const std::string& fullString, const std::string& ending) 
{
    if (fullString.length() >= ending.length()) 
    {
        return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
    } 
    else 
    {
        return false;
    }
}


template <int Dims>
class Mesh
{
private:
    std::array<int, Dims> nums;
    std::array<double, Dims> lengths;
    std::array<double, Dims> steps;
    int number_of_nodes;

    std::vector<MaterialInfo> material_infos;
    std::vector<int> material_indexes;

    int number_of_submeshes;

    void CalculateNumberOfSubmeshes()
    {
        std::set<int> submeshTypes;
        for (int index : material_indexes)
            submeshTypes.insert(index);
        number_of_submeshes = submeshTypes.size();
    }

public:

    Mesh() = default;
    Mesh(const Mesh&) = default;
    Mesh(std::array<int, Dims> nums, std::array<double, Dims> lengths, const std::vector<MaterialInfo>& material_infos, const std::vector<int>& material_indexes)
        : nums(nums), lengths(lengths), material_infos(material_infos), material_indexes(material_indexes)
    {
        for (int i = 0; i < Dims; ++i)
            steps[i] = lengths[i] / (nums[i] - 1);
        number_of_nodes = std::accumulate(nums.begin(), nums.end(), 1, std::multiplies<int>());

        CalculateNumberOfSubmeshes();
    }

    Mesh(std::string mesh_file, const std::vector<MaterialInfo>& material_infos)
        : material_infos(material_infos)
    {
        ReadMeshFile(mesh_file);
        CalculateNumberOfSubmeshes();
    }

    std::array<int, Dims> GetIds(int global_id) const;
    int GetGlobalId(std::array<int, Dims> ids) const;

    const MaterialInfo& GetMaterialInfo(std::array<int, Dims> ids) const
    {
        return material_infos[material_indexes[GetGlobalId(ids)]];
    }

    const MaterialInfo& GetMaterialInfo(int global_id) const
    {
        return material_infos[material_indexes[global_id]];
    }

    int GetMaterialIndex(std::array<int, Dims> ids) const
    {
        return material_indexes[GetGlobalId(ids)];
    }

    int GetMaterialIndex(int global_id) const
    {
        return material_indexes[global_id];
    }

    void SetMaterialIndex(std::array<int, Dims> ids, int new_material_index)
    {
        material_indexes[GetGlobalId(ids)] = new_material_index;
    }

    void SetMaterialIndex(int global_id, int new_material_index)
    {
        material_indexes[global_id] = new_material_index;
    }

    inline std::array<double, Dims> GetLengths() const {return lengths;}
    inline std::array<double, Dims> GetSteps() const {return steps;}
    inline std::array<int, Dims> GetNums() const {return nums;}
    inline int GetNumberOfNodes() const {return number_of_nodes;}
    inline int GetNumberOfSubmeshes() const {return number_of_submeshes;}


    void ReadMeshFileAscii(std::string filename);
    void ReadMeshFileBinary(std::string filename);
    void ReadMeshFile(std::string filename);
};


template <int Dims>
void Mesh<Dims>::ReadMeshFileAscii(std::string filename)
{
    std::ifstream file(filename, std::ios::in);
    if (file.is_open())
    {
        for (int d = 0; d < Dims; d++)
        {
            file >> nums[d];
        }
        for (int d = 0; d < Dims; d++)
        {
            file >> lengths[d];
        }
        int num_of_nodes = std::accumulate(nums.begin(), nums.end(), 1, std::multiplies<int>());
        material_indexes.resize(num_of_nodes);
        for (int i = 0; i < num_of_nodes; i++)
        {
            file >> material_indexes[i];
        }
        file.close();
    }
    else 
    {
        std::cerr << "Error: Unable to open ascii mesh file " << filename << "\n";
        std::exit(EXIT_FAILURE);
    }
}


template <int Dims>
void Mesh<Dims>::ReadMeshFileBinary(std::string filename)
{
    std::ifstream file(filename, std::ios::binary);
    if (file.is_open())
    {
        file.read((char*)nums.data(), Dims * sizeof(int));
        file.read((char*)lengths.data(), Dims * sizeof(double));
        int num_of_nodes = std::accumulate(nums.begin(), nums.end(), 1, std::multiplies<int>());
        material_indexes.resize(num_of_nodes);
        file.read((char*)material_indexes.data(), num_of_nodes * sizeof(int));
        file.close();
    }
    else 
    {
        std::cerr << "Error: Unable to open binary mesh file " << filename << "\n";
        std::exit(EXIT_FAILURE);
    }
}


template <int Dims>
void Mesh<Dims>::ReadMeshFile(std::string filename)
{
    if (EndsWith(filename, ".mesh"))
    {
        ReadMeshFileAscii(filename);
    }
    else if (EndsWith(filename, ".bmesh"))
    {
        ReadMeshFileBinary(filename);
    }
    else
    {
        std::cerr << "Error: Mesh file" << filename << " has wrong extension\n";
        std::exit(EXIT_FAILURE);
    }

    for (int i = 0; i < Dims; ++i)
        steps[i] = lengths[i] / (nums[i] - 1);
    number_of_nodes = std::accumulate(nums.begin(), nums.end(), 1, std::multiplies<int>());
}



template<>
std::array<int, 2> Mesh<2>::GetIds(int global_id) const
{

    return {global_id % nums[0], global_id / nums[0]};
}

template<>
std::array<int, 3> Mesh<3>::GetIds(int global_id) const
{
    return {global_id % nums[0], (global_id / nums[0]) % nums[1],  global_id / (nums[0] * nums[1])};
}

template<>
int Mesh<2>::GetGlobalId(std::array<int, 2> ids) const
{
    return ids[0] + ids[1]*nums[0];
}

template<>
int Mesh<3>::GetGlobalId(std::array<int, 3> ids) const
{
    return ids[0] + ids[1]*nums[0] + ids[2]*nums[0]*nums[1];
}

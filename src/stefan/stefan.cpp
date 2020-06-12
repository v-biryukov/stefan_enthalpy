#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <math.h>
#include <vector>
#include <string>
#include <iomanip>

#include "../utils/Utils.h"

#include "material_info.h"
#include "solver.h"
#include "settings_parser.h"



void PrintHelpInfo()
{
    std::cout << "Help info" << std::endl;
}

template <int Dims>
void ReadInitialStateFile(std::string filename, std::array<int, Dims> nums, std::vector<double>& initial_temperatures)
{
    std::ifstream file (filename, std::ios::binary);
    if (file.is_open())
    {
        int num_of_nodes = std::accumulate(nums.begin(), nums.end(), 1, std::multiplies<int>());
        initial_temperatures.resize(num_of_nodes);
        file.read ((char*)initial_temperatures.data(), num_of_nodes * sizeof(double));
        file.close();
    }
    else
    {
        std::cerr << "Error: Unable to open initial temperatures file \n";
        std::exit(EXIT_FAILURE);
    }
}

template <int Dims>
int Run(const Settings& settings)
{
    std::vector<MaterialInfo> material_infos;
    for (auto mi : settings.mesh_settings.medium_params_info)
    {
        material_infos.push_back(MaterialInfo(mi.T_phase, mi.T_phase, mi.density_L, mi.density_S,
                                              mi.thermal_conductivity_L, mi.thermal_conductivity_S, 
                                              mi.specific_heat_fusion, mi.specific_heat_capacity_L, mi.specific_heat_capacity_S));
    }

    Mesh<Dims> mesh = Mesh<Dims>(settings.mesh_settings.mesh_file, material_infos);

    std::vector<double> initial_temperatures;
    if (settings.task_settings.initial_state_settings.type == Settings::TaskSettings::InitialStateSettings::InitialStateSettingsType::PerNode)
    {
        ReadInitialStateFile<Dims>(settings.task_settings.initial_state_settings.initial_state_data_file, mesh.GetNums(), initial_temperatures);
    }
    else if (settings.task_settings.initial_state_settings.type == Settings::TaskSettings::InitialStateSettings::InitialStateSettingsType::PerSubmesh)
    {
        for (int i = 0; i < mesh.GetNumberOfNodes(); i++)
            initial_temperatures.push_back(settings.task_settings.initial_state_settings.initial_temperatures_by_submesh[mesh.GetMaterialIndex(i)]);
    }

    std::array<std::array<double, 3>, 2*Dims> boundary_conditions;
    for (int i = 0; i < boundary_conditions.size(); ++i)
        boundary_conditions[i] = settings.mesh_settings.boundary_settings_info[i].params;

    Solver<Dims> sol = Solver<Dims>(mesh, boundary_conditions, initial_temperatures);
    sol.AddFields(settings);
    for (int time_step_num = 0; time_step_num < settings.task_settings.number_of_steps; ++time_step_num)
    {
        if (time_step_num % settings.snapshot_settings_info.period_frames == 0)
        {
            std::stringstream ss;
            ss << std::setfill('0') << std::setw(5) << time_step_num;
            std::string filename = settings.snapshot_settings_info.snaphot_data_path;
            ReplaceSubstring(filename, "<step>", ss.str());
            sol.SaveToVtk(filename, settings.snapshot_settings_info);
        }
        std::cout << "step: " << time_step_num << std::endl;

        sol.Step(settings.task_settings.time_step);
    }
    return 0;
}


int main(int argc, char** argv)
{
    if (argc != 2)
    {
        PrintHelpInfo();
        std::exit(EXIT_FAILURE);
    }
    else if (!std::string(argv[1]).compare("-h") || !std::string(argv[1]).compare("--help"))
    {
        PrintHelpInfo();
        std::exit(EXIT_SUCCESS);
    }
    else
    {
        Settings settings = ParseFile(std::string(argv[1]));
        if (settings.dims_count == 2)
        {
            Run<2>(settings);
        }
        else if (settings.dims_count == 3)
        {
            Run<3>(settings);
        }
    }
    return 0;
}

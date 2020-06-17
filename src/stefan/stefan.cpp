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
int Run(const Settings<Dims>& settings)
{
    Solver<Dims> sol = Solver<Dims>(settings);
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
        std::string config_filename = argv[1];
        int dims_count = ReadDimsCountFromConfig(config_filename);

        if (dims_count == 2)
        {
            Settings<2> settings;
            settings.ParseFile(config_filename);
            Run<2>(settings);
        }
        else if (dims_count == 3)
        {
            Settings<3> settings;
            settings.ParseFile(config_filename);
            Run<3>(settings);
        }
        
    }
    return 0;
}

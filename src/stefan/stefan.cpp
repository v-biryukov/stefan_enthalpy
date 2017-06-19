#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <math.h>
#include <vector>
#include <string>

#include "material_info.h"
#include "solver.h"
#include "settings_parser.h"



void PrintHelpInfo()
{
	std::cout << "Help info" << std::endl;
}


template <int Dims>
void ReadMeshFile(std::string file_name, std::array<int, Dims> & nums, 
	std::array<double, Dims> & lengths, std::vector<int> & material_indexes)
{
	std::ifstream file (file_name, std::ios::binary);
	if (file.is_open())
	{
		file.read ((char*)nums.data(), Dims * sizeof(int));
		file.read ((char*)lengths.data(), Dims * sizeof(double));
		int num_of_nodes = std::accumulate(nums.begin(), nums.end(), 1, std::multiplies<int>());
		material_indexes.resize(num_of_nodes);
		file.read ((char*)material_indexes.data(), num_of_nodes * sizeof(int));
		file.close();
	}
	else 
	{
		std::cerr << "Error: Unable to open mesh file \n";
		std::exit(EXIT_FAILURE);
	}
}

template <int Dims>
void ReadInitialStateFile(std::string file_name, std::array<int, Dims> nums, std::vector<double> & initial_temperatures)
{
	std::ifstream file (file_name, std::ios::binary);
	if (file.is_open())
	{
		int num_of_nodes = std::accumulate(nums.begin(), nums.end(), 1, std::multiplies<int>());
		initial_temperatures.reserve(num_of_nodes);
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
int Run(const Settings & settings)
{
	std::array<int, Dims> nums;
	std::array<double, Dims> lengths;
	std::vector<int> material_indexes;
	ReadMeshFile<Dims>(settings.mesh_settings.mesh_file, nums, lengths, material_indexes);

	std::vector<MaterialInfo> material_infos;
	for (auto mi : settings.mesh_settings.medium_params_info)
	{
		material_infos.push_back(MaterialInfo(mi.T_phase-2.5, mi.T_phase+2.5, mi.rho_L, mi.rho_S,
			mi.thermal_conductivity_L, mi.thermal_conductivity_S, 
			mi.specific_heat_fusion, mi.specific_heat_capacity_L, mi.specific_heat_capacity_S));
	}

	Mesh<Dims> mesh = Mesh<Dims>(nums, lengths, material_infos, material_indexes);

	std::vector<double> initial_temperatures;
	if (settings.task_settings.initial_state_settings.type == Settings::TaskSettings::InitialStateSettings::InitialStateSettingsType::PerNode)
	{
		ReadInitialStateFile<Dims>(settings.task_settings.initial_state_settings.initial_state_data_file, nums, initial_temperatures);
	}
	else if (settings.task_settings.initial_state_settings.type == Settings::TaskSettings::InitialStateSettings::InitialStateSettingsType::PerSubmesh)
	{
		for (auto el : material_indexes)
			initial_temperatures.push_back(settings.task_settings.initial_state_settings.initial_temperatures_by_submesh[el]);
	}

	std::array<std::array<double, 3>, 2*Dims> boundary_conditions;
		for (int i = 0; i < boundary_conditions.size(); ++i)
		boundary_conditions[i] = settings.mesh_settings.boundary_settings_info[i].params;
	Solver<Dims> sol = Solver<Dims>(mesh, boundary_conditions, initial_temperatures);
	for (int time_step_num = 0; time_step_num < settings.task_settings.number_of_steps; ++time_step_num)
	{
		if (time_step_num % settings.snapshot_settings_info.period_frames == 0)
		{
			std::stringstream ss;
			ss << time_step_num;
			sol.SaveToVtk("out/result_" + ss.str() + ".vtk");
		}
		std::cout << "step: " << time_step_num << std::endl;

		sol.Step(settings.task_settings.time_step);
	}
	return 0;
}


int main(int argc, char ** argv)
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


/*

int main()
{
	int num_x = 241;
	int num_y = 121;
	double Lx = 12.0;
	double Ly = 6.0;

	std::vector<material_info> mis;
	mis.push_back(material_info(267, 267, 1000.0, 920.0, 0.591, 2.22, 334000.0, 4200.0, 1100.0));
	mis.push_back(material_info(270, 270, 1000.0, 920.0, 0.591, 2.22, 334000.0, 4200.0, 2100.0));
	mis.push_back(material_info(273, 273, 1000.0, 920.0, 0.591, 2.22, 334000.0, 4200.0, 2100.0));
	mis.push_back(material_info(50, 50, 1.3, 1.3, 0.0243, 0.0243, 1e9, 1005.0, 1005.0));
	mis.push_back(material_info(5000, 5000, 2837.0, 2837.0, 1.4, 1.4, 1e9, 1480, 1480));
	

	std::function<int(double, double)> mat_idx = []( double x, double y )
	{ 
		if (y <= 1)
			return 4;
		else if (y < 2)
			return 2;
		else if (y < 3)
			return 1;
		else if (y < 4 || ((y<4.5) && (x>4) && (x<8)))
			return 0;
		else
			return 3;
	};

	mesh2d mesh = mesh2d(num_x, num_y, Lx, Ly, mis, mat_idx);

	double * td = new double [num_x * num_y];
	for (int n = 0; n < num_x; ++n)
		for (int i = 0; i < num_y; ++i)
		{
			double x = n*Lx/(num_x-1);
			double y = i*Ly/(num_y-1);
			if ((x > 4 && x < 8 && y > 1 && y < 4.5))
			{
				td[n + i*num_x] = 263.0;
			}
			else if (y <= 1)
			{
				td[n + i*num_x] = 273.0;
			}
			else
			{
				td[n + i*num_x] = 303.0;
			}
		}

	solver2d sol = solver2d(mesh, td);
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
	delete [] td;
	return 0;
}

	/*
	std::ifstream file(file_name, std::ios::binary);
	for (int i = 0; i < Dims; ++i)
		file >> nums[i];
	for (int i = 0; i < Dims; ++i)
		file >> lengths[i];

	mat_idx.reserve(fileSize);
	
	std::copy(std::istream_iterator<int>(file),
		std::istream_iterator<int>(),
		std::back_inserter(mat_idx));
	*/

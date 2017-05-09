#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <math.h>
#include <vector>
#include <functional>

#include "material_info.hpp"
#include "solver2d.hpp"
#include "solver3d.hpp"


using Scalar = double;

int main()
{
	int num_x = 61;
	int num_y = 61;
	int num_z = 61;
	Scalar Lx = 6.0;
	Scalar Ly = 6.0;
	Scalar Lz = 6.0;

	std::vector<material_info<Scalar>> mis;
	mis.push_back(material_info<Scalar>(267, 267, 1000.0, 920.0, 0.591, 2.22, 334000.0, 4200.0, 1100.0));
	//mis.push_back(material_info(270, 270, 1000.0, 920.0, 0.591, 2.22, 334000.0, 4200.0, 2100.0));
	//mis.push_back(material_info(273, 273, 1000.0, 920.0, 0.591, 2.22, 334000.0, 4200.0, 2100.0));
	//mis.push_back(material_info(50, 50, 1.3, 1.3, 0.0243, 0.0243, 1e9, 1005.0, 1005.0));
	//mis.push_back(material_info(5000, 5000, 2837.0, 2837.0, 1.4, 1.4, 1e9, 1480, 1480));
	

	std::function<int(Scalar, Scalar, Scalar)> mat_idx = []( Scalar /*x*/, Scalar /*y*/, Scalar /*z*/ )
	{ 
		return 0;
	};

	auto mesh = mesh3d<Scalar>(num_x, num_y, num_z, Lx, Ly, Lz, mis, mat_idx);

	Scalar* td = new Scalar [num_x * num_y * num_z];
	for (int n = 0; n < num_x; ++n)
		for (int i = 0; i < num_y; ++i)
			for (int l = 0; l < num_y; ++l)
			{
				Scalar x = n*Lx/(num_x-1);
				Scalar y = i*Ly/(num_y-1);
				Scalar z = l*Lz/(num_z-1);
				if (x > 2 && x < 4 && y > 2 && y < 4 && z > 2 && z < 4)
				{
					td[n + i*num_x + l*num_x*num_y] = 272.0;
				}
				else
				{
					td[n + i*num_x + l*num_x*num_y] = 303.0;
				}
			}

	auto sol = solver3d<Scalar>(mesh, td);
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
	delete [] td;
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
*/
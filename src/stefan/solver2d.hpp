#pragma once
#include <math.h>
#include <vector>
#include <functional>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "SolverBase.h"

template <typename Space>
class solver2d : public SolverBase<Space>
{
public:
	SPACE_TYPEDEFS;

	using SolverBase<Space>::step;

	void iterate_tridiagonal(int axis,
		const std::vector<Scalar>& enthalpy_data,
		const std::vector<Scalar>& current_iteration_data,
		std::vector<Scalar>& next_iteration_data,
		Scalar dt) override
	{
		int idx1_p, idx1_c, idx1_m;
		int idx2_p, idx2_c, idx2_m;
		int num_1, num_2;
		size_t ix, iy;
		Scalar h1, h2;

		if (axis == 0)
		{
			idx1_p = 1;
			idx1_c = 0;
			idx1_m = -1;
			idx2_p = mesh.getStride().x;
			idx2_c = 0;
			idx2_m = -(int)(mesh.getStride().x);
			num_1 = mesh.getStride().x; num_2 = mesh.getStride().y;
			h1 = mesh.get_h().x; h2 = mesh.get_h().y;
		}
		else /*if (axis == 1)*/
		{
			idx1_p = mesh.getStride().x;
			idx1_c = 0;
			idx1_m = -(int)(mesh.getStride().x);
			idx2_p = 1;
			idx2_c = 0;
			idx2_m = -1;
			num_1 = mesh.getStride().y; num_2 = mesh.getStride().x;
			h1 =  mesh.get_h().y; h2 = mesh.get_h().x;
		}

		material_info<Scalar> mi, min1, mip1, min2, mip2;
		
		for (int i2 = 1; i2 < num_2-1; ++i2)
		{
			a[0] = 0.0; b[0] = -1.0; c[0] = 1.0; f[0] = 0.0;
			a[num_1-1] = -1.0; b[num_1-1] = 1.0; c[num_1-1] = 0.0; f[num_1-1] = 0.0;

			for (size_t i1 = 1; i1 < num_1 - 1; ++i1)
			{
				if (axis == 0)
				{
					ix = i1; iy = i2;
					mi = mesh.get_material({ ix, iy });
					min1 = mesh.get_material({ ix - 1, iy });
					mip1 = mesh.get_material({ ix + 1, iy });
					min2 = mesh.get_material({ ix, iy - 1 });
					mip2 = mesh.get_material({ ix, iy + 1 });
				}
				else
				{
					ix = i2; iy = i1;
					mi = mesh.get_material({ ix, iy });
					min1 = mesh.get_material({ ix, iy - 1 });
					mip1 = mesh.get_material({ ix, iy + 1 });
					min2 = mesh.get_material({ ix - 1, iy });
					mip2 = mesh.get_material({ ix + 1, iy });
				}
				
				const Scalar* E_iter = current_iteration_data.data() + id({ ix, iy });
				const Scalar* E_step = enthalpy_data.data() + id({ ix, iy });

				Scalar kp = (mip1.get_thermal_conductivity_by_E(E_iter[idx1_p]) + mi.get_thermal_conductivity_by_E(E_iter[idx1_c])) / 2.0;
				Scalar km = (mi.get_thermal_conductivity_by_E(E_iter[idx1_c]) + min1.get_thermal_conductivity_by_E(E_iter[idx1_m])) / 2.0;

				Scalar coef = dt / (2.0 * h1 * h1);
				a[i1] = -coef * min1.get_alpha(E_iter[idx1_m]) * km;
				b[i1] = 1.0 + coef * mi.get_alpha(E_iter[idx1_c]) * (km + kp);
				c[i1] = -coef * mip1.get_alpha(E_iter[idx1_p]) * kp;
				f[i1] = E_step[idx1_c] + coef * (min1.get_beta(E_iter[idx1_m])*km - mi.get_beta(E_iter[idx1_c])*(km+kp) + mip1.get_beta(E_iter[idx1_p])*kp);

				Scalar k2p = (mip2.get_thermal_conductivity_by_E(E_step[idx2_p]) + mi.get_thermal_conductivity_by_E(E_step[idx2_c])) / 2.0;
				Scalar k2n = (min2.get_thermal_conductivity_by_E(E_step[idx2_m]) + mi.get_thermal_conductivity_by_E(E_step[idx2_c])) / 2.0;
				Scalar T2p = mip2.get_T_by_enthalpy(E_step[idx2_p]) - mi.get_T_by_enthalpy(E_step[idx2_c]);
				Scalar T2n = min2.get_T_by_enthalpy(E_step[idx2_m]) - mi.get_T_by_enthalpy(E_step[idx2_c]);
				f[i1] += (k2p * T2p + k2n * T2n) * dt/(2 * h2 * h2);
			}

			solve_tridiagonal(num_1, a, b, c, f, x);
			for (int k = 0; k < num_1; ++k)
				next_iteration_data[idx1_p * k + idx2_p * i2] = x[k];
		}
	}

public:
	solver2d(const Mesh<Space>& mesh, const Scalar* temperature_data) : SolverBase(mesh, temperature_data)
	{
	}

	void save_to_vtk(const std::string& name)
	{
		int num_x = mesh.getStride().x;
		int num_y = mesh.getStride().y;

		std::ofstream vtk_file(name.c_str());
		vtk_file << "# vtk DataFile Version 3.0\nVx data\nASCII\n\n";
		vtk_file << "DATASET STRUCTURED_POINTS\nDIMENSIONS " << num_x << " " << num_y << " 1" << "\n";
		vtk_file << "ORIGIN 0 0 0\n" << "SPACING " << mesh.get_h().x << " " << mesh.get_h().y << " 1" << "\n";

		vtk_file << "\nPOINT_DATA " << mesh.get_number_of_nodes() << "\n";
		vtk_file << "SCALARS E FLOAT\n";
		vtk_file << "LOOKUP_TABLE default\n";
		
		for (int i = 0; i < mesh.get_number_of_nodes(); ++i)
		{
			vtk_file <<  enthalpy_data[i] << "\n";
		}
		vtk_file << "SCALARS T FLOAT\n";
		vtk_file << "LOOKUP_TABLE default\n";
		for (size_t j = 0; j < num_y; j++)
			for (size_t i = 0; i < num_x; i++)
			{
				auto mi = mesh.get_material({ i, j });
				vtk_file <<  mi.get_T_by_enthalpy(enthalpy_data[i + j*num_x]) << "\n";
			}
		vtk_file << "SCALARS K FLOAT\n";
		vtk_file << "LOOKUP_TABLE default\n";
		for (size_t j = 0; j < num_y; j++)
			for (size_t i = 0; i < num_x; i++)
			{
				auto mi = mesh.get_material({ i, j });
				vtk_file <<  mi.get_thermal_conductivity_by_E(enthalpy_data[i + j*num_x]) << "\n";
			}
		vtk_file << "SCALARS state FLOAT\n";
		vtk_file << "LOOKUP_TABLE default\n";
		for (size_t j = 0; j < num_y; j++)
			for (size_t i = 0; i < num_x; i++)
			{
				auto mi = mesh.get_material({ i, j });
				Scalar T = mi.get_T_by_enthalpy(enthalpy_data[i + j*num_x]);
				Scalar state;
				if (T <= mi.T1)
					state = 0.0;
				else if (T <= mi.T2)
					state = (T-mi.T1)/(mi.T2-mi.T1);
				else
					state = 1.0;
				vtk_file <<  state << "\n";
			}
	}
};

#pragma once
#include <math.h>
#include <vector>
#include <functional>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <assert.h>
#include "SolverBase.h"

template <typename Space>
class solver3d : public SolverBase<Space>
{
public:
	SPACE_TYPEDEFS;

	void iterate_tridiagonal(int axis,
		const std::vector<Scalar>& enthalpy_data,
		const std::vector<Scalar>& current_iteration_data,
		std::vector<Scalar>& next_iteration_data,
		Scalar dt) override
	{
		int idx1_p, idx1_m;
		int idx2_p, idx2_m;
		int idx3_p, idx3_m;
		int num_1, num_2, num_3;
		size_t ix, iy, iz;
		Scalar h1, h2, h3;

		int idx_c = 0;
		if (axis == 0)
		{
			idx1_p = 1;
			idx1_m = -1;
			idx2_p = mesh.getStride().x;
			idx2_m = -idx2_p;
			idx3_p = mesh.getStride().x * mesh.getStride().y;
			idx3_m = -idx3_p;
			num_1 = mesh.getStride().x; num_2 = mesh.getStride().y; num_3 = mesh.getStride().z;
			h1 = mesh.get_h().x; h2 = mesh.get_h().y; h3 = mesh.get_h().z;
		}
		else if (axis == 1)
		{
			idx1_p = mesh.getStride().x;
			idx1_m = -idx1_p;
			idx2_p = mesh.getStride().x * mesh.getStride().y;
			idx2_m = -idx2_p;
			idx3_p = 1;
			idx3_m = -1;
			num_1 = mesh.getStride().y; num_2 = mesh.getStride().z; num_3 = mesh.getStride().x;
			h1 =  mesh.get_h().y; h2 = mesh.get_h().z;  h3 = mesh.get_h().x;
		}
		else if (axis == 2)
		{
			idx1_p = mesh.getStride().x * mesh.getStride().y;
			idx1_m = -idx1_p;
			idx2_p = 1;
			idx2_m = -1;
			idx3_p = mesh.getStride().x;
			idx3_m = -idx3_p;
			num_1 = mesh.getStride().z; num_2 = mesh.getStride().x; num_3 = mesh.getStride().y;
			h1 =  mesh.get_h().z; h2 = mesh.get_h().x;  h3 = mesh.get_h().y;
		}
		else
		{
			assert(false);
		}

		material_info<Scalar> mi, min1, mip1, min2, mip2, min3, mip3;
		
		for (int i3 = 1; i3 < num_3-1; ++i3)
		{
			for (int i2 = 1; i2 < num_2-1; ++i2)
			{
				a[0] = 0.0; b[0] = -1.0; c[0] = 1.0; f[0] = 0.0;
				a[num_1-1] = -1.0; b[num_1-1] = 1.0; c[num_1-1] = 0.0; f[num_1-1] = 0.0;

				for (int i1 = 1; i1 < num_1-1; ++i1)
				{
					if (axis == 0)
					{
						ix = i1; iy = i2; iz = i3;
						mi = mesh.get_material({ ix, iy, iz });
						min1 = mesh.get_material({ ix - 1, iy, iz });
						mip1 = mesh.get_material({ ix + 1, iy, iz });
						min2 = mesh.get_material({ ix, iy - 1, iz });
						mip2 = mesh.get_material({ ix, iy + 1, iz });
						min3 = mesh.get_material({ ix, iy, iz - 1 });
						mip3 = mesh.get_material({ ix, iy, iz + 1 });
					}
					else if (axis == 1)
					{
						ix = i3; iy = i1; iz = i2;
						mi = mesh.get_material({ ix, iy, iz });
						min1 = mesh.get_material({ ix, iy - 1, iz });
						mip1 = mesh.get_material({ ix, iy + 1, iz });
						min2 = mesh.get_material({ ix, iy, iz - 1 });
						mip2 = mesh.get_material({ ix, iy, iz + 1 });
						min3 = mesh.get_material({ ix - 1, iy, iz });
						mip3 = mesh.get_material({ ix + 1, iy, iz });
					}
					else
					{
						ix = i2; iy = i3; iz = i1;
						mi = mesh.get_material({ ix, iy, iz });
						min1 = mesh.get_material({ ix, iy, iz - 1 });
						mip1 = mesh.get_material({ ix, iy, iz + 1 });
						min2 = mesh.get_material({ ix - 1, iy, iz });
						mip2 = mesh.get_material({ ix + 1, iy, iz });
						min3 = mesh.get_material({ ix, iy - 1, iz });
						mip3 = mesh.get_material({ ix, iy + 1, iz });
					}
					
					const Scalar* E_iter = current_iteration_data.data() + id({ ix, iy, iz });
					const Scalar* E_step = enthalpy_data.data() + id({ ix, iy, iz });

					Scalar kp = (mip1.get_thermal_conductivity_by_E(E_iter[idx1_p]) + mi.get_thermal_conductivity_by_E(E_iter[idx_c])) / 2.0;
					Scalar km = (mi.get_thermal_conductivity_by_E(E_iter[idx_c]) + min1.get_thermal_conductivity_by_E(E_iter[idx1_m])) / 2.0;

					Scalar coef = dt / (3.0 * h1 * h1);
					a[i1] = -coef * min1.get_alpha(E_iter[idx1_m]) * km;
					b[i1] = 1.0 + coef * mi.get_alpha(E_iter[idx_c]) * (km + kp);
					c[i1] = -coef * mip1.get_alpha(E_iter[idx1_p]) * kp;
					f[i1] = E_step[idx_c] + coef * (min1.get_beta(E_iter[idx1_m])*km - mi.get_beta(E_iter[idx_c])*(km+kp) + mip1.get_beta(E_iter[idx1_p])*kp);

					Scalar k2p = (mip2.get_thermal_conductivity_by_E(E_step[idx2_p]) + mi.get_thermal_conductivity_by_E(E_step[idx_c])) / 2.0;
					Scalar k2n = (min2.get_thermal_conductivity_by_E(E_step[idx2_m]) + mi.get_thermal_conductivity_by_E(E_step[idx_c])) / 2.0;
					Scalar T2p = mip2.get_T_by_enthalpy(E_step[idx2_p]) - mi.get_T_by_enthalpy(E_step[idx_c]);
					Scalar T2n = min2.get_T_by_enthalpy(E_step[idx2_m]) - mi.get_T_by_enthalpy(E_step[idx_c]);
					f[i1] += k2p * T2p * dt/(3*h2*h2) + k2n * T2n * dt/(3*h2*h2);

					Scalar k3p = (mip3.get_thermal_conductivity_by_E(E_step[idx3_p]) + mi.get_thermal_conductivity_by_E(E_step[idx_c])) / 2.0;
					Scalar k3n = (min3.get_thermal_conductivity_by_E(E_step[idx3_m]) + mi.get_thermal_conductivity_by_E(E_step[idx_c])) / 2.0;
					Scalar T3p = mip3.get_T_by_enthalpy(E_step[idx3_p]) - mi.get_T_by_enthalpy(E_step[idx_c]);
					Scalar T3n = min3.get_T_by_enthalpy(E_step[idx3_m]) - mi.get_T_by_enthalpy(E_step[idx_c]);
					f[i1] += k3p * T3p * dt/(3*h3*h3) + k3n * T3n * dt/(3*h3*h3);
				}

				solve_tridiagonal(num_1, a, b, c, f, x);
				for (int k = 0; k < num_1; ++k)
					next_iteration_data[idx1_p * k + idx2_p * i2 + idx3_p * i3] = x[k];
			}
		}
	}

	using SolverBase<Space>::step;

public:
	solver3d(const Mesh<Space>& mesh, const Scalar* temperature_data) : 
		SolverBase(mesh, temperature_data)
	{
	}

	

	void save_to_vtk(const std::string& name)
	{
		int num_x = mesh.getStride().x;
		int num_y = mesh.getStride().y;
		int num_z = mesh.getStride().z;

		std::ofstream vtk_file(name.c_str());
		vtk_file << "# vtk DataFile Version 3.0\nVx data\nASCII\n\n";
		vtk_file << "DATASET STRUCTURED_POINTS\nDIMENSIONS " << num_x << " " << num_y << " " << num_z << "\n";
		vtk_file << "ORIGIN 0 0 0\n" << "SPACING " << mesh.get_h().x << " " << mesh.get_h().y << " " << mesh.get_h().z << "\n";

		vtk_file << "\nPOINT_DATA " << mesh.get_number_of_nodes() << "\n";
		vtk_file << "SCALARS E FLOAT\n";
		vtk_file << "LOOKUP_TABLE default\n";

		for (IndexType l = 0; l < num_z; l++)
			for (IndexType i = 0; i < num_y; i++)
				for (IndexType n = 0; n < num_x; n++)
				{
					vtk_file << enthalpy_data[id({ n, i, l })] << "\n";
				}
		vtk_file << "SCALARS T FLOAT\n";
		vtk_file << "LOOKUP_TABLE default\n";
		for (IndexType l = 0; l < num_z; l++)
			for (IndexType i = 0; i < num_y; i++)
				for (IndexType n = 0; n < num_x; n++)
				{
					auto mi = mesh.get_material({ n, i, l });
					vtk_file << mi.get_T_by_enthalpy(enthalpy_data[id({ n, i, l })]) << "\n";
				}
		vtk_file << "SCALARS K FLOAT\n";
		vtk_file << "LOOKUP_TABLE default\n";
		for (size_t l = 0; l < num_z; l++)
			for (size_t i = 0; i < num_y; i++)
				for (size_t n = 0; n < num_x; n++)
				{
					auto mi = mesh.get_material({ n, i, l });
					vtk_file << mi.get_thermal_conductivity_by_E(enthalpy_data[id({ n, i, l })]) << "\n";
				}
		vtk_file << "SCALARS state FLOAT\n";
		vtk_file << "LOOKUP_TABLE default\n";
		for (size_t l = 0; l < num_z; l++)
			for (size_t i = 0; i < num_y; i++)
				for (size_t n = 0; n < num_x; n++)
				{
					auto mi = mesh.get_material({ n, i, l });
					Scalar T = mi.get_T_by_enthalpy(enthalpy_data[id({ n, i, l })]);
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

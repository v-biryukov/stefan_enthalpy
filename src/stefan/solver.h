#pragma once
#include <algorithm>
#include <cassert>
#include <string>
#include <iostream>
#include <fstream>
#include <omp.h>

#include "mesh.h"
#include "material_info.h"
#include "settings_parser.h"

template <int Dims>
class Solver
{
	Mesh<Dims> mesh;

	std::vector<double> enthalpy_data;
	std::vector<double> current_iteration_data;
	std::vector<double> next_iteration_data;

	// Boundary conditions
	std::array<std::array<double, 3>, 2 * Dims> boundary_conditions;
	const double boundary_condition_eps = 1e-15;

	void SetEnthalpyByTData(const std::vector<double>& temperature_data)
	{
		for (int i = 0; i < mesh.GetNumberOfNodes(); ++i)
			enthalpy_data[i] = mesh.GetMaterialInfo(i).GetEnthalpyByT(temperature_data[i]);
	}

	static void SolveTridiagonal (int n, std::vector<double>& a,
		std::vector<double>& b, std::vector<double>& c,
		std::vector<double>& f, std::vector<double>& x)
	{
		double m;
		for (int i = 1; i < n; i++)
		{
			m = a[i]/b[i-1];

			b[i] = b[i] - m*c[i-1];
			f[i] = f[i] - m*f[i-1];
		}

		x[n-1] = f[n-1]/b[n-1];
		for (int i = n - 2; i >= 0; i--)
		{
			x[i]=(f[i]-c[i]*x[i+1])/b[i];
		}
	}

	static double CalculateLinfNorm(const std::vector<double>& a, const std::vector<double>& b)
	{
		assert(a.size() == b.size());
		double Linf_norm = 0.0;
		for (int i = 0; i < a.size(); ++i)
		{
			double current_norm = fabs(a[i] - b[i]);
			if (current_norm > Linf_norm)
				Linf_norm = current_norm;
		}
		return Linf_norm;
	}

	static double CalculateLinfNorm(const std::vector<double>& a)
	{
		double Linf_norm = 0.0;
		for (int i = 0; i < a.size(); ++i)
		{
			double current_norm = fabs(a[i]);
			if (current_norm > Linf_norm)
				Linf_norm = current_norm;
		}
		return Linf_norm;
	}

	static double CalculateL2Norm(const std::vector<double>& a, const std::vector<double>& b)
	{
		assert(a.size() == b.size());
		double L2_norm = 0.0;
		for (int i = 0; i < a.size(); ++i)
		{
			L2_norm += (a[i] - b[i])*(a[i] - b[i]);
		}
		L2_norm = sqrt(L2_norm);
		return L2_norm;
	}

	static double CalculateL2Norm(const std::vector<double>& a)
	{
		double L2_norm = 0.0;
		for (int i = 0; i < a.size(); ++i)
		{
			L2_norm += a[i]*a[i];
		}
		L2_norm = sqrt(L2_norm);
		return L2_norm;
	}

	static void Assign(std::vector<double>& to, const std::vector<double>& from)
	{
		assert(to.size() == from.size());
		std::copy_n(from.begin(), to.size(), to.begin());
	}

	void Step(int axis, double dt)
	{
		double Linf_norm = 1.0;
		Assign(next_iteration_data, enthalpy_data);
		while (Linf_norm > 1e-3)
		{
			Assign(current_iteration_data, next_iteration_data);
			IterateTridiagonal(axis, enthalpy_data, current_iteration_data, next_iteration_data, dt);
			Linf_norm = CalculateLinfNorm(current_iteration_data, next_iteration_data) / CalculateLinfNorm(current_iteration_data);
		}
		Assign(enthalpy_data, next_iteration_data);
	}


	void SetBoundaryConditions(int axis, std::array<int, Dims> index,
		std::vector<double> & a, std::vector<double> & b,
		std::vector<double> & c, std::vector<double> & f)
	{
		int num = mesh.GetNums()[axis];
		double h = mesh.GetSteps()[axis];
		std::rotate(index.begin(), index.begin() + axis, index.end());
		index[axis] = 0;

		double alpha = mesh.GetMaterialInfo(index).GetAlpha(enthalpy_data[mesh.GetGlobalId(index)]);
		if (fabs(alpha) < boundary_condition_eps)
			alpha = boundary_condition_eps;
		double beta = mesh.GetMaterialInfo(index).GetBeta(enthalpy_data[mesh.GetGlobalId(index)]);

		a[0] = 0.0;
		b[0] = alpha*(boundary_conditions[2*axis][0] - boundary_conditions[2*axis][1]/h);
		c[0] = boundary_conditions[2*axis][1] * alpha / h ;
		f[0] = (boundary_conditions[2*axis][2] - boundary_conditions[2*axis][0] * beta);
		index[axis] = num-1;
		alpha = mesh.GetMaterialInfo(index).GetAlpha(enthalpy_data[mesh.GetGlobalId(index)]);
		if (fabs(alpha) < boundary_condition_eps)
			alpha = boundary_condition_eps;
		beta = mesh.GetMaterialInfo(index).GetBeta(enthalpy_data[mesh.GetGlobalId(index)]);
		a[num-1] = -boundary_conditions[2*axis+1][1] * alpha / h ;
		b[num-1] = alpha*(boundary_conditions[2*axis+1][0] + boundary_conditions[2*axis+1][1]/h);
		c[num-1] = 0.0;
		f[num-1] = (boundary_conditions[2*axis+1][2] - boundary_conditions[2*axis+1][0] * beta);
	}

	const bool isInfoToWrite(int global_id, std::vector<int>& stride) const
	{
		std::array<int, Dims> ids = mesh.GetIds(global_id);
		bool write_info = true;
		for (int dim = 0; dim < Dims; ++dim)
		{
			write_info = write_info && (ids[dim] % stride[dim] == 0);
		}
		return write_info;
	}


	void IterateTridiagonal(int axis, 
		const std::vector<double>& enthalpy_data,
		const std::vector<double>& current_iteration_data,
		std::vector<double>& next_iteration_data,
		double dt);

public:



	Solver(Mesh<Dims> & mesh, std::array<std::array<double, 3>, 2*Dims> boundary_conditions, const std::vector<double> & temperature_data)
		: mesh(mesh), boundary_conditions(boundary_conditions)
	{
		enthalpy_data.resize(mesh.GetNumberOfNodes());
		current_iteration_data.resize(mesh.GetNumberOfNodes());
		next_iteration_data.resize(mesh.GetNumberOfNodes());

		SetEnthalpyByTData(temperature_data);
		auto nums = mesh.GetNums();
		auto max_n = *std::max_element(nums.begin(), nums.end());
	}

	void Step(double dt)
	{
		for (int i = 0; i < Dims; ++i)
			Step(i, dt);
	}

	void SaveToVtk(const std::string name, const Settings::SnapshotSettings & snapshot_settings) const
	{
		std::array<int, Dims> nums = mesh.GetNums();
		std::vector<int> stride = snapshot_settings.stride;
		int number_of_nodes_to_save = 1;
		for (int dim = 0; dim < Dims; ++dim)
			number_of_nodes_to_save *= ((nums[dim] - 1) / stride[dim] + 1);

		std::ofstream vtk_file(name.c_str());
		vtk_file << "# vtk DataFile Version 3.0\nStefan data\nASCII\n\n";
		vtk_file << "DATASET STRUCTURED_POINTS\nDIMENSIONS " << (nums[0] - 1) / stride[0] + 1 << " " << (nums[1] - 1) / stride[1] + 1;
		if (Dims == 2)
			vtk_file << " " << 1 << "\n";
		else if (Dims == 3)
			vtk_file << " " << (nums[2] - 1)  / stride[2] + 1 << "\n";

		vtk_file << "ORIGIN 0 0 0\n" << "SPACING " << mesh.GetSteps()[0] * stride[0] << " " << mesh.GetSteps()[1]  * stride[1];
		if (Dims == 2)
			vtk_file << " " << 1 << "\n";
		else if (Dims == 3)
			vtk_file << " " << mesh.GetSteps()[2] * stride[2] << "\n";
		vtk_file << "\nPOINT_DATA " << number_of_nodes_to_save << "\n";
		if (snapshot_settings.write_enthalpy)
		{
			vtk_file << "SCALARS E FLOAT\n";
			vtk_file << "LOOKUP_TABLE default\n";
			for(std::size_t i = 0; i < enthalpy_data.size() ; ++i)
			{
				if (isInfoToWrite(i, stride))
					vtk_file <<  enthalpy_data[i] << "\n";
			}
		}
		if (snapshot_settings.write_temperature)
		{
			vtk_file << "SCALARS T FLOAT\n";
			vtk_file << "LOOKUP_TABLE default\n";
			for(std::size_t i = 0; i < enthalpy_data.size(); ++i)
			{
				if (isInfoToWrite(i, stride))
					vtk_file <<  mesh.GetMaterialInfo(i).GetTByEnthalpy(enthalpy_data[i]) << "\n";
			}
		}

		if (snapshot_settings.write_thermal_elasticity)
		{
			vtk_file << "SCALARS K FLOAT\n";
			vtk_file << "LOOKUP_TABLE default\n";
			for(std::size_t i = 0; i < enthalpy_data.size(); ++i)
			{
				if (isInfoToWrite(i, stride))
					vtk_file <<  mesh.GetMaterialInfo(i).GetThermalConductivityByE(enthalpy_data[i]) << "\n";
			}
		}

		if (snapshot_settings.write_state)
		{
			vtk_file << "SCALARS state FLOAT\n";
			vtk_file << "LOOKUP_TABLE default\n";
			for(std::size_t i = 0; i < enthalpy_data.size(); ++i)
			{
				if (isInfoToWrite(i, stride))
				{
					const MaterialInfo& mi = mesh.GetMaterialInfo(i);
					double T = mi.GetTByEnthalpy(enthalpy_data[i]);
					double state;
					if (T <= mi.T1)
						state = 0.0;
					else if (T <= mi.T2)
						state = (T-mi.T1)/(mi.T2-mi.T1);
					else
						state = 1.0;
					vtk_file <<  state << "\n";
				}
			}
		}
		vtk_file.close();
	}

};



template<>
void Solver<2>::IterateTridiagonal(int axis,
	const std::vector<double>& enthalpy_data,
	const std::vector<double>& current_iteration_data,
	std::vector<double>& next_iteration_data,
	double dt)
{
	#pragma omp parallel firstprivate(axis, dt)
	{
		int idx1_p,  idx1_m;
		int idx2_p,  idx2_m;
		int num_1, num_2;
		int ix, iy;
		double h1, h2;
		int ort1x, ort1y, ort2x, ort2y;

		if (axis == 0)
		{
			idx1_p = 1;
			idx1_m = -1;
			idx2_p = mesh.GetNums()[0];
			idx2_m = -mesh.GetNums()[0];
			num_1 = mesh.GetNums()[0]; num_2 = mesh.GetNums()[1];
			h1 = mesh.GetSteps()[0]; h2 = mesh.GetSteps()[1];
			ort1x = 1; ort1y = 0; ort2x = 0; ort2y = 1;
		}
		else if (axis == 1)
		{
			idx1_p = mesh.GetNums()[0];
			idx1_m = -mesh.GetNums()[0];
			idx2_p = 1;
			idx2_m = -1;
			num_1 = mesh.GetNums()[1]; num_2 = mesh.GetNums()[0];
			h1 =  mesh.GetSteps()[1]; h2 = mesh.GetSteps()[0];
			ort1x = 0; ort1y = 1; ort2x = 1; ort2y = 0;
		}

		auto nums = mesh.GetNums();
		auto max_n = *std::max_element(nums.begin(), nums.end());
		std::vector<double> a(max_n);
		std::vector<double> b(max_n);
		std::vector<double> c(max_n);
		std::vector<double> f(max_n);
		std::vector<double> x(max_n);

		const double coef1 = dt / (2.0 * h1 * h1);
		const double coef2 = dt / (2.0 * h2 * h2);

		#pragma omp for schedule(dynamic)
		for (int i2 = 1; i2 < num_2-1; ++i2)
		{
			SetBoundaryConditions(axis, {0, i2}, a, b, c, f);

			for (int i1 = 1; i1 < num_1-1; ++i1)
			{
				if (axis == 0)
				{
					ix = i1; iy = i2;
				}
				else
				{
					ix = i2; iy = i1;
				}
				const MaterialInfo & mi = mesh.GetMaterialInfo({ix, iy});
				const MaterialInfo & min1 = mesh.GetMaterialInfo({ix-ort1x, iy-ort1y});
				const MaterialInfo & mip1 = mesh.GetMaterialInfo({ix+ort1x, iy+ort1y});
				const MaterialInfo & min2 = mesh.GetMaterialInfo({ix-ort2x, iy-ort2y});
				const MaterialInfo & mip2 = mesh.GetMaterialInfo({ix+ort2x, iy+ort2y});

				const double * E_iter = current_iteration_data.data() + mesh.GetGlobalId({ix, iy});
				const double * E_step = enthalpy_data.data() + mesh.GetGlobalId({ix, iy});
				const double central_thermal_conductivity_iter = mi.GetThermalConductivityByE(E_iter[0]);
				const double central_thermal_conductivity_step = mi.GetThermalConductivityByE(E_step[0]);
				const double kp = (mip1.GetThermalConductivityByE(E_iter[idx1_p]) + central_thermal_conductivity_iter) / 2.0;
				const double km = (central_thermal_conductivity_iter + min1.GetThermalConductivityByE(E_iter[idx1_m])) / 2.0;

				a[i1] = -coef1 * min1.GetAlpha(E_iter[idx1_m]) * km;
				b[i1] = 1.0 + coef1 * mi.GetAlpha(E_iter[0]) * (km + kp);
				c[i1] = -coef1 * mip1.GetAlpha(E_iter[idx1_p]) * kp;
				const double k2p = (mip2.GetThermalConductivityByE(E_step[idx2_p]) + central_thermal_conductivity_step) / 2.0;
				const double k2n = (min2.GetThermalConductivityByE(E_step[idx2_m]) + central_thermal_conductivity_step) / 2.0;
				const double T2p = mip2.GetTByEnthalpy(E_step[idx2_p]) - mi.GetTByEnthalpy(E_step[0]);
				const double T2n = min2.GetTByEnthalpy(E_step[idx2_m]) - mi.GetTByEnthalpy(E_step[0]);
				f[i1] = E_step[0] + coef1 * (min1.GetBeta(E_iter[idx1_m])*km - mi.GetBeta(E_iter[0])*(km+kp) + mip1.GetBeta(E_iter[idx1_p])*kp) + (k2p * T2p + k2n * T2n) * coef2;
			}
			SolveTridiagonal(num_1, a, b, c, f, x);
			for (int k = 0; k < num_1; ++k)
				next_iteration_data[idx1_p * k + idx2_p * i2] = x[k];
		}
	}
}


template<>
void Solver<3>::IterateTridiagonal(int axis, 
	const std::vector<double>& enthalpy_data,
	const std::vector<double>& current_iteration_data,
	std::vector<double>& next_iteration_data,
	double dt)
{
	#pragma omp parallel firstprivate(axis, dt)
	{
		int idx1_p, idx1_m;
		int idx2_p, idx2_m;
		int idx3_p, idx3_m;
		int num_1, num_2, num_3;
		int ix, iy, iz;
		double h1, h2, h3;
		int ort1x, ort1y, ort1z, ort2x, ort2y, ort2z, ort3x, ort3y, ort3z;

		if (axis == 0)
		{
			idx1_p = 1;
			idx1_m = -1;
			idx2_p = mesh.GetNums()[0];
			idx2_m = -mesh.GetNums()[0];
			idx3_p = mesh.GetNums()[0]*mesh.GetNums()[1];
			idx3_m = -mesh.GetNums()[0]*mesh.GetNums()[1];
			num_1 = mesh.GetNums()[0]; num_2 = mesh.GetNums()[1]; num_3 = mesh.GetNums()[2];
			h1 = mesh.GetSteps()[0]; h2 = mesh.GetSteps()[1]; h3 = mesh.GetSteps()[2];
			ort1x = 1; ort1y = 0; ort1z = 0;
			ort2x = 0; ort2y = 1; ort2z = 0;
			ort3x = 0; ort3y = 0; ort3z = 1;
		}
		else if (axis == 1)
		{
			idx1_p = mesh.GetNums()[0];
			idx1_m = -mesh.GetNums()[0];
			idx2_p = mesh.GetNums()[0]*mesh.GetNums()[1];
			idx2_m = -mesh.GetNums()[0]*mesh.GetNums()[1];
			idx3_p = 1;
			idx3_m = -1;
			num_1 = mesh.GetNums()[1]; num_2 = mesh.GetNums()[2];; num_3 = mesh.GetNums()[0];
			h1 =  mesh.GetSteps()[1]; h2 = mesh.GetSteps()[2];  h3 = mesh.GetSteps()[0];
			ort1x = 0; ort1y = 1; ort1z = 0;
			ort2x = 0; ort2y = 0; ort2z = 1;
			ort3x = 1; ort3y = 0; ort3z = 0;
		}
		else if (axis == 2)
		{
			idx1_p = mesh.GetNums()[0]*mesh.GetNums()[1];
			idx1_m = -mesh.GetNums()[0]*mesh.GetNums()[1];
			idx2_p = 1;
			idx2_m = -1;
			idx3_p = mesh.GetNums()[0];
			idx3_m = -mesh.GetNums()[0];
			num_1 = mesh.GetNums()[2]; num_2 = mesh.GetNums()[0]; num_3 = mesh.GetNums()[1];
			h1 =  mesh.GetSteps()[2]; h2 = mesh.GetSteps()[0];  h3 = mesh.GetSteps()[1];
			ort1x = 0; ort1y = 0; ort1z = 1;
			ort2x = 1; ort2y = 0; ort2z = 0;
			ort3x = 0; ort3y = 1; ort3z = 0;
		}

		auto nums = mesh.GetNums();
		auto max_n = *std::max_element(nums.begin(), nums.end());
		std::vector<double> a(max_n);
		std::vector<double> b(max_n);
		std::vector<double> c(max_n);
		std::vector<double> f(max_n);
		std::vector<double> x(max_n);

		MaterialInfo mi, min1, mip1, min2, mip2, min3, mip3;
		const double coef1 = dt / (3.0 * h1 * h1);
		const double coef2 = dt / (3.0 * h2 * h2);
		const double coef3 = dt / (3.0 * h3 * h3);

		#pragma omp for schedule(dynamic)
		for (int i3 = 1; i3 < num_3-1; ++i3)
		{
			for (int i2 = 1; i2 < num_2-1; ++i2)
			{
				SetBoundaryConditions(axis, {0, i2, i3}, a, b, c, f);
				for (int i1 = 1; i1 < num_1-1; ++i1)
				{
					if (axis == 0)
					{
						ix = i1; iy = i2; iz = i3;
					}
					else if (axis == 1)
					{
						ix = i3; iy = i1; iz = i2;
					}
					else
					{
						ix = i2; iy = i3; iz = i1;
					}
					const MaterialInfo & mi = mesh.GetMaterialInfo({ix, iy});
					const MaterialInfo & min1 = mesh.GetMaterialInfo({ix-ort1x, iy-ort1y, iz-ort1z});
					const MaterialInfo & mip1 = mesh.GetMaterialInfo({ix+ort1x, iy+ort1y, iz+ort1z});
					const MaterialInfo & min2 = mesh.GetMaterialInfo({ix-ort2x, iy-ort2y, iz-ort2z});
					const MaterialInfo & mip2 = mesh.GetMaterialInfo({ix+ort2x, iy+ort2y, iz+ort2z});
					const MaterialInfo & min3 = mesh.GetMaterialInfo({ix-ort3x, iy-ort3y, iz-ort3z});
					const MaterialInfo & mip3 = mesh.GetMaterialInfo({ix+ort3x, iy+ort3y, iz+ort3z});

					const double * E_iter = current_iteration_data.data() + mesh.GetGlobalId({ix, iy, iz});
					const double * E_step = enthalpy_data.data() + mesh.GetGlobalId({ix, iy, iz});
					const double central_thermal_conductivity_iter = mi.GetThermalConductivityByE(E_iter[0]);
					const double central_thermal_conductivity_step = mi.GetThermalConductivityByE(E_step[0]);
					const double central_T = mi.GetTByEnthalpy(E_step[0]);

					const double kp = (mip1.GetThermalConductivityByE(E_iter[idx1_p]) + central_thermal_conductivity_iter) / 2.0;
					const double km = (central_thermal_conductivity_iter + min1.GetThermalConductivityByE(E_iter[idx1_m])) / 2.0;

					a[i1] = -coef1 * min1.GetAlpha(E_iter[idx1_m]) * km;
					b[i1] = 1.0 + coef1 * mi.GetAlpha(E_iter[0]) * (km + kp);
					c[i1] = -coef1 * mip1.GetAlpha(E_iter[idx1_p]) * kp;

					double k2p = (mip2.GetThermalConductivityByE(E_step[idx2_p]) + central_thermal_conductivity_step) / 2.0;
					double k2n = (min2.GetThermalConductivityByE(E_step[idx2_m]) + central_thermal_conductivity_step) / 2.0;
					double T2p = mip2.GetTByEnthalpy(E_step[idx2_p]) - central_T;
					double T2n = min2.GetTByEnthalpy(E_step[idx2_m]) - central_T;
					double k3p = (mip3.GetThermalConductivityByE(E_step[idx3_p]) + central_thermal_conductivity_step) / 2.0;
					double k3n = (min3.GetThermalConductivityByE(E_step[idx3_m]) + central_thermal_conductivity_step) / 2.0;
					double T3p = mip3.GetTByEnthalpy(E_step[idx3_p]) - central_T;
					double T3n = min3.GetTByEnthalpy(E_step[idx3_m]) - central_T;

					f[i1] = E_step[0] + coef1 * (min1.GetBeta(E_iter[idx1_m])*km - mi.GetBeta(E_iter[0])*(km+kp) + mip1.GetBeta(E_iter[idx1_p])*kp)
							+ (k2p * T2p + k2n * T2n) * coef2
							+ (k3p * T3p + k3n * T3n) * coef3;
				}

				SolveTridiagonal(num_1, a, b, c, f, x);
				for (int k = 0; k < num_1; ++k)
						next_iteration_data[idx1_p * k + idx2_p * i2 + idx3_p * i3] = x[k];
			}
		}
	}
}


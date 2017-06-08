#include <algorithm>
#include <cassert>
#include <string>
#include <iostream>
#include <fstream>

#include "mesh.h"
#include "material_info.h"

template <int Dims>
class Solver
{
	Mesh<Dims> mesh;

	std::vector<double> enthalpy_data;
	std::vector<double> current_iteration_data;
	std::vector<double> next_iteration_data;

	// Auxiliary vectors for tridiagonal method
	std::vector<double> a, b, c, f, x;

	// Boundary conditions
	std::array<std::array<double, 3>, 2*Dims> boundary_conditions;

	void SetEnthalpyByTData(const std::vector<double> & temperature_data)
	{
		for (int i = 0; i < mesh.GetNumberOfNodes(); ++i)
			enthalpy_data[i] = mesh.GetMaterialInfo(i).GetEnthalpyByT(temperature_data[i]);
	}

	void SolveTridiagonal (int n, std::vector<double> & a, 
		std::vector<double> & b, std::vector<double> & c, 
		std::vector<double> & f, std::vector<double> & x)
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
			x[i]=(f[i]-c[i]*x[i+1])/b[i];
	}

	double CalculateLinfNorm(const std::vector<double>& a, const std::vector<double>& b)
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

	double CalculateLinfNorm(const std::vector<double>& a)
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
		double beta = mesh.GetMaterialInfo(index).GetBeta(enthalpy_data[mesh.GetGlobalId(index)]);
		a[0] = 0.0; 
		b[0] = alpha*(boundary_conditions[2*axis][0] - boundary_conditions[2*axis][1]/h); 
		c[0] = boundary_conditions[2*axis][1] * alpha / h; 
		f[0] = boundary_conditions[2*axis][2] - boundary_conditions[2*axis][0] * beta;
		index[axis] = num-1;
		alpha = mesh.GetMaterialInfo(index).GetAlpha(enthalpy_data[mesh.GetGlobalId(index)]);
		beta = mesh.GetMaterialInfo(index).GetBeta(enthalpy_data[mesh.GetGlobalId(index)]);
		a[num-1] = -boundary_conditions[2*axis+1][1] * alpha / h; 
		b[num-1] = alpha*(boundary_conditions[2*axis+1][0] + boundary_conditions[2*axis+1][1]/h); 
		c[num-1] = 0.0; 
		f[num-1] = boundary_conditions[2*axis+1][2] - boundary_conditions[2*axis+1][0] * beta;
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
		a.resize(max_n);
		b.resize(max_n);
		c.resize(max_n);
		f.resize(max_n);
		x.resize(max_n);
	}

	void Step(double dt)
	{
		for (int i = 0; i < Dims; ++i)
			Step(i, dt);
	}

	void SaveToVtk(std::string name)
	{
		std::ofstream vtk_file(name.c_str());
		vtk_file << "# vtk DataFile Version 3.0\nVx data\nASCII\n\n";
		vtk_file << "DATASET STRUCTURED_POINTS\nDIMENSIONS " << mesh.GetNums()[0] << " " << mesh.GetNums()[1];
		if (Dims == 2)
			vtk_file << " " << 1 << "\n";
		else if (Dims == 3)
			vtk_file << " " << mesh.GetNums()[2] << "\n";

		vtk_file << "ORIGIN 0 0 0\n" << "SPACING " << mesh.GetSteps()[0] << " " << mesh.GetSteps()[1];
		if (Dims == 2)
			vtk_file << " " << 1 << "\n";
		else if (Dims == 3)
			vtk_file << " " << mesh.GetSteps()[2] << "\n";

		vtk_file << "\nPOINT_DATA " << mesh.GetNumberOfNodes() << "\n";
		vtk_file << "SCALARS E FLOAT\n";
		vtk_file << "LOOKUP_TABLE default\n";
		
		for(std::size_t i = 0; i < enthalpy_data.size(); ++i)
		{
	    	vtk_file <<  enthalpy_data[i] << "\n";
		}

		vtk_file << "SCALARS T FLOAT\n";
		vtk_file << "LOOKUP_TABLE default\n";
		for(std::size_t i = 0; i < enthalpy_data.size(); ++i)
		{
			vtk_file <<  mesh.GetMaterialInfo(i).GetTByEnthalpy(enthalpy_data[i]) << "\n";
		}

		vtk_file << "SCALARS K FLOAT\n";
		vtk_file << "LOOKUP_TABLE default\n";
		for(std::size_t i = 0; i < enthalpy_data.size(); ++i)
		{
			vtk_file <<  mesh.GetMaterialInfo(i).GetThermalConductivityByE(enthalpy_data[i]) << "\n";
		}

		vtk_file << "SCALARS state FLOAT\n";
		vtk_file << "LOOKUP_TABLE default\n";
		for(std::size_t i = 0; i < enthalpy_data.size(); ++i)
		{
			MaterialInfo & mi = mesh.GetMaterialInfo(i);
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

};



template<>
void Solver<2>::IterateTridiagonal(int axis, 
	const std::vector<double>& enthalpy_data,
	const std::vector<double>& current_iteration_data,
	std::vector<double>& next_iteration_data,
	double dt)
{
	int idx1_p, idx1_c, idx1_m;
	int idx2_p, idx2_c, idx2_m;
	int num_1, num_2;
	int ix, iy;
	double h1, h2;

	if (axis == 0)
	{
		idx1_p = 1;
		idx1_c = 0;
		idx1_m = -1;
		idx2_p = mesh.GetNums()[0];
		idx2_c = 0;
		idx2_m = -mesh.GetNums()[0];
		num_1 = mesh.GetNums()[0]; num_2 = mesh.GetNums()[1];
		h1 = mesh.GetSteps()[0]; h2 = mesh.GetSteps()[1];
	}
	else if (axis == 1)
	{
		idx1_p = mesh.GetNums()[0];
		idx1_c = 0;
		idx1_m = -mesh.GetNums()[0];
		idx2_p = 1;
		idx2_c = 0;
		idx2_m = -1;
		num_1 = mesh.GetNums()[1]; num_2 = mesh.GetNums()[0];
		h1 =  mesh.GetSteps()[1]; h2 = mesh.GetSteps()[0];
	}

	MaterialInfo mi, min1, mip1, min2, mip2;
	
	for (int i2 = 1; i2 < num_2-1; ++i2)
	{
		/*
		a[0] = 0.0; b[0] = -1.0; c[0] = 1.0; f[0] = 0.0;
		a[num_1-1] = -1.0; b[num_1-1] = 1.0; c[num_1-1] = 0.0; f[num_1-1] = 0.0;
		*/

		// Setting boundary conditions
		SetBoundaryConditions(axis, {0, i2}, a, b, c, f);
		
		for (int i1 = 1; i1 < num_1-1; ++i1)
		{
			if (axis == 0)
			{
				ix = i1; iy = i2;
				mi = mesh.GetMaterialInfo({ix, iy});
				min1 = mesh.GetMaterialInfo({ix-1, iy});
				mip1 = mesh.GetMaterialInfo({ix+1, iy});
				min2 = mesh.GetMaterialInfo({ix, iy-1});
				mip2 = mesh.GetMaterialInfo({ix, iy+1});
			}
			else
			{
				ix = i2; iy = i1;
				mi = mesh.GetMaterialInfo({ix, iy});
				min1 = mesh.GetMaterialInfo({ix, iy-1});
				mip1 = mesh.GetMaterialInfo({ix, iy+1});
				min2 = mesh.GetMaterialInfo({ix-1, iy});
				mip2 = mesh.GetMaterialInfo({ix+1, iy});
			}
			
			const double * E_iter = current_iteration_data.data() + mesh.GetGlobalId({ix, iy});
			const double * E_step = enthalpy_data.data() + mesh.GetGlobalId({ix, iy});

			double kp = (mip1.GetThermalConductivityByE(E_iter[idx1_p]) + mi.GetThermalConductivityByE(E_iter[idx1_c])) / 2.0;
			double km = (mi.GetThermalConductivityByE(E_iter[idx1_c]) + min1.GetThermalConductivityByE(E_iter[idx1_m])) / 2.0;

			double coef = dt / (2.0 * h1 * h1);
			a[i1] = -coef * min1.GetAlpha(E_iter[idx1_m]) * km;
			b[i1] = 1.0 + coef * mi.GetAlpha(E_iter[idx1_c]) * (km + kp);
			c[i1] = -coef * mip1.GetAlpha(E_iter[idx1_p]) * kp;
			f[i1] = E_step[idx1_c] + coef * (min1.GetBeta(E_iter[idx1_m])*km - mi.GetBeta(E_iter[idx1_c])*(km+kp) + mip1.GetBeta(E_iter[idx1_p])*kp);

			double k2p = (mip2.GetThermalConductivityByE(E_step[idx2_p]) + mi.GetThermalConductivityByE(E_step[idx2_c])) / 2.0;
			double k2n = (min2.GetThermalConductivityByE(E_step[idx2_m]) + mi.GetThermalConductivityByE(E_step[idx2_c])) / 2.0;
			double T2p = mip2.GetTByEnthalpy(E_step[idx2_p]) - mi.GetTByEnthalpy(E_step[idx2_c]);
			double T2n = min2.GetTByEnthalpy(E_step[idx2_m]) - mi.GetTByEnthalpy(E_step[idx2_c]);
			f[i1] += k2p * T2p * dt/(2*h2*h2) + k2n * T2n * dt/(2*h2*h2);
		}

		SolveTridiagonal(num_1, a, b, c, f, x);
		for (int k = 0; k < num_1; ++k)
				next_iteration_data[idx1_p * k + idx2_p * i2] = x[k];
	}
}


template<>
void Solver<3>::IterateTridiagonal(int axis, 
	const std::vector<double>& enthalpy_data,
	const std::vector<double>& current_iteration_data,
	std::vector<double>& next_iteration_data,
	double dt)
{
	int idx1_p, idx1_m;
	int idx2_p, idx2_m;
	int idx3_p, idx3_m;
	int num_1, num_2, num_3;
	int ix, iy, iz;
	double h1, h2, h3;

	int idx_c = 0;
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
	}

	MaterialInfo mi, min1, mip1, min2, mip2, min3, mip3;
	
	for (int i3 = 1; i3 < num_3-1; ++i3)
	{
		for (int i2 = 1; i2 < num_2-1; ++i2)
		{
			/*
			a[0] = 0.0; b[0] = -1.0; c[0] = 1.0; f[0] = 0.0;
			a[num_1-1] = -1.0; b[num_1-1] = 1.0; c[num_1-1] = 0.0; f[num_1-1] = 0.0;
			*/

			// Setting boundary conditions
			SetBoundaryConditions(axis, {0, i2, i3}, a, b, c, f);


			for (int i1 = 1; i1 < num_1-1; ++i1)
			{
				if (axis == 0)
				{
					ix = i1; iy = i2; iz = i3;
					mi = mesh.GetMaterialInfo({ix, iy, iz});
					min1 = mesh.GetMaterialInfo({ix-1, iy, iz});
					mip1 = mesh.GetMaterialInfo({ix+1, iy, iz});
					min2 = mesh.GetMaterialInfo({ix, iy-1, iz});
					mip2 = mesh.GetMaterialInfo({ix, iy+1, iz});
					min3 = mesh.GetMaterialInfo({ix, iy, iz-1});
					mip3 = mesh.GetMaterialInfo({ix, iy, iz+1});
				}
				else if (axis == 1)
				{
					ix = i3; iy = i1; iz = i2;
					mi = mesh.GetMaterialInfo({ix, iy, iz});
					min1 = mesh.GetMaterialInfo({ix, iy-1, iz});
					mip1 = mesh.GetMaterialInfo({ix, iy+1, iz});
					min2 = mesh.GetMaterialInfo({ix, iy, iz-1});
					mip2 = mesh.GetMaterialInfo({ix, iy, iz+1});
					min3 = mesh.GetMaterialInfo({ix-1, iy, iz});
					mip3 = mesh.GetMaterialInfo({ix+1, iy, iz});
				}
				else
				{
					ix = i2; iy = i3; iz = i1;
					mi = mesh.GetMaterialInfo({ix, iy, iz});
					min1 = mesh.GetMaterialInfo({ix, iy, iz-1});
					mip1 = mesh.GetMaterialInfo({ix, iy, iz+1});
					min2 = mesh.GetMaterialInfo({ix-1, iy, iz});
					mip2 = mesh.GetMaterialInfo({ix+1, iy, iz});
					min3 = mesh.GetMaterialInfo({ix, iy-1, iz});
					mip3 = mesh.GetMaterialInfo({ix, iy+1, iz});
				}
				
				const double * E_iter = current_iteration_data.data() + mesh.GetGlobalId({ix, iy, iz});
				const double * E_step = enthalpy_data.data() + mesh.GetGlobalId({ix, iy, iz});

				double kp = (mip1.GetThermalConductivityByE(E_iter[idx1_p]) + mi.GetThermalConductivityByE(E_iter[idx_c])) / 2.0;
				double km = (mi.GetThermalConductivityByE(E_iter[idx_c]) + min1.GetThermalConductivityByE(E_iter[idx1_m])) / 2.0;

				double coef = dt / (3.0 * h1 * h1);
				a[i1] = -coef * min1.GetAlpha(E_iter[idx1_m]) * km;
				b[i1] = 1.0 + coef * mi.GetAlpha(E_iter[idx_c]) * (km + kp);
				c[i1] = -coef * mip1.GetAlpha(E_iter[idx1_p]) * kp;
				f[i1] = E_step[idx_c] + coef * (min1.GetBeta(E_iter[idx1_m])*km - mi.GetBeta(E_iter[idx_c])*(km+kp) + mip1.GetBeta(E_iter[idx1_p])*kp);

				double k2p = (mip2.GetThermalConductivityByE(E_step[idx2_p]) + mi.GetThermalConductivityByE(E_step[idx_c])) / 2.0;
				double k2n = (min2.GetThermalConductivityByE(E_step[idx2_m]) + mi.GetThermalConductivityByE(E_step[idx_c])) / 2.0;
				double T2p = mip2.GetTByEnthalpy(E_step[idx2_p]) - mi.GetTByEnthalpy(E_step[idx_c]);
				double T2n = min2.GetTByEnthalpy(E_step[idx2_m]) - mi.GetTByEnthalpy(E_step[idx_c]);
				f[i1] += k2p * T2p * dt/(3*h2*h2) + k2n * T2n * dt/(3*h2*h2);

				double k3p = (mip3.GetThermalConductivityByE(E_step[idx3_p]) + mi.GetThermalConductivityByE(E_step[idx_c])) / 2.0;
				double k3n = (min3.GetThermalConductivityByE(E_step[idx3_m]) + mi.GetThermalConductivityByE(E_step[idx_c])) / 2.0;
				double T3p = mip3.GetTByEnthalpy(E_step[idx3_p]) - mi.GetTByEnthalpy(E_step[idx_c]);
				double T3n = min3.GetTByEnthalpy(E_step[idx3_m]) - mi.GetTByEnthalpy(E_step[idx_c]);
				f[i1] += k3p * T3p * dt/(3*h3*h3) + k3n * T3n * dt/(3*h3*h3);
			}

			SolveTridiagonal(num_1, a, b, c, f, x);
			for (int k = 0; k < num_1; ++k)
					next_iteration_data[idx1_p * k + idx2_p * i2 + idx3_p * i3] = x[k];
		}
	}
}
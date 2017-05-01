#pragma once
#include <math.h>
#include <vector>
#include <functional>
#include <algorithm>
#include <iostream>
#include <fstream>

class mesh2d
{
	int num_x;
	int num_y;
	int number_of_nodes;
	double Lx;
	double Ly;
	double hx;
	double hy;


	std::vector<material_info> material_infos;
	std::function<int(double, double)> material_index;

public:

	mesh2d()
	{
	}
	mesh2d(int num_x, int num_y, double Lx, double Ly, std::vector<material_info> material_infos, std::function<int(double, double)> material_index) 
		: num_x(num_x), num_y(num_y), Lx(Lx), Ly(Ly), material_infos(material_infos), material_index(material_index)
	{
		hx = Lx / (num_x-1);
		hy = Ly / (num_y-1);
		number_of_nodes = num_x * num_y;
	}

	material_info & get_material(int n, int i)
	{
		return material_infos[material_index(n*hx, i*hy)];
	}

	inline double get_Lx() {return Lx;}
	inline double get_Ly() {return Ly;}
	inline double get_hx() {return hx;}
	inline double get_hy() {return hy;}
	inline int get_num_x() {return num_x;}
	inline int get_num_y() {return num_y;}
	inline int get_number_of_nodes() {return number_of_nodes;}
};


class solver2d
{
	mesh2d mesh;

	double * enthalpy_data;
	double * current_iteration_data;
	double * next_iteration_data;

	// Auxiliary arrays for tridiagonal method
	double *a, *b, *c, *f, *x;

	void set_enthalpy_by_T_data(const double * temperature_data)
	{
		int num_x = mesh.get_num_x();
		int num_y = mesh.get_num_y();
		for (int i = 0; i < num_y; ++i)
			for (int n = 0; n < num_x; ++n)
				enthalpy_data[n + i * num_x] = mesh.get_material(n, i).get_enthalpy_by_T(temperature_data[n + i * num_x]);
	}

	void iterate_tridiagonal(int axis, const double * enthalpy_data, const double * current_iteration_data, double * next_iteration_data, double dt)
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
			idx2_p = mesh.get_num_x();
			idx2_c = 0;
			idx2_m = -mesh.get_num_x();
			num_1 = mesh.get_num_x(); num_2 = mesh.get_num_y();
			h1 = mesh.get_hx(); h2 = mesh.get_hy();
		}
		else if (axis == 1)
		{
			idx1_p = mesh.get_num_x();
			idx1_c = 0;
			idx1_m = -mesh.get_num_x();
			idx2_p = 1;
			idx2_c = 0;
			idx2_m = -1;
			num_1 = mesh.get_num_y(); num_2 = mesh.get_num_x();
			h1 =  mesh.get_hy(); h2 = mesh.get_hx();
		}

		material_info mi, min1, mip1, min2, mip2;
		
		for (int i2 = 1; i2 < num_2-1; ++i2)
		{
			a[0] = 0.0; b[0] = -1.0; c[0] = 1.0; f[0] = 0.0;
			a[num_1-1] = -1.0; b[num_1-1] = 1.0; c[num_1-1] = 0.0; f[num_1-1] = 0.0;

			for (int i1 = 1; i1 < num_1-1; ++i1)
			{
				if (axis == 0)
				{
					ix = i1; iy = i2;
					mi = mesh.get_material(ix, iy);
					min1 = mesh.get_material(ix-1, iy);
					mip1 = mesh.get_material(ix+1, iy);
					min2 = mesh.get_material(ix, iy-1);
					mip2 = mesh.get_material(ix, iy+1);
				}
				else
				{
					ix = i2; iy = i1;
					mi = mesh.get_material(ix, iy);
					min1 = mesh.get_material(ix, iy-1);
					mip1 = mesh.get_material(ix, iy+1);
					min2 = mesh.get_material(ix-1, iy);
					mip2 = mesh.get_material(ix+1, iy);
				}
				
				const double * E_iter = current_iteration_data + ix + iy*mesh.get_num_x();
				const double * E_step = enthalpy_data + ix + iy*mesh.get_num_x();

				double kp = (mip1.get_thermal_conductivity_by_E(E_iter[idx1_p]) + mi.get_thermal_conductivity_by_E(E_iter[idx1_c])) / 2.0;
				double km = (mi.get_thermal_conductivity_by_E(E_iter[idx1_c]) + min1.get_thermal_conductivity_by_E(E_iter[idx1_m])) / 2.0;

				double coef = dt / (2.0 * h1 * h1);
				a[i1] = -coef * min1.get_alpha(E_iter[idx1_m]) * km;
				b[i1] = 1.0 + coef * mi.get_alpha(E_iter[idx1_c]) * (km + kp);
				c[i1] = -coef * mip1.get_alpha(E_iter[idx1_p]) * kp;
				f[i1] = E_step[idx1_c] + coef * (min1.get_beta(E_iter[idx1_m])*km - mi.get_beta(E_iter[idx1_c])*(km+kp) + mip1.get_beta(E_iter[idx1_p])*kp);

				double k2p = (mip2.get_thermal_conductivity_by_E(E_step[idx2_p]) + mi.get_thermal_conductivity_by_E(E_step[idx2_c])) / 2.0;
				double k2n = (min2.get_thermal_conductivity_by_E(E_step[idx2_m]) + mi.get_thermal_conductivity_by_E(E_step[idx2_c])) / 2.0;
				double T2p = mip2.get_T_by_enthalpy(E_step[idx2_p]) - mi.get_T_by_enthalpy(E_step[idx2_c]);
				double T2n = min2.get_T_by_enthalpy(E_step[idx2_m]) - mi.get_T_by_enthalpy(E_step[idx2_c]);
				f[i1] += k2p * T2p * dt/(2*h2*h2) + k2n * T2n * dt/(2*h2*h2);
			}

			solve_tridiagonal(num_1, a, b, c, f, x);
			for (int k = 0; k < num_1; ++k)
					next_iteration_data[idx1_p * k + idx2_p * i2] = x[k];
		}
	}

	void solve_tridiagonal (int n, double *a, double *b, double *c, double *f, double *x)
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

	double calculate_L2_norm(const double * A, const double * B)
	{
		int num_x = mesh.get_num_x();
		int num_y = mesh.get_num_y();
		double L2 = 0.0;
		for (int n = 0; n < num_x; ++n)
			for (int i = 0; i < num_y; ++i)
			{
				L2 += (A[n + i*num_x] - B[n + i*num_x]) * (A[n + i*num_x] - B[n + i*num_x]);
			}
		L2 = sqrt(L2);
		return L2;
	}

	double calculate_L2_norm(const double * A)
	{
		int num_x = mesh.get_num_x();
		int num_y = mesh.get_num_y();
		double L2 = 0.0;
		for (int n = 0; n < num_x; ++n)
			for (int i = 0; i < num_y; ++i)
			{
				L2 += A[n + i*num_x] * A[n + i*num_x];
			}
		L2 = sqrt(L2);
		return L2;
	}

	double calculate_Linf_norm(const double * A, const double * B)
	{
		int num_x = mesh.get_num_x();
		int num_y = mesh.get_num_y();
		double Linf = 0.0;
		for (int n = 0; n < num_x; ++n)
			for (int i = 0; i < num_y; ++i)
			{
				double current_norm = fabs(A[n + i*num_x] - B[n + i*num_x]);
				if (current_norm > Linf)
					Linf = current_norm;
			}
		return Linf;
	}

	double calculate_Linf_norm(const double * A)
	{
		int num_x = mesh.get_num_x();
		int num_y = mesh.get_num_y();
		double Linf = 0.0;
		for (int n = 0; n < num_x; ++n)
			for (int i = 0; i < num_y; ++i)
			{
				double current_norm = fabs(A[n + i*num_x]);
				if (current_norm > Linf)
					Linf = current_norm;
			}
		return Linf;
	}

	double assign(double * A, const double * B)
	{
		int num_x = mesh.get_num_x();
		int num_y = mesh.get_num_y();
		for (int n = 0; n < num_x; ++n)
			for (int i = 0; i < num_y; ++i)
			{
				A[n + i*num_x] = B[n + i*num_x];
			}
	}

	void step(int axis, double dt)
	{
		double L2_norm = 1.0;
		assign(next_iteration_data, enthalpy_data);
		while (L2_norm > 1e-3)
		{
			assign(current_iteration_data, next_iteration_data);
			iterate_tridiagonal(axis, enthalpy_data, current_iteration_data, next_iteration_data, dt);
			L2_norm = calculate_Linf_norm(current_iteration_data, next_iteration_data) / calculate_Linf_norm(current_iteration_data);
		}
		assign(enthalpy_data, next_iteration_data);
	}



public:
	solver2d(mesh2d & mesh, const double * temperature_data) 
	{
		this->mesh = mesh;
		enthalpy_data = new double[mesh.get_number_of_nodes()];
		current_iteration_data = new double [mesh.get_number_of_nodes()];
		next_iteration_data = new double [mesh.get_number_of_nodes()];
		set_enthalpy_by_T_data(temperature_data);

		int max_n = std::max(mesh.get_num_x(), mesh.get_num_y());
		a = new double[max_n];
		b = new double[max_n];
		c = new double[max_n];
		f = new double[max_n];
		x = new double[max_n];
	}

	~solver2d()
	{
		delete [] enthalpy_data;
		delete [] current_iteration_data;
		delete [] next_iteration_data;

		delete [] a;
		delete [] b;
		delete [] c;
		delete [] f;
		delete [] x;
	}

	void step(double dt)
	{
		step(0, dt);
		step(1, dt);
	}

	void save_to_vtk(std::string name)
	{
		int num_x = mesh.get_num_x();
		int num_y = mesh.get_num_y();

		std::ofstream vtk_file(name.c_str());
		vtk_file << "# vtk DataFile Version 3.0\nVx data\nASCII\n\n";
		vtk_file << "DATASET STRUCTURED_POINTS\nDIMENSIONS " << num_x << " " << num_y << " 1" << "\n";
		vtk_file << "ORIGIN 0 0 0\n" << "SPACING " << mesh.get_hx() << " " << mesh.get_hy() << " 1" << "\n";

		vtk_file << "\nPOINT_DATA " << num_x * num_y << "\n";
		vtk_file << "SCALARS E FLOAT\n";
		vtk_file << "LOOKUP_TABLE default\n";
		
		for (int j = 0; j < num_y; j++)
			for (int i = 0; i < num_x; i++)
			{
				vtk_file <<  enthalpy_data[i + j*num_x] << "\n";
			}
		vtk_file << "SCALARS T FLOAT\n";
		vtk_file << "LOOKUP_TABLE default\n";
		for (int j = 0; j < num_y; j++)
			for (int i = 0; i < num_x; i++)
			{
				material_info & mi = mesh.get_material(i, j);
				vtk_file <<  mi.get_T_by_enthalpy(enthalpy_data[i + j*num_x]) << "\n";
			}
		vtk_file << "SCALARS K FLOAT\n";
		vtk_file << "LOOKUP_TABLE default\n";
		for (int j = 0; j < num_y; j++)
			for (int i = 0; i < num_x; i++)
			{
				material_info & mi = mesh.get_material(i, j);
				vtk_file <<  mi.get_thermal_conductivity_by_E(enthalpy_data[i + j*num_x]) << "\n";
			}
		vtk_file << "SCALARS state FLOAT\n";
		vtk_file << "LOOKUP_TABLE default\n";
		for (int j = 0; j < num_y; j++)
			for (int i = 0; i < num_x; i++)
			{
				material_info & mi = mesh.get_material(i, j);
				double T = mi.get_T_by_enthalpy(enthalpy_data[i + j*num_x]);
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

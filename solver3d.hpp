#pragma once
#include <math.h>
#include <vector>
#include <functional>
#include <algorithm>
#include <iostream>
#include <fstream>

class mesh3d
{
	int num_x;
	int num_y;
	int num_z;
	int number_of_nodes;
	double Lx;
	double Ly;
	double Lz;
	double hx;
	double hy;
	double hz;


	std::vector<material_info> material_infos;
	std::function<int(double, double, double)> material_index;

public:

	mesh3d()
	{
	}
	mesh3d(int num_x, int num_y, int num_z, double Lx, double Ly, double Lz, std::vector<material_info> material_infos, 
		std::function<int(double, double, double)> material_index) 
		: num_x(num_x), num_y(num_y), num_z(num_z), Lx(Lx), Ly(Ly), Lz(Lz), material_infos(material_infos), material_index(material_index)
	{
		hx = Lx / (num_x-1);
		hy = Ly / (num_y-1);
		hz = Lz / (num_z-1);
		number_of_nodes = num_x * num_y * num_z;
	}

	material_info & get_material(int n, int i, int l)
	{
		return material_infos[material_index(n*hx, i*hy, l*hz)];
	}

	inline double get_Lx() {return Lx;}
	inline double get_Ly() {return Ly;}
	inline double get_Lz() {return Lz;}
	inline double get_hx() {return hx;}
	inline double get_hy() {return hy;}
	inline double get_hz() {return hz;}
	inline int get_num_x() {return num_x;}
	inline int get_num_y() {return num_y;}
	inline int get_num_z() {return num_z;}
	inline int get_number_of_nodes() {return number_of_nodes;}
};


class solver3d
{
	mesh3d mesh;

	double * enthalpy_data;
	double * current_iteration_data;
	double * next_iteration_data;

	// Auxiliary arrays for tridiagonal method
	double *a, *b, *c, *f, *x;

	inline int id(int n, int i, int l) {return n + i*mesh.get_num_x() + l * mesh.get_num_x() * mesh.get_num_y();}

	void set_enthalpy_by_T_data(const double * temperature_data)
	{
		int num_x = mesh.get_num_x();
		int num_y = mesh.get_num_y();
		int num_z = mesh.get_num_z();
		for (int i = 0; i < num_y; ++i)
			for (int n = 0; n < num_x; ++n)
				for (int l = 0; l < num_z; ++l)
					enthalpy_data[id(n,i,l)] = mesh.get_material(n, i, l).get_enthalpy_by_T(temperature_data[id(n,i,l)]);
	}

	void iterate_tridiagonal(int axis, const double * enthalpy_data, const double * current_iteration_data, double * next_iteration_data, double dt)
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
			idx2_p = mesh.get_num_x();
			idx2_m = -mesh.get_num_x();
			idx3_p = mesh.get_num_x()*mesh.get_num_y();
			idx3_m = -mesh.get_num_x()*mesh.get_num_y();
			num_1 = mesh.get_num_x(); num_2 = mesh.get_num_y(); num_3 = mesh.get_num_z();
			h1 = mesh.get_hx(); h2 = mesh.get_hy(); h3 = mesh.get_hz();
		}
		else if (axis == 1)
		{
			idx1_p = mesh.get_num_x();
			idx1_m = -mesh.get_num_x();
			idx2_p = mesh.get_num_x()*mesh.get_num_y();
			idx2_m = -mesh.get_num_x()*mesh.get_num_y();
			idx3_p = 1;
			idx3_m = -1;
			num_1 = mesh.get_num_y(); num_2 = mesh.get_num_z(); num_3 = mesh.get_num_x();
			h1 =  mesh.get_hy(); h2 = mesh.get_hz();  h3 = mesh.get_hx();
		}
		else if (axis == 2)
		{
			idx1_p = mesh.get_num_x()*mesh.get_num_y();
			idx1_m = -mesh.get_num_x()*mesh.get_num_y();
			idx2_p = 1;
			idx2_m = -1;
			idx3_p = mesh.get_num_x();
			idx3_m = -mesh.get_num_x();
			num_1 = mesh.get_num_z(); num_2 = mesh.get_num_x(); num_3 = mesh.get_num_y();
			h1 =  mesh.get_hz(); h2 = mesh.get_hx();  h3 = mesh.get_hy();
		}

		material_info mi, min1, mip1, min2, mip2, min3, mip3;
		
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
						mi = mesh.get_material(ix, iy, iz);
						min1 = mesh.get_material(ix-1, iy, iz);
						mip1 = mesh.get_material(ix+1, iy, iz);
						min2 = mesh.get_material(ix, iy-1, iz);
						mip2 = mesh.get_material(ix, iy+1, iz);
						min3 = mesh.get_material(ix, iy, iz-1);
						mip3 = mesh.get_material(ix, iy, iz+1);
					}
					else if (axis == 1)
					{
						ix = i3; iy = i1; iz = i2;
						mi = mesh.get_material(ix, iy, iz);
						min1 = mesh.get_material(ix, iy-1, iz);
						mip1 = mesh.get_material(ix, iy+1, iz);
						min2 = mesh.get_material(ix, iy, iz-1);
						mip2 = mesh.get_material(ix, iy, iz+1);
						min3 = mesh.get_material(ix-1, iy, iz);
						mip3 = mesh.get_material(ix+1, iy, iz);
					}
					else
					{
						ix = i2; iy = i3; iz = i1;
						mi = mesh.get_material(ix, iy, iz);
						min1 = mesh.get_material(ix, iy, iz-1);
						mip1 = mesh.get_material(ix, iy, iz+1);
						min2 = mesh.get_material(ix-1, iy, iz);
						mip2 = mesh.get_material(ix+1, iy, iz);
						min3 = mesh.get_material(ix, iy-1, iz);
						mip3 = mesh.get_material(ix, iy+1, iz);
					}
					
					const double * E_iter = current_iteration_data + id(ix, iy, iz);
					const double * E_step = enthalpy_data + id(ix, iy, iz);

					double kp = (mip1.get_thermal_conductivity_by_E(E_iter[idx1_p]) + mi.get_thermal_conductivity_by_E(E_iter[idx_c])) / 2.0;
					double km = (mi.get_thermal_conductivity_by_E(E_iter[idx_c]) + min1.get_thermal_conductivity_by_E(E_iter[idx1_m])) / 2.0;

					double coef = dt / (3.0 * h1 * h1);
					a[i1] = -coef * min1.get_alpha(E_iter[idx1_m]) * km;
					b[i1] = 1.0 + coef * mi.get_alpha(E_iter[idx_c]) * (km + kp);
					c[i1] = -coef * mip1.get_alpha(E_iter[idx1_p]) * kp;
					f[i1] = E_step[idx_c] + coef * (min1.get_beta(E_iter[idx1_m])*km - mi.get_beta(E_iter[idx_c])*(km+kp) + mip1.get_beta(E_iter[idx1_p])*kp);

					double k2p = (mip2.get_thermal_conductivity_by_E(E_step[idx2_p]) + mi.get_thermal_conductivity_by_E(E_step[idx_c])) / 2.0;
					double k2n = (min2.get_thermal_conductivity_by_E(E_step[idx2_m]) + mi.get_thermal_conductivity_by_E(E_step[idx_c])) / 2.0;
					double T2p = mip2.get_T_by_enthalpy(E_step[idx2_p]) - mi.get_T_by_enthalpy(E_step[idx_c]);
					double T2n = min2.get_T_by_enthalpy(E_step[idx2_m]) - mi.get_T_by_enthalpy(E_step[idx_c]);
					f[i1] += k2p * T2p * dt/(3*h2*h2) + k2n * T2n * dt/(3*h2*h2);

					double k3p = (mip3.get_thermal_conductivity_by_E(E_step[idx3_p]) + mi.get_thermal_conductivity_by_E(E_step[idx_c])) / 2.0;
					double k3n = (min3.get_thermal_conductivity_by_E(E_step[idx3_m]) + mi.get_thermal_conductivity_by_E(E_step[idx_c])) / 2.0;
					double T3p = mip3.get_T_by_enthalpy(E_step[idx3_p]) - mi.get_T_by_enthalpy(E_step[idx_c]);
					double T3n = min3.get_T_by_enthalpy(E_step[idx3_m]) - mi.get_T_by_enthalpy(E_step[idx_c]);
					f[i1] += k3p * T3p * dt/(3*h3*h3) + k3n * T3n * dt/(3*h3*h3);
				}

				solve_tridiagonal(num_1, a, b, c, f, x);
				for (int k = 0; k < num_1; ++k)
						next_iteration_data[idx1_p * k + idx2_p * i2 + idx3_p * i3] = x[k];
			}
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

	double calculate_Linf_norm(const double * A, const double * B)
	{
		double Linf = 0.0;
		for (int n = 0; n < mesh.get_num_x(); ++n)
			for (int i = 0; i < mesh.get_num_y(); ++i)
				for (int l = 0; l < mesh.get_num_z(); ++l)
				{
					double current_norm = fabs(A[id(n, i, l)] - B[id(n, i, l)]);
					if (current_norm > Linf)
						Linf = current_norm;
				}
		return Linf;
	}

	double calculate_Linf_norm(const double * A)
	{
		double Linf = 0.0;
		for (int n = 0; n < mesh.get_num_x(); ++n)
			for (int i = 0; i < mesh.get_num_y(); ++i)
				for (int l = 0; l < mesh.get_num_z(); ++l)
				{
					double current_norm = fabs(A[id(n, i, l)]);
					if (current_norm > Linf)
						Linf = current_norm;
				}
		return Linf;
	}

	double assign(double * A, const double * B)
	{
		for (int n = 0; n < mesh.get_num_x(); ++n)
			for (int i = 0; i < mesh.get_num_y(); ++i)
				for (int l = 0; l < mesh.get_num_z(); ++l)
				{
					A[id(n, i, l)] = B[id(n, i, l)];
				}
	}

	void step(int axis, double dt)
	{
		double norm = 1.0;
		assign(next_iteration_data, enthalpy_data);
		while (norm > 1e-3)
		{
			assign(current_iteration_data, next_iteration_data);
			iterate_tridiagonal(axis, enthalpy_data, current_iteration_data, next_iteration_data, dt);
			norm = calculate_Linf_norm(current_iteration_data, next_iteration_data) / calculate_Linf_norm(current_iteration_data);
		}
		assign(enthalpy_data, next_iteration_data);
	}



public:
	solver3d(mesh3d & mesh, const double * temperature_data) 
	{
		this->mesh = mesh;
		enthalpy_data = new double[mesh.get_number_of_nodes()];
		current_iteration_data = new double [mesh.get_number_of_nodes()];
		next_iteration_data = new double [mesh.get_number_of_nodes()];
		set_enthalpy_by_T_data(temperature_data);

		int max_n = std::max(std::max(mesh.get_num_x(), mesh.get_num_y()), mesh.get_num_z());
		a = new double[max_n];
		b = new double[max_n];
		c = new double[max_n];
		f = new double[max_n];
		x = new double[max_n];
	}

	~solver3d()
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
		step(2, dt);
	}

	void save_to_vtk(std::string name)
	{
		int num_x = mesh.get_num_x();
		int num_y = mesh.get_num_y();
		int num_z = mesh.get_num_z();

		std::ofstream vtk_file(name.c_str());
		vtk_file << "# vtk DataFile Version 3.0\nVx data\nASCII\n\n";
		vtk_file << "DATASET STRUCTURED_POINTS\nDIMENSIONS " << num_x << " " << num_y << " " << num_z << "\n";
		vtk_file << "ORIGIN 0 0 0\n" << "SPACING " << mesh.get_hx() << " " << mesh.get_hy() << " " << mesh.get_hz() << "\n";

		vtk_file << "\nPOINT_DATA " << num_x * num_y * num_z << "\n";
		vtk_file << "SCALARS E FLOAT\n";
		vtk_file << "LOOKUP_TABLE default\n";

		for (int l = 0; l < num_z; l++)
			for (int i = 0; i < num_y; i++)
				for (int n = 0; n < num_x; n++)
				{
					vtk_file <<  enthalpy_data[id(n, i, l)] << "\n";
				}
		vtk_file << "SCALARS T FLOAT\n";
		vtk_file << "LOOKUP_TABLE default\n";
		for (int l = 0; l < num_z; l++)
			for (int i = 0; i < num_y; i++)
				for (int n = 0; n < num_x; n++)
				{
					material_info & mi = mesh.get_material(n, i, l);
					vtk_file <<  mi.get_T_by_enthalpy(enthalpy_data[id(n, i, l)]) << "\n";
				}
		vtk_file << "SCALARS K FLOAT\n";
		vtk_file << "LOOKUP_TABLE default\n";
		for (int l = 0; l < num_z; l++)
			for (int i = 0; i < num_y; i++)
				for (int n = 0; n < num_x; n++)
				{
					material_info & mi = mesh.get_material(n, i, l);
					vtk_file <<  mi.get_thermal_conductivity_by_E(enthalpy_data[id(n, i, l)]) << "\n";
				}
		vtk_file << "SCALARS state FLOAT\n";
		vtk_file << "LOOKUP_TABLE default\n";
		for (int l = 0; l < num_z; l++)
			for (int i = 0; i < num_y; i++)
				for (int n = 0; n < num_x; n++)
				{
					material_info & mi = mesh.get_material(n, i, l);
					double T = mi.get_T_by_enthalpy(enthalpy_data[id(n, i, l)]);
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

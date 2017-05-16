#pragma once
#include <math.h>
#include <vector>
#include <functional>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "../math/Spaces.h"


template <typename Space>
class mesh2d
{
	SPACE_TYPEDEFS;

	IndexVector stride;
	AABB aabb;

	Scalar hx;
	Scalar hy;

	using MaterialParamsGetter = std::function<material_info<Scalar>(const Vector&)>;
	MaterialParamsGetter materialParamsGetter;

public:

	mesh2d() = default;

	mesh2d(const IndexVector& stride, const AABB& aabb,
		MaterialParamsGetter materialParamsGetter)
		: stride(stride), aabb(aabb), materialParamsGetter(materialParamsGetter)
	{
		auto size = aabb.Size();
		hx = size.x / (stride.x - 1);
		hy = size.y / (stride.y - 1);
	}

	material_info<Scalar> get_material(int n, int i) const
	{
		return materialParamsGetter(aabb.boxPoint1 + Vector(n * hx, i * hy));
	}


	IndexVector getStride() const { return stride; }
	inline Scalar get_hx() const {return hx;}
	inline Scalar get_hy() const {return hy;}
	inline int get_num_x() const { return stride.x; }
	inline int get_num_y() const { return stride.y; }
	inline IndexType get_number_of_nodes() const {return stride.GetVolume();}
};

template <typename Space>
class solver2d
{
	SPACE_TYPEDEFS;

	mesh2d<Space> mesh;

	Scalar* enthalpy_data;
	Scalar* current_iteration_data;
	Scalar* next_iteration_data;

	// Auxiliary arrays for tridiagonal method
	Scalar *a, *b, *c, *f, *x;


	inline IndexType id(IndexType n, IndexType i) { return idx<Space>({ n, i }, mesh.getStride()); }

	void set_enthalpy_by_T_data(const Scalar * temperature_data)
	{
		for (int i = 0; i < mesh.get_num_y(); ++i)
			for (int n = 0; n < mesh.get_num_x(); ++n)
			{
				const auto index = id(n, i);
				enthalpy_data[index] = mesh.get_material(n, i).get_enthalpy_by_T(temperature_data[index]);
			}
	}

	void iterate_tridiagonal(int axis, 
		const Scalar * enthalpy_data, 
		const Scalar * current_iteration_data, 
		Scalar * next_iteration_data, 
		Scalar dt)
	{
		int idx1_p, idx1_c, idx1_m;
		int idx2_p, idx2_c, idx2_m;
		int num_1, num_2;
		int ix, iy;
		Scalar h1, h2;

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

		material_info<Scalar> mi, min1, mip1, min2, mip2;
		
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
				
				const Scalar * E_iter = current_iteration_data + id(ix, iy);
				const Scalar * E_step = enthalpy_data + id(ix, iy);

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

	void solve_tridiagonal (int n, Scalar *a, Scalar *b, Scalar *c, Scalar *f, Scalar *x)
	{
		Scalar m;
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

	Scalar calculate_L2_norm(const Scalar * A, const Scalar * B)
	{
		int num_x = mesh.get_num_x();
		int num_y = mesh.get_num_y();
		Scalar L2 = 0.0;
		for (int n = 0; n < num_x; ++n)
			for (int i = 0; i < num_y; ++i)
			{
				L2 += Sqr(A[id(n, i)] - B[id(n, i)]);
			}
		L2 = sqrt(L2);
		return L2;
	}

	Scalar calculate_L2_norm(const Scalar * A)
	{
		int num_x = mesh.get_num_x();
		int num_y = mesh.get_num_y();
		Scalar L2 = 0.0;
		for (int n = 0; n < num_x; ++n)
			for (int i = 0; i < num_y; ++i)
			{
				L2 += Sqr(A[id(n, i)]);
			}
		L2 = sqrt(L2);
		return L2;
	}

	Scalar calculate_Linf_norm(const Scalar * A, const Scalar * B)
	{
		int num_x = mesh.get_num_x();
		int num_y = mesh.get_num_y();
		Scalar Linf = 0.0;
		for (int n = 0; n < num_x; ++n)
			for (int i = 0; i < num_y; ++i)
			{
				Scalar current_norm = fabs(A[id(n, i)] - B[id(n, i)]);
				if (current_norm > Linf)
					Linf = current_norm;
			}
		return Linf;
	}

	Scalar calculate_Linf_norm(const Scalar * A)
	{
		int num_x = mesh.get_num_x();
		int num_y = mesh.get_num_y();
		Scalar Linf = 0.0;
		for (int n = 0; n < num_x; ++n)
			for (int i = 0; i < num_y; ++i)
			{
				Scalar current_norm = fabs(A[id(n, i)]);
				if (current_norm > Linf)
					Linf = current_norm;
			}
		return Linf;
	}

	void assign(Scalar * A, const Scalar * B)
	{
		int num_x = mesh.get_num_x();
		int num_y = mesh.get_num_y();
		for (int n = 0; n < num_x; ++n)
			for (int i = 0; i < num_y; ++i)
			{
				A[n + i*num_x] = B[n + i*num_x];
			}
	}

	void step(int axis, Scalar dt)
	{
		Scalar L2_norm = 1.0;
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
	solver2d(mesh2d<Space>& mesh, const Scalar * temperature_data) : mesh(mesh)
	{
		enthalpy_data = new Scalar[mesh.get_number_of_nodes()];
		current_iteration_data = new Scalar [mesh.get_number_of_nodes()];
		next_iteration_data = new Scalar [mesh.get_number_of_nodes()];
		set_enthalpy_by_T_data(temperature_data);

		int max_n = std::max(mesh.get_num_x(), mesh.get_num_y());
		a = new Scalar[max_n];
		b = new Scalar[max_n];
		c = new Scalar[max_n];
		f = new Scalar[max_n];
		x = new Scalar[max_n];
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

	void step(Scalar dt)
	{
		step(0, dt);
		step(1, dt);
	}

	void save_to_vtk(const std::string& name)
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
				auto mi = mesh.get_material(i, j);
				vtk_file <<  mi.get_T_by_enthalpy(enthalpy_data[i + j*num_x]) << "\n";
			}
		vtk_file << "SCALARS K FLOAT\n";
		vtk_file << "LOOKUP_TABLE default\n";
		for (int j = 0; j < num_y; j++)
			for (int i = 0; i < num_x; i++)
			{
				auto mi = mesh.get_material(i, j);
				vtk_file <<  mi.get_thermal_conductivity_by_E(enthalpy_data[i + j*num_x]) << "\n";
			}
		vtk_file << "SCALARS state FLOAT\n";
		vtk_file << "LOOKUP_TABLE default\n";
		for (int j = 0; j < num_y; j++)
			for (int i = 0; i < num_x; i++)
			{
				auto mi = mesh.get_material(i, j);
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

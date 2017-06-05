#pragma once
#include "../math/Spaces.h"
#include "mesh.h"

template <typename Space>
class SolverBase
{
	SPACE_TYPEDEFS;
public:
	SolverBase(const Mesh<Space>& mesh, const Scalar* temperature_data) : mesh(mesh)
	{
		enthalpy_data.resize(mesh.get_number_of_nodes());
		current_iteration_data.resize(mesh.get_number_of_nodes());
		next_iteration_data.resize(mesh.get_number_of_nodes());

		set_enthalpy_by_T_data(temperature_data);

		const auto max_n = mesh.getStride().ComponentMax();
		a.resize(max_n);
		b.resize(max_n);
		c.resize(max_n);
		f.resize(max_n);
		x.resize(max_n);
	}

	void step(Scalar dt)
	{
		for (int stepIndex = 0; stepIndex < Space::Dimension; ++stepIndex)
		{
			step(stepIndex, dt);
		}
	}

protected:

	void set_enthalpy_by_T_data(const Scalar* temperature_data)
	{
		for (const auto& node : AllNodes<Space>(mesh.getStride()))
		{
			const auto index = id(node);
			enthalpy_data[index] = mesh.get_material(node).get_enthalpy_by_T(temperature_data[index]);
		}
	}

	IndexType id(const IndexVector& coord) { return idx<Space>(coord, mesh.getStride()); }

	static void solve_tridiagonal(int n, 
		const std::vector<Scalar>& a, 
		std::vector<Scalar>& b, 
		const std::vector<Scalar>& c, 
		std::vector<Scalar>& f, 
		std::vector<Scalar>& x)
	{
		Scalar m;
		for (int i = 1; i < n; i++)
		{
			m = a[i] / b[i - 1];
			b[i] = b[i] - m * c[i - 1];
			f[i] = f[i] - m * f[i - 1];
		}

		x[n - 1] = f[n - 1] / b[n - 1];

		for (int i = n - 2; i >= 0; i--)
			x[i] = (f[i] - c[i] * x[i + 1]) / b[i];
	}

	static Scalar calculate_Linf_norm(const std::vector<Scalar>& a, const std::vector<Scalar>& b)
	{
		assert(a.size() == b.size());
		Scalar Linf = 0.0;
		for (IndexType i = 0; i < a.size(); ++i)
		{
			Scalar current_norm = fabs(a[i] - b[i]);
			if (current_norm > Linf)
				Linf = current_norm;
		}
		return Linf;
	}

	static Scalar calculate_Linf_norm(const std::vector<Scalar>& a)
	{
		Scalar Linf = 0.0;
		for (int i = 0; i < a.size(); ++i)
		{
			Scalar current_norm = fabs(a[i]);
			if (current_norm > Linf)
				Linf = current_norm;
		}
		return Linf;
	}

	static void assign(std::vector<Scalar>& to, const std::vector<Scalar>& from)
	{
		assert(to.size() == from.size());
		std::copy_n(from.begin(), to.size(), to.begin());
	}

	static Scalar calculate_L2_norm(const std::vector<Scalar>& a, 
		const std::vector<Scalar>& b)
	{
		assert(a.size() == b.size());
		Scalar L2 = 0.0;
		for (int i = 0; i < a.size(); ++i)
		{
			L2 += Sqr(a[i] - b[i]);
		}
		L2 = sqrt(L2);
		return L2;
	}

	static Scalar calculate_L2_norm(const std::vector<Scalar>& a)
	{
		Scalar L2 = 0.0;
		for (int i = 0; i < a.size(); ++i)
		{
			L2 += Sqr(a[i]);
		}
		L2 = sqrt(L2);
		return L2;
	}

	virtual void iterate_tridiagonal(int axis,
		const std::vector<Scalar>& enthalpy_data,
		const std::vector<Scalar>& current_iteration_data,
		std::vector<Scalar>& next_iteration_data,
		Scalar dt) = 0;

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

	Mesh<Space> mesh;

	std::vector<Scalar> enthalpy_data;
	std::vector<Scalar> current_iteration_data;
	std::vector<Scalar> next_iteration_data;

	// Auxiliary arrays for tridiagonal method
	std::vector<Scalar> a, b, c, f, x;
};
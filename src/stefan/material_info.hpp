#pragma once

#include <assert.h>

template <typename Scalar>
struct material_info
{
	// material melts from T1 to T2; T1~T2
	Scalar T1 = 0;
	Scalar T2 = 0;

	Scalar rho_L = 0;
	Scalar rho_S = 0;
	Scalar thermal_conductivity_L = 0;
	Scalar thermal_conductivity_S = 0;
	Scalar specific_heat_fusion = 0;
	Scalar specific_heat_capacity_L = 0;
	Scalar specific_heat_capacity_S = 0;

	Scalar get_enthalpy_by_T(Scalar T) const 
	{
		Scalar E;
		if (T < T1)
			E = get_b() * T;
		else if (T < T2)
			E = (T - get_p2()) / get_p1();
		else
			E = get_r1() * T + get_r2();
		return E;
	}

	Scalar get_T_by_enthalpy(Scalar enthalpy) const 
	{
		Scalar E1 = get_b() * T1;
		Scalar E2 = get_r1() * T2 + get_r2();
		assert(E1 < E2);

		Scalar t = 0;

		if (enthalpy < E1)
			t = enthalpy / get_b();
		else if (enthalpy < E2)
			t = enthalpy * get_p1() + get_p2();
		else 
			t = (enthalpy - get_r2()) / get_r1();
		return t;
	}

	Scalar get_thermal_conductivity_by_E(Scalar enthalpy) const 
	{
		Scalar E1 = get_b() * T1;
		Scalar E2 = get_r1() * T2 + get_r2();
		if (enthalpy < E1)
			return thermal_conductivity_S;
		else if (enthalpy < E2)
			return ((thermal_conductivity_L - thermal_conductivity_S) * enthalpy + thermal_conductivity_S * E2 - thermal_conductivity_L * E1) / (E2 - E1);
		else
			return thermal_conductivity_L;
	}

	Scalar get_alpha(Scalar enthalpy) const 
	{
		Scalar E1 = get_b() * T1;
		Scalar E2 = get_r1() * T2 + get_r2();
		if (enthalpy < E1)
			return 1.0 / get_b();
		else if (enthalpy < E2)
			return get_p1();
		else
			return 1.0 / get_r1();
	}


	Scalar get_beta(Scalar enthalpy) const 
	{
		Scalar E1 = get_b() * T1;
		Scalar E2 = get_r1() * T2 + get_r2();
		if (enthalpy < E1)
			return Scalar(0.0);
		else if (enthalpy < E2)
			return get_p2();
		else
			return -get_r2() / get_r1();
	}

private:
	Scalar get_b() const 
	{
		return rho_S * specific_heat_capacity_S;
	}

	Scalar get_p1() const 
	{
		return (T2 - T1) / (rho_S * specific_heat_fusion);
	}

	Scalar get_p2() const 
	{
		return T1 * (1.0 - specific_heat_capacity_S * (T2 - T1) / specific_heat_fusion);
	}

	Scalar get_r1() const 
	{
		return rho_L * specific_heat_capacity_L;
	}

	Scalar get_r2() const
	{
		return rho_S * specific_heat_fusion + rho_S * specific_heat_capacity_S * T1 - rho_L * specific_heat_capacity_L * T2;
	}
};
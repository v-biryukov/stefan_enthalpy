#pragma once


struct material_info
{
	// material melts from T1 to T2; T1~T2
	double T1;
	double T2;

	double rho_L;
	double rho_S;
	double thermal_conductivity_L;
	double thermal_conductivity_S;
	double specific_heat_fusion;
	double specific_heat_capacity_L;
	double specific_heat_capacity_S;


	material_info()
	{
	}
	material_info(double T1, double T2, double rho_L, double rho_S, double thermal_conductivity_L, double thermal_conductivity_S,
		double specific_heat_fusion, double specific_heat_capacity_L, double specific_heat_capacity_S)
		: T1(T1), T2(T2), rho_L(rho_L), rho_S(rho_S), thermal_conductivity_L(thermal_conductivity_L), thermal_conductivity_S(thermal_conductivity_S),
		specific_heat_fusion(specific_heat_fusion), specific_heat_capacity_L(specific_heat_capacity_L), specific_heat_capacity_S(specific_heat_capacity_S)
	{
	}

	double get_enthalpy_by_T(double T)
	{
		double E;
		if (T < T1)
			E = get_b() * T;
		else if (T < T2)
			E = (T - get_p2()) / get_p1();
		else
			E = get_r1() * T + get_r2();
		return E;
	}

	double get_T_by_enthalpy(double enthalpy)
	{
		double E1 = get_b() * T1;
		double E2 = get_r1() * T2 + get_r2();
		if (enthalpy < E1)
			return enthalpy / get_b();
		else if (enthalpy < E2)
			return enthalpy * get_p1() + get_p2();
		else 
			return (enthalpy - get_r2()) / get_r1();
	}

	double get_thermal_conductivity_by_E(double enthalpy)
	{
		double E1 = get_b() * T1;
		double E2 = get_r1() * T2 + get_r2();
		if (enthalpy < E1)
			return thermal_conductivity_S;
		else if (enthalpy < E2)
			return ((thermal_conductivity_L - thermal_conductivity_S) * enthalpy + thermal_conductivity_S * E2 - thermal_conductivity_L * E1) / (E2 - E1);
		else
			return thermal_conductivity_L;
	}

	double get_alpha(double enthalpy)
	{
		double E1 = get_b() * T1;
		double E2 = get_r1() * T2 + get_r2();
		if (enthalpy < E1)
			return 1.0 / get_b();
		else if (enthalpy < E2)
			return get_p1();
		else
			return 1.0 / get_r1();
	}


	double get_beta(double enthalpy)
	{
		double E1 = get_b() * T1;
		double E2 = get_r1() * T2 + get_r2();
		if (enthalpy < E1)
			return 0.0;
		else if (enthalpy < E2)
			return get_p2();
		else
			return -get_r2() / get_r1();
	}

private:
	double get_b()
	{
		return rho_S * specific_heat_capacity_S;
	}

	double get_p1()
	{
		return (T2-T1)/(rho_S*specific_heat_fusion);
	}

	double get_p2()
	{
		return T1*(1.0 - specific_heat_capacity_S*(T2-T1)/specific_heat_fusion);
	}

	double get_r1()
	{
		return rho_L * specific_heat_capacity_L;
	}

	double get_r2()
	{
		return rho_S * specific_heat_fusion + rho_S * specific_heat_capacity_S * T1 - rho_L * specific_heat_capacity_L * T2;
	}

};
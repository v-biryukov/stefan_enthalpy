#pragma once


struct MaterialInfo
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


	MaterialInfo()
	{
	}
	MaterialInfo(double T1, double T2, double rho_L, double rho_S, double thermal_conductivity_L, double thermal_conductivity_S,
		double specific_heat_fusion, double specific_heat_capacity_L, double specific_heat_capacity_S)
		: T1(T1), T2(T2), rho_L(rho_L), rho_S(rho_S), thermal_conductivity_L(thermal_conductivity_L), thermal_conductivity_S(thermal_conductivity_S),
		specific_heat_fusion(specific_heat_fusion), specific_heat_capacity_L(specific_heat_capacity_L), specific_heat_capacity_S(specific_heat_capacity_S)
	{
	}

    double GetEnthalpyByT(double T) const
	{
		double E;
		if (T < T1)
			E = GetB() * T;
		else if (T < T2)
			E = (T - GetP2()) / GetP1();
		else
			E = GetR1() * T + GetR2();
		return E;
	}

    double GetTByEnthalpy(double enthalpy) const
	{
		double E1 = GetB() * T1;
		double E2 = GetR1() * T2 + GetR2();
		if (enthalpy < E1)
			return enthalpy / GetB();
		else if (enthalpy < E2)
			return enthalpy * GetP1() + GetP2();
		else 
			return (enthalpy - GetR2()) / GetR1();
	}

    double GetThermalConductivityByE(double enthalpy) const
	{
		double E1 = GetB() * T1;
		double E2 = GetR1() * T2 + GetR2();
		if (enthalpy < E1)
			return thermal_conductivity_S;
		else if (enthalpy < E2)
			return ((thermal_conductivity_L - thermal_conductivity_S) * enthalpy + thermal_conductivity_S * E2 - thermal_conductivity_L * E1) / (E2 - E1);
		else
			return thermal_conductivity_L;
	}

    double GetAlpha(double enthalpy) const
	{
		double E1 = GetB() * T1;
		double E2 = GetR1() * T2 + GetR2();
		if (enthalpy < E1)
			return 1.0 / GetB();
		else if (enthalpy < E2)
			return GetP1();
		else
			return 1.0 / GetR1();
	}


    double GetBeta(double enthalpy) const
	{
		double E1 = GetB() * T1;
		double E2 = GetR1() * T2 + GetR2();
		if (enthalpy < E1)
			return 0.0;
		else if (enthalpy < E2)
			return GetP2();
		else
			return -GetR2() / GetR1();
	}

private:
    double GetB() const
	{
		return rho_S * specific_heat_capacity_S;
	}

    double GetP1() const
	{
        return (T2-T1)/(rho_S*specific_heat_fusion);
	}

    double GetP2() const
	{
        return T1*(1.0 - specific_heat_capacity_S*(T2-T1)/specific_heat_fusion);
	}

    double GetR1() const
	{
		return rho_L * specific_heat_capacity_L;
	}

    double GetR2() const
	{
		return rho_S * specific_heat_fusion + rho_S * specific_heat_capacity_S * T1 - rho_L * specific_heat_capacity_L * T2;
	}

};

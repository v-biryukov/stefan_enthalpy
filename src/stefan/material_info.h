#pragma once


struct MaterialInfo
{
	// material melts from T1 to T2; T1~T2
	const double T1;
	const double T2;

	const double rho_L;
	const double rho_S;
	const double thermal_conductivity_L;
	const double thermal_conductivity_S;
	const double specific_heat_fusion;
	const double specific_heat_capacity_L;
	const double specific_heat_capacity_S;

	MaterialInfo()
        : T1(0.0), T2(0.0), rho_L(0.0), rho_S(0.0), thermal_conductivity_L(0.0), thermal_conductivity_S(0.0),
        specific_heat_fusion(0.0), specific_heat_capacity_L(0.0), specific_heat_capacity_S(0.0)
	{
	}
    MaterialInfo(double T1, double T2, double rho_L, double rho_S, double thermal_conductivity_L, double thermal_conductivity_S,
		double specific_heat_fusion, double specific_heat_capacity_L, double specific_heat_capacity_S)
        : T1(T1), T2(T2), rho_L(rho_L), rho_S(rho_S), thermal_conductivity_L(thermal_conductivity_L), thermal_conductivity_S(thermal_conductivity_S),
		specific_heat_fusion(specific_heat_fusion), specific_heat_capacity_L(specific_heat_capacity_L), specific_heat_capacity_S(specific_heat_capacity_S)
	{
		b = 1.0/(rho_S * specific_heat_capacity_S);
        p1 = (T2-T1)/(rho_S*specific_heat_fusion);
        p2 = T1*(1.0 - specific_heat_capacity_S*(T2-T1)/specific_heat_fusion);
		r1 = 1.0/(rho_L * specific_heat_capacity_L);
        r2 = rho_S * specific_heat_fusion + rho_S * specific_heat_capacity_S * T1 - rho_L * specific_heat_capacity_L * T2;
		E1 = T1 / b;
		E2 = T2 / r1 + r2;
	}

	double GetEnthalpyByT(double T) const
	{
		double E;
		if (T < T1)
			E = T / b;
		else if (T < T2)
            E = (T - p2) / p1;
		else
			E = T / r1 + r2;
		return E;
	}

	double GetTByEnthalpy(double enthalpy) const
	{
		if (enthalpy < E1)
			return enthalpy * b;
		else if (enthalpy < E2)
            return enthalpy * p1 + p2;
		else
			return (enthalpy - r2) * r1;
	}

	double GetThermalConductivityByE(double enthalpy) const
	{
		if (enthalpy < E1)
			return thermal_conductivity_S;
		else if (enthalpy < E2)
			return ((thermal_conductivity_L - thermal_conductivity_S) * enthalpy + thermal_conductivity_S * E2 - thermal_conductivity_L * E1) / (E2 - E1);
		else
			return thermal_conductivity_L;
	}

	double GetAlpha(double enthalpy) const
	{
		if (enthalpy < E1)
			return b;
		else if (enthalpy < E2)
            return p1;
		else
			return r1;
	}


	double GetBeta(double enthalpy) const
	{
		if (enthalpy < E1)
			return 0.0;
		else if (enthalpy < E2)
            return p2;
		else
			return -r2 * r1;
	}

private:


    double b;
	double p1, p2;
	double r1, r2;
    double E1, E2;

    /*
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
    */

};

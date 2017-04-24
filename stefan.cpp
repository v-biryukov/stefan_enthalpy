#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>





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


	material_info(double T1, double T2, double rho_L, double rho_S, double thermal_conductivity_L, double thermal_conductivity_S,
		double specific_heat_fusion, double specific_heat_capacity_L, double specific_heat_capacity_S)
		: T1(T1), T2(T2), rho_L(rho_L), rho_S(rho_S), thermal_conductivity_L(thermal_conductivity_L), thermal_conductivity_S(thermal_conductivity_S),
		specific_heat_fusion(specific_heat_fusion), specific_heat_capacity_L(specific_heat_capacity_L), specific_heat_capacity_S(specific_heat_capacity_S)
	{
	}

	double get_enthalpy_by_T(double T)
	{
		if (T < T1)
			return get_b() * T;
		else if (T < T2)
			return get_p1() * T - get_p2();
		else
			return get_r1() * T + get_r2();
	}

	double get_T_by_enthalpy(double enthalpy)
	{
		double E1 = get_b() * T1;
		double E2 = get_r1() * T2 + get_r2();
		if (enthalpy < E1)
			return enthalpy / get_b();
		else if (enthalpy < E2)
			return (enthalpy + get_p2()) / get_p1();
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
			return ((thermal_conductivity_L - thermal_conductivity_S) * enthalpy + 
				thermal_conductivity_S * E2 - thermal_conductivity_L * E1) / (E2 - E1);
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
			return 1.0 / get_p1();
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
			return get_p2() / get_p1();
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
		return rho_S * specific_heat_capacity_S + rho_S * specific_heat_fusion / (T2-T1);
	}

	double get_p2()
	{
		return rho_S * specific_heat_fusion * (T2 + T1) / 2.0 / (T2-T1);
	}

	double get_r1()
	{
		return rho_L * specific_heat_capacity_L;
	}

	double get_r2()
	{
		return rho_S * specific_heat_fusion + (rho_S * specific_heat_capacity_S - rho_L * specific_heat_capacity_L) * T2;
	}

};


class solver2d
{
	int num_x;
	int num_y;
	int number_of_nodes;
	double Lx;
	double Ly;
	double hx;
	double hy;
	material_info mi;

	double * enthalpy_data;
	double * current_iteration_data;
	double * next_iteration_data;



	void set_enthalpy_by_T_data(material_info & mi, const double * temperature_data)
	{
		for (int n = 0; n < num_x; ++n)
			for (int i = 0; i < num_y; ++i)
				enthalpy_data[n + i * num_x] = mi.get_enthalpy_by_T(temperature_data[n + i * num_x]);
	}


	void iterate_jacobi(int axis, const double * enthalpy_data, double * current_iteration_data, double * next_iteration_data, double dt)
	{
		const double * Ex;
		const double * Ey;
		if (axis == 0)
		{
			Ex = current_iteration_data;
			Ey = enthalpy_data;
		}
		else if (axis == 1)
		{
			Ex = enthalpy_data;
			Ey = current_iteration_data;
		}

		for (int n = 1; n < num_x-1; ++n)
			for (int i = 1; i < num_y-1; ++i)
			{
				double k1 = (mi.get_thermal_conductivity_by_E(Ex[n+1 + i*num_x]) + mi.get_thermal_conductivity_by_E(Ex[n + i*num_x])) / 2.0;
				double k2 = (mi.get_thermal_conductivity_by_E(Ex[n + i*num_x]) + mi.get_thermal_conductivity_by_E(Ex[n-1 + i*num_x])) / 2.0;
				double Tn1 = (mi.get_T_by_enthalpy(Ex[n+1 + i*num_x]) - mi.get_T_by_enthalpy(Ex[n + i*num_x])) / hx;
				double Tn2 = (mi.get_T_by_enthalpy(Ex[n-1 + i*num_x]) - mi.get_T_by_enthalpy(Ex[n + i*num_x])) / hx;
				double Lamda1 = (k1 * Tn1 * hy + k2 * Tn2 * hy);
				k1 = (mi.get_thermal_conductivity_by_E(Ey[n + (i+1)*num_x]) + mi.get_thermal_conductivity_by_E(Ey[n + i*num_x])) / 2.0;
				k2 = (mi.get_thermal_conductivity_by_E(Ey[n + i*num_x]) + mi.get_thermal_conductivity_by_E(Ey[n + (i-1)*num_x])) / 2.0;
				Tn1 = (mi.get_T_by_enthalpy(Ey[n + (i+1)*num_x]) - mi.get_T_by_enthalpy(Ey[n + i*num_x])) / hx;
				Tn2 = (mi.get_T_by_enthalpy(Ey[n + (i-1)*num_x]) - mi.get_T_by_enthalpy(Ey[n + i*num_x])) / hx;
				double Lamda2 = (k1 * Tn1 * hy + k2 * Tn2 * hy);
				next_iteration_data[n + i*num_x] = current_iteration_data[n + i*num_x] + (Lamda1 + Lamda2) * dt / (2.0 * hx * hy);
			}
	}


	void iterate_tridiagonal(int axis, const double * enthalpy_data, double * current_iteration_data, double * next_iteration_data, double dt)
	{
		double * a = new double[num_x];
		double * b = new double[num_x];
		double * c = new double[num_x];
		double * f = new double[num_x];
		double * x = new double[num_x];

		int idx1_p, idx1_c, idx1_m;
		int idx2_p, idx2_c, idx2_m;
		int num_1, num_2;
		double h1, h2;

		double * E = current_iteration_data;
		if (axis == 0)
		{
			idx1_p = 1;
			idx1_c = 0;
			idx1_m = -1;
			idx2_p = num_x;
			idx2_c = 0;
			idx2_m = -num_x;
			num_1 = num_x; num_2 = num_y;
			h1 = hx; h2 = hy;
		}
		else if (axis == 1)
		{
			idx1_p = num_x;
			idx1_c = 0;
			idx1_m = -num_x;
			idx2_p = 1;
			idx2_c = 0;
			idx2_m = -1;
			num_1 = num_y; num_2 = num_x;
			h1 = hy; h2 = hx;
		}

		for (int i2 = 1; i2 < num_2-1; ++i2)
		{
			a[0] = 0.0; b[0] = -1.0; c[0] = 1.0; f[0] = 0.0;
			a[num_1-1] = -1.0; b[num_1-1] = 1.0; c[num_1-1] = 0.0; f[num_1-1] = 0.0;
			

			for (int i1 = 1; i1 < num_1-1; ++i1)
			{
				int ix = axis ? i2 : i1;
				int iy = axis ? i1 : i2;
				double * E_iter = current_iteration_data + ix + iy*num_x;
				const double * E_step = enthalpy_data + ix + iy*num_x;

				double kp = (mi.get_thermal_conductivity_by_E(E_iter[idx1_p]) + mi.get_thermal_conductivity_by_E(E_iter[idx1_c])) / 2.0;
				double km = (mi.get_thermal_conductivity_by_E(E_iter[idx1_c]) + mi.get_thermal_conductivity_by_E(E_iter[idx1_m])) / 2.0;


				double coef = -dt / (2.0 * h1 * h1);
				a[i1] = coef * mi.get_alpha(E_iter[idx1_m]) * km;
				b[i1] = 1.0 - coef * mi.get_alpha(E_iter[idx1_c]) * (km + kp);
				c[i1] = coef * mi.get_alpha(E_iter[idx1_p]) * kp;
				f[i1] = E_step[idx1_c] - coef * (mi.get_beta(E_iter[idx1_m])*km - mi.get_beta(E_iter[idx1_c])*(km+kp) + mi.get_beta(E_iter[idx1_p])*kp);


				double k2p = (mi.get_thermal_conductivity_by_E(E_step[idx2_p]) + mi.get_thermal_conductivity_by_E(E_step[idx2_c])) / 2.0;
				double k2n = (mi.get_thermal_conductivity_by_E(E_step[idx2_m]) + mi.get_thermal_conductivity_by_E(E_step[idx2_c])) / 2.0;
				double T2p = mi.get_T_by_enthalpy(E_step[idx2_p]) - mi.get_T_by_enthalpy(E_step[idx2_c]);
				double T2n = mi.get_T_by_enthalpy(E_step[idx2_m]) - mi.get_T_by_enthalpy(E_step[idx2_c]);
				f[i1] += k2p * T2p * dt/(2*h2*h2) + k2n * T2n * dt/(2*h2*h2);
			}

			solve_triadiagonal(num_1, a, b, c, f, x);

			if (axis == 0)
				for (int k = 0; k < num_1; ++k)
					next_iteration_data[k + i2*num_x] = x[k];
			else if (axis == 1)
				for (int k = 0; k < num_1; ++k)
					next_iteration_data[i2 + k*num_x] = x[k];
		}

		delete [] a;
		delete [] b;
		delete [] c;
		delete [] f;
		delete [] x;
	}

	void solve_triadiagonal (int n, double *a, double *b, double *c, double *f, double *x)
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
		double L2 = 0.0;
		for (int n = 0; n < num_x; ++n)
			for (int i = 0; i < num_y; ++i)
			{
				L2 += (A[n + i*num_x] - B[n + i*num_x]) * (A[n + i*num_x] - B[n + i*num_x]);
			}
		L2 = sqrt(L2 / number_of_nodes);
		return L2;
	}

	double calculate_L2_norm(const double * A)
	{
		double L2 = 0.0;
		for (int n = 0; n < num_x; ++n)
			for (int i = 0; i < num_y; ++i)
			{
				L2 += A[n + i*num_x] * A[n + i*num_x];
			}
		L2 = sqrt(L2 / number_of_nodes);
		return L2;
	}

	double assign(double * A, const double * B)
	{
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
		while (L2_norm > 5e-3)
		{
			assign(current_iteration_data, next_iteration_data);
			iterate_tridiagonal(axis, enthalpy_data, current_iteration_data, next_iteration_data, dt);
			L2_norm = calculate_L2_norm(current_iteration_data, next_iteration_data) / calculate_L2_norm(current_iteration_data);
			//std::cout << "L2 = " << L2_norm << std::endl;
		}
		assign(enthalpy_data, next_iteration_data);
	}



public:
	solver2d(int num_x, int num_y, double Lx, double Ly, material_info & mi, const double * temperature_data) 
		: num_x(num_x), num_y(num_y), Lx(Lx), Ly(Ly), mi(mi)
	{
		hx = Lx / (num_x-1);
		hy = Ly / (num_y-1);
		number_of_nodes = num_x * num_y;
		enthalpy_data = new double[number_of_nodes];
		current_iteration_data = new double [number_of_nodes];
		next_iteration_data = new double [number_of_nodes];
		set_enthalpy_by_T_data(mi, temperature_data);
	}

	~solver2d()
	{
		delete [] enthalpy_data;
		delete [] current_iteration_data;
		delete [] next_iteration_data;
	}

	void step(double dt)
	{
		step(0, dt);
		step(1, dt);
	}

	void save_to_vtk(std::string name)
	{
	    std::ofstream vtk_file(name.c_str());
	    vtk_file << "# vtk DataFile Version 3.0\nVx data\nASCII\n\n";
	    vtk_file << "DATASET POLYDATA\nPOINTS " << number_of_nodes << " float\n";
	    for ( int i = 0; i < number_of_nodes; i++ )
	        vtk_file << (i/num_x) * hx << " " << (i%num_x) * hy << " "  << 0.0 << "\n";
	    vtk_file << "\nPOLYGONS " << (num_x-1)*(num_y-1) << " " << (num_x-1)*(num_y-1)*5 << "\n";
	    for (int i = 0; i < num_x-1; i++)
			for (int j = 0; j < num_y-1; j++)
	        	vtk_file << 4 << " " << i+j*num_x << " " << i+1+j*num_x << " " << i+1+(j+1)*num_x << " " << i+(j+1)*num_x << "\n";
	    vtk_file << "\nPOINT_DATA " << number_of_nodes << "\n";
	    vtk_file << "SCALARS E FLOAT\n";
	    vtk_file << "LOOKUP_TABLE default\n";
	    for (int i = 0; i < num_x; i++)
			for (int j = 0; j < num_y; j++)
			{
	        	vtk_file <<  enthalpy_data[i + j*num_x] << "\n";
	    	}
	    vtk_file << "SCALARS T FLOAT\n";
	    vtk_file << "LOOKUP_TABLE default\n";
	    for (int i = 0; i < num_x; i++)
			for (int j = 0; j < num_y; j++)
			{
	        	vtk_file <<  mi.get_T_by_enthalpy(enthalpy_data[i + j*num_x]) << "\n";
	    	}
	   	vtk_file << "SCALARS K FLOAT\n";
	    vtk_file << "LOOKUP_TABLE default\n";
	    for (int i = 0; i < num_x; i++)
			for (int j = 0; j < num_y; j++)
			{
	        	vtk_file <<  mi.get_thermal_conductivity_by_E(enthalpy_data[i + j*num_x]) << "\n";
	    	}
	    vtk_file << "SCALARS state FLOAT\n";
	    vtk_file << "LOOKUP_TABLE default\n";
	    for (int i = 0; i < num_x; i++)
			for (int j = 0; j < num_y; j++)
			{
				double T = mi.get_T_by_enthalpy(enthalpy_data[i + j*num_x]);
				double state;
				if (T < mi.T1)
					state = 0.0;
				else if (T < mi.T2)
					state = (T-mi.T1)/(mi.T2-mi.T1);
				else
					state = 1.0;
	        	vtk_file <<  state << "\n";
	    	}
	}


};

int main()
{
	int num_x = 100;
	int num_y = 100;
	material_info mi = material_info(273.0, 274.0, 1000.0, 920.0, 0.591, 2.22, 334000.0, 4200.0, 2100.0);
	double * td = new double [num_x * num_y];
	for (int n = 0; n < num_x; ++n)
		for (int i = 0; i < num_y; ++i)
		{
			//if ((n > 15 && n < 25 && i > 10 && i < 90) || (n > 45 && n < 55 && i > 10 && i < 90) || (n > 75 && n < 85 && i > 10 && i < 90)
			//	|| (n > 10 && n < 90 && i > 15 && i < 25) || (n > 10 && n < 90 && i > 45 && i < 55) || (n > 10 && n < 90 && i > 75 && i < 85))
			if ((n > 30 && n < 70 && i > 30 && i < 70))
			{
				td[n + i*num_x] = 70.0;
			}
			else
			{
				td[n + i*num_x] = 300.0;
			}
		}
	solver2d sol = solver2d(num_x, num_y, 5.0, 5.0, mi, td);
	for (int i = 0; i < 100000; ++i)
	{
		if (i%10 == 0)
		{
			std::cout << "step: " << i << std::endl;
			std::stringstream ss;
	        ss << i;
	        sol.save_to_vtk("out/result_" + ss.str() + ".vtk");
    	}
		sol.step(1000);	
	}
	delete [] td;
	return 0;
}

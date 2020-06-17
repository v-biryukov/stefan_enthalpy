#pragma once
#include <cmath>

template <int Dims>
class Solver;

template <int Dims>
struct SubmeshEffect
{
    virtual void ApplyEffect(Solver<Dims>& sol, int global_id)
    {
    }
};


template <int Dims>
struct FixedTemperatureEffect : public SubmeshEffect<Dims>
{
    double temperature;

    void ApplyEffect(Solver<Dims>& sol, int global_id)
    {
        sol.enthalpy_data[global_id] = sol.mesh.GetMaterialInfo(global_id).GetEnthalpyByT(temperature);
    }
};


template <int Dims>
struct SineTemperatureEffect : public SubmeshEffect<Dims>
{
    double average;
    double amplitude;
    double period;
    double phase_shift;

    void ApplyEffect(Solver<Dims>& sol, int global_id)
    {
        double temperature = average + amplitude * sin(sol.time / period + phase_shift);
        sol.enthalpy_data[global_id] = sol.mesh.GetMaterialInfo(global_id).GetEnthalpyByT(temperature);
    }
};



template <int Dims>
struct ChangeIndexOnMeltingEffect : public SubmeshEffect<Dims>
{
    int index;
    void ApplyEffect(Solver<Dims>& sol, int global_id)
    {
        double current_temperature = sol.mesh.GetMaterialInfo(global_id).GetTByEnthalpy(sol.enthalpy_data[global_id]);
        double melting_temperature = sol.GetMesh().GetMaterialInfo(global_id).T2;

        if (current_temperature > melting_temperature + 1e-6)
        {
            sol.mesh.SetMaterialIndex(global_id, index);
        }
    }
};

template <int Dims>
struct ChangeIndexOnFreezingEffect : public SubmeshEffect<Dims>
{
    int index;
    void ApplyEffect(Solver<Dims>& sol, int global_id)
    {
        double current_temperature = sol.mesh.GetMaterialInfo(global_id).GetTByEnthalpy(sol.enthalpy_data[global_id]);
        double freezing_temperature = sol.GetMesh().GetMaterialInfo(global_id).T1;

        if (current_temperature < freezing_temperature - 1e-6)
        {
            sol.mesh.SetMaterialIndex(global_id, index);
        }
    }
};
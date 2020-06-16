#pragma once


struct SubmeshEffect
{
    enum class EffectType
    {
        FixedTemperature,
        ChangeIndexOnMelting
    };
    EffectType type;

    virtual void ApplyEffect()
    {
    }
};
struct FixedTemperatureEffect : public SubmeshEffect
{
    double temperature;
};
struct ChangeIndexOnMeltingEffect : public SubmeshEffect
{
    int index;
};
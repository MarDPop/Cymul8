#pragma once

#include "Constants.h"
#include "Coordinates.h"

struct Air
{
    static constexpr double R_GAS = 8.31446261815324; // J / K mol exact
    static constexpr double ATM = 101325; // Pa
    static constexpr double SL_TEMPERATURE = CONSTANTS::ABSOLUTE_ZERO_KELVIN + 15; // Kelvin
    static constexpr double DRY_AIR_MOLAR_MASS = 0.0289652; // kg/ mol
    static constexpr double SPECIFIC_HEAT_RATIO = 1.4;
    static constexpr double R_DRY_AIR = R_GAS / DRY_AIR_MOLAR_MASS;

    double pressure = ATM; // Pa
    double inv_sound_speed = 1.0 / sqrt(SPECIFIC_HEAT_RATIO*R_DRY_AIR*SL_TEMPERATURE); // s/m
    double dynamic_viscosity = 1.803e-7; // Pa s
    double temperature = SL_TEMPERATURE; // K
    double density = ATM / (R_DRY_AIR * SL_TEMPERATURE); // kg/m3
    double gamma = SPECIFIC_HEAT_RATIO; // ratio of specific heats
    double molar_mass = DRY_AIR_MOLAR_MASS; // kg/ m3 dry air
    double R_specific = R_DRY_AIR;
};

class Atmosphere
{

public:

    static double geometric2geopotential(double z, const double R0 = 6371)
    {
        return R0 * z / (R0 + z);
    }

    static double geopotential2geometric(double H, const double R0 = 6371)
    {
        return R0 * H / (R0 - H);
    }

    virtual void set_time(double time) {}

    virtual void init(Air& air) const {}

    virtual void update(const Coordinate::Geodetic& lla,
        Air& air, 
        Eigen::Vector3d& wind) const {}

};

class Atmosphere_Const_Temperature_Const_Gas : public virtual Atmosphere
{

public:    

    const double GAMMA;

    const double MW;

    const double R_SPECIFIC;

    const double PRES_SL;

    const double TEMP_SL;

    const double DENS_SL;

    const double G_SL;

    const double INV_SOUND_SPEED;

    const double SCALE_FACTOR;

    const double DYNAMIC_VISC;

    Atmosphere_Const_Temperature_Const_Gas(double pres_sl,
        double temp_sl,
        double gamma = 1.4,
        double mw = 0.02886,
        double g0 = 9.806,
        double mu_vis = 1.803e-7) :
            PRES_SL(pres_sl),
            GAMMA(gamma),
            MW(mw),
            R_SPECIFIC(Air::R_GAS / mw),
            TEMP_SL(temp_sl),
            DENS_SL(pres_sl / (R_SPECIFIC * temp_sl)),
            G_SL(g0),
            INV_SOUND_SPEED(1.0 / sqrt(gamma * R_SPECIFIC * temp_sl)),
            SCALE_FACTOR(-g0/(R_SPECIFIC*temp_sl)),
            DYNAMIC_VISC(mu_vis)
    {}

    void init(Air& air) const
    {
        air.temperature = TEMP_SL;
        air.inv_sound_speed = INV_SOUND_SPEED;
        air.gamma = GAMMA;
        air.molar_mass = MW;
        air.dynamic_viscosity = DYNAMIC_VISC;
        air.R_specific = R_SPECIFIC;
    }

    void update(const Coordinate::Geodetic& lla,
        Air& air,
        Eigen::Vector3d& wind) const override
    {
        double factor = exp(SCALE_FACTOR * lla.altitude);
        air.pressure = PRES_SL * factor;
        air.density = DENS_SL * factor;
    }

};

class Atmosphere_ISA : public virtual Atmosphere
{

};
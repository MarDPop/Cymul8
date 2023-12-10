#pragma once

#include "Constants.h"
#include "Coordinates.h"
#include "../util/Table.h"
#include <array>
#include <vector>

struct Air
{
    static constexpr double R_GAS = 8.31446261815324; // J / K mol exact
    static constexpr double ATM = 101325; // Pa
    static constexpr double SL_TEMPERATURE = CONSTANTS::ABSOLUTE_ZERO_KELVIN + 15; // Kelvin
    static constexpr double DRY_AIR_MOLAR_MASS = 0.0289652; // kg/ mol
    static constexpr double SPECIFIC_HEAT_RATIO = 1.4;
    static constexpr double R_DRY_AIR = R_GAS / DRY_AIR_MOLAR_MASS;
    static constexpr double WATER_MOLAR_MASS = 0.01801; // kg /mol

    double pressure = ATM; // Pa
    double inv_sound_speed = 1.0 / sqrt(SPECIFIC_HEAT_RATIO*R_DRY_AIR*SL_TEMPERATURE); // s/m
    double dynamic_viscosity = 1.803e-5; // Pa s
    double temperature = SL_TEMPERATURE; // K
    double density = ATM / (R_DRY_AIR * SL_TEMPERATURE); // kg/m3
    double gamma = SPECIFIC_HEAT_RATIO; // ratio of specific heats
    double molar_mass = DRY_AIR_MOLAR_MASS; // kg/ m3 dry air
    double R_specific = R_DRY_AIR;

    /**
     * Calculates the isothermal pressure ratio from base geopotential height
     * @param kelvin
     * @return dynamic pressure in Pa s
     */
    static double dynamic_viscosity_sutherland(double kelvin)
    {
        return (5.058e-8 - 1.295e-11 * kelvin) * kelvin + 4.371e-6;
    }

    /**
     * Gets the partial pressure of water at temperature
     * @param kelvin
     * @return pressure in Pa
     */
    static double partial_pressure_H20(double kelvin)
    {
        double celsius = kelvin - CONSTANTS::ABSOLUTE_ZERO_KELVIN;
        return 610.76*exp(17.27*celsius/(celsius+237.5));
    }
};

class Wind
{
    virtual void update(Eigen::Vector3d& wind) const = 0;
};

class Atmosphere
{

public:

    /**
     * Calculates the isothermal pressure ratio from base geopotential height
     * @param z
     * @param R0
     * @return
     */
    static double geometric2geopotential(double z, const double R0 = 6371)
    {
        return R0 * z / (R0 + z);
    }

    /**
     * Calculates the isothermal pressure ratio from base geopotential height
     * @param H
     * @param R0
     * @return
     */
    static double geopotential2geometric(double H, const double R0 = 6371)
    {
        return R0 * H / (R0 - H);
    }


    /**
     * Calculates the isothermal pressure ratio from base geopotential height
     * @param deltaH
     * @param temperature
     * @return
     */
    static double isothermalPressureRatio(double deltaH, 
        double temperature, 
        double R_GAS, 
        double g0)
    {
        return exp(-g0 * deltaH / (R_GAS * temperature));
    }

    /**
     * Calculates the isothermal pressure ratio from base geopotential height
     * @param deltaH
     * @param temperature
     * @return
     */
    static double gradientPressureRatio(double lapseRate, 
        double T1, 
        double T2,
        double R_GAS,
        double g0)
    {
        return pow(T2 / T1, -g0 / (R_GAS*lapseRate));
    }

    virtual void set_time(double time) {}

    virtual void init(Air& air) const {}

    virtual void update_air(const Coordinate::Geodetic& lla,
        Air& air) const {}

};

class Atmosphere_Const_Temperature_Const_Gas final : public virtual Atmosphere
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

    void init(Air& air) const override
    {
        air.temperature = TEMP_SL;
        air.inv_sound_speed = INV_SOUND_SPEED;
        air.gamma = GAMMA;
        air.molar_mass = MW;
        air.dynamic_viscosity = DYNAMIC_VISC;
        air.R_specific = R_SPECIFIC;
    }

    void update_air(const Coordinate::Geodetic& lla,
        Air& air) const override
    {
        double factor = exp(SCALE_FACTOR * geometric2geopotential(lla.altitude));
        air.pressure = PRES_SL * factor;
        air.density = DENS_SL * factor;
    }

};

class Atmosphere_Table final : public virtual Atmosphere
{
    BasicTable<double, 5> _table;

    const double _gamma;

    const double _mw;

    const double _R;

public:

    Atmosphere_Table(const std::vector<double>& height,
        const std::vector<std::array<double, 5>>& values,
        double gamma,
        double mw,
        double R) :
            _table(height, values),
            _gamma(gamma),
            _mw(mw),
            _R(R) {}

    void init(Air& air) const override
    {
        air.gamma = _gamma;
        air.molar_mass = _mw;
        air.R_specific = _R;
    }

    void update_air(const Coordinate::Geodetic& lla,
        Air& air) const override
    {
        _table.get(lla.altitude, &air.pressure);
    }

};


class Atmosphere_Table_VariableMixture final : public virtual Atmosphere
{
    BasicTable<double, 8> _table;

public:

    Atmosphere_Table_VariableMixture(const std::vector<double>& height,
        const std::vector<std::array<double, 8>>& values) :
            _table(height, values) {}

    void update_air(const Coordinate::Geodetic& lla,
        Air& air) const override
    {
        _table.get(lla.altitude, &air.pressure);
    }

};

class Atmosphere_MILSTD210 final : public virtual Atmosphere
{
    static constexpr unsigned LOWER_TABLE_SIZE = 101;
    static constexpr unsigned UPPER_TABLE_SIZE = 200;
    
    Fixed_Incement_Table_Zero_Based<double, LOWER_TABLE_SIZE, 8> _lower_table;

    Fixed_Incement_Table_Zero_Based<double, UPPER_TABLE_SIZE, 8> _upper_table;

    double _scale_factor;

public:

    static constexpr double LOWER_HEIGHT_INCREMENT = 304.8; // meters
    static constexpr double LOWER_TABLE_HEIGHT = 30480.0; // meters
    static constexpr double UPPER_HEIGHT_INCREMENT = LOWER_HEIGHT_INCREMENT * 2;

    static const double STANDARD[LOWER_TABLE_SIZE];
    static const double HOT[LOWER_TABLE_SIZE];
    static const double COLD[LOWER_TABLE_SIZE];
    static const double POLAR[LOWER_TABLE_SIZE];
    static const double TROPICAL[LOWER_TABLE_SIZE];

    static const double STANDARD_HUMIDITY[LOWER_TABLE_SIZE];
    static const double HOT_HUMIDITY[LOWER_TABLE_SIZE];
    static const double COLD_HUMIDITY[LOWER_TABLE_SIZE];
    static const double POLAR_HUMIDITY[LOWER_TABLE_SIZE];
    static const double TROPICAL_HUMIDITY[LOWER_TABLE_SIZE];

    inline static constexpr double ATMOSPHERE_LAYER_HEIGHTS[8] = {11000, 20000, 32000, 470000, 51000, 71000, 84000, 90000};

    static constexpr double TROPOSPHERE_LAPSE_RATE = -0.0065;
    static constexpr double STRATOSPHERE_LAPSE_RATE1 = 0.001;
    static constexpr double STRATOSPHERE_LAPSE_RATE2 = 0.0028;
    static constexpr double MESOSPHERE_LAPSE_RATE1 = -0.0028;
    static constexpr double MESOSPHERE_LAPSE_RATE2 = -0.002;
    static constexpr double THERMOSPHERE_LAPSE_RATE2 = 0.004;

    Atmosphere_MILSTD210() :
        _lower_table(LOWER_HEIGHT_INCREMENT),
        _upper_table(UPPER_HEIGHT_INCREMENT) {}

    enum class DAY_TYPE
    {
        STANDARD = 0,
        HOT,
        COLD,
        POLAR,
        TROPICAL
    };

    void set(DAY_TYPE type,
        double temp_offset,
        double humidity_offset,
        double p0,
        double g0);

    void update_air(const Coordinate::Geodetic& lla,
        Air& air) const override
    {
        double geopotential_h = geometric2geopotential(lla.altitude);
        if (geopotential_h < LOWER_TABLE_HEIGHT)
        {
            _lower_table.get(geopotential_h, &air.pressure);
        }
        else
        {
            double h = geopotential_h - LOWER_TABLE_HEIGHT;
            _upper_table.get(h, &air.pressure);
        }
        
    }

};

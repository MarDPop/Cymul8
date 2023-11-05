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

    virtual void update(const Coordinate::Geodetic& lla,
        double time, 
        Air& air, 
        Eigen::Vector3d& wind) {}

};
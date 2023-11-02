#pragma once

#include <memory>
#include "../physics/Coordinates.h"
#include "../physics/SolarSystem.h"

class Environment
{
public:

    enum STATE
    {
        NONE = -1,
        LAUNCH_LANDING = 0,
        ATMOSPHERIC,
        LOW_ORBIT,
        HIGH_ORBIT,
        INTERPLANETARY,
        COASTING
    };

private:

    SolarSystemBody* _current_planet;

    double _TALO;

    Coordinate::GeocentricInertial _ECI;

    Coordinate::GeocentricFixed _ECEF;

    Coordinate::Geodetic _LLA;

    Coordinate::ENU _ENU;

    Coordinate::Spherical _RTP;

    EpochTime _launch_time;

    Coordinate::Geodetic _launch_lla;

    Eigen::Vector3d _frame_acceleration; // includes gravity

    Air* _air;

    double _radiation_pressure;

    STATE _current_state;

    void update_launch(const double* state,
        double talo);

    void update_atm(const double* state,
        double talo);

    void update_low_orbit(const double* state,
        double talo);

    void update_high_orbit(const double* state,
        double talo);

    void update_interplanetary(const double* state,
        double talo);

    void update_coast(const double* state,
        double talo);

public:
    
    void set_launch(EpochTime launch_time,
                    Coordinate::Geodetic launch)
    {
        _launch_time = launch_time;
        _launch_lla = launch;
    }

    void set_air_reference(Air* air)
    {
        _air = air;
    }

    void set_radation_pressure(double radiation_pressure)
    {
        _radiation_pressure = radiation_pressure;
    }

    void set_state(STATE state);

    void (Environment::* update)(const double*, double);

};
#pragma once

#include <memory>
#include "../physics/Coordinates.h"
#include "../physics/SolarSystem.h"

struct AeroData
{
    Air air;
    double airspeed;
    double mach;
    double dynamic_pressure;
    double impact_pressure;
};

struct LaunchData
{
    Coordinate::Geodetic launch_lla;
    EpochTime launch_time;
    SolarSystemBody* launch_body;
};

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

    Eigen::Vector3d _frame_acceleration; // includes gravity

    AeroData _aero;

    double _radiation_pressure;

    std::unique_ptr<LaunchData> _launch;

    STATE _current_state;

    void update_launch( const Eigen::Vector3d& position, 
                        const Eigen::Vector3d& velocity,
                        double talo);

    void update_atm(const Eigen::Vector3d& position,
        const Eigen::Vector3d& velocity,
        double talo);

    void update_low_orbit(const Eigen::Vector3d& position,
        const Eigen::Vector3d& velocity,
        double talo);

    void update_high_orbit(const Eigen::Vector3d& position,
        const Eigen::Vector3d& velocity,
        double talo);

    void update_interplanetary(const Eigen::Vector3d& position,
        const Eigen::Vector3d& velocity,
        double talo);

    void update_coast(const Eigen::Vector3d& position,
        const Eigen::Vector3d& velocity,
        double talo);

public:
    
    void set_launch(EpochTime launch_time,
                    Coordinate::Geodetic launch)
    {
        _launch_time = launch_time;
        _launch_lla = launch;
    }

    void set_radation_pressure(double radiation_pressure)
    {
        _radiation_pressure = radiation_pressure;
    }

    void set_state(STATE state);

    void (Environment::*update)(const Eigen::Vector3d&, const Eigen::Vector3d&, double);

};
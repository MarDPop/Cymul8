#pragma once

#include <memory>
#include "../physics/Coordinates.h"
#include "../physics/SolarSystem.h"
#include "../physics/Geometry.h"

struct AeroData
{
    Eigen::Vector3d air_velocity_ecef;
    double airspeed;
    double mach;
    double dynamic_pressure;
    double impact_pressure;
};

class Environment
{
public:

    enum class STATE
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

    EpochTime _reference_time;

    Coordinate::GeocentricFixed _reference_ground;

    double _TALO;

    Coordinate::GeocentricInertial _ECI;

    Coordinate::GeocentricFixed _ECEF;

    Coordinate::Geodetic _LLA;

    Coordinate::ENU _ENU;

    Coordinate::Spherical _RTP;

    Eigen::Vector3d _frame_acceleration; // includes gravity

    Air _air;

    AeroData _aero_data;

    double _radiation_pressure = 0.0;

    STATE _current_state = STATE::NONE;

    bool _use_aero = true;

    void update_launch_landing( const Eigen::Vector3d& position, 
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
        _reference_time = launch_time;
        _reference_ground = WGS84::LLA2ECEF(launch); // in future make some function for geodetic to ecef for each planet
    }

    void set_radation_pressure(double radiation_pressure)
    {
        _radiation_pressure = radiation_pressure;
    }

    const Air& get_air() const
    {
        return _air;
    }

    const AeroData& get_aero_data() const
    {
        return _aero_data;
    }

    const Eigen::Vector3d get_frame_acceleration() const
    {
        return _frame_acceleration;
    }

    bool in_air() const
    {
        return _use_aero;
    }

    void set_state(STATE state);

    void (Environment::*update)(const Eigen::Vector3d&, const Eigen::Vector3d&, double);

};
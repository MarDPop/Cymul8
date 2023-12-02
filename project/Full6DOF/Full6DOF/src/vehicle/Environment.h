#pragma once

#include <memory>
#include "../physics/Coordinates.h"
#include "../physics/SolarSystem.h"
#include "../physics/Geometry.h"


struct AeroData
{
    Air air;
    /**
    * @brief direction of the vehicle against the air in the ecef frame
    */
    Eigen::Vector3d air_velocity_ecef_unit;
    double airspeed;
    double mach;
    double dynamic_pressure;
    double impact_pressure;

    void compute(const Air& air, 
        const Eigen::Vector3d& air_ecef_velocity);
};

struct GroundAndTimeReference
{
    double jd2000_utc_launch;

    long unix_ns_launch;

    double TALO;

    double jd2000_ut1;

    Coordinate::GeocentricFixed ground;
};

struct PlanetCoordinates
{
    Coordinate::GeocentricInertial PCI;

    Coordinate::GeocentricFixed PCF;

    Coordinate::Geodetic LLA;

    Coordinate::ENU ENU;

    Coordinate::Spherical RTP;
};

struct Environment
{
    Eigen::Vector3d frame_acceleration; // includes gravity

    AeroData aero_data;
};

class LocalFrame
{
    SolarSystemBody* _current_planet;

    Coordinate::GeocentricFixed _ref;

    long _unix_time;

    PlanetCoordinates _coordinates;

public:

    void update(double* pos_vel, double talo,
        Environment& environment);

};


class System
{
private:

    Planetary_System* _planetary_system;

    SolarSystemBody* _current_planet;


public:   

    void set_launch(Time::UNIX_TIMESTAMP launch_time,
        Coordinate::Geodetic launch)
    {
        Time::EpochDate jd2000_date = Time::to_epoch_date_j2000(launch_time);
        _ref.jd2000_utc_launch = jd2000_date.get_day_number() + jd2000_date.get_day_fraction();
        _ref.unix_ns_launch = launch_time;
        _ref.TALO = 0;
        _ref.ground = WGS84::LLA2ECEF(launch); // in future make some function for geodetic to ecef for each planet
    }

    void set_current_planet(SolarSystemBody* current_planet)
    {
        _current_planet = current_planet;
    }

    const SolarSystemBody& get_current_planet() const
    {
        return *_current_planet;
    }

    const GroundAndTimeReference& get_references() const
    {
        return _ref;
    }

    const Eigen::Vector3d& get_frame_acceleration() const
    {
        return _frame_acceleration;
    }

    void update(const Eigen::Vector3d& position, 
        const Eigen::Vector3d& velocity, 
        double time,
        bool near_body = true);

};


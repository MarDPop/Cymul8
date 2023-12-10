#pragma once

#include <memory>
#include "../physics/Coordinates.h"
#include "../physics/SolarSystem.h"
#include "../physics/Geometry.h"
#include "../util/Table.h"

struct FrameReference
{
    double jd2000_utc_launch;

    long unix_ns;

    Coordinate::GeocentricFixed point;

    static FrameReference from_launch(const Coordinate::Geodetic& lla,
        const Time::Gregorian& launch_time);
};

struct PlanetCoordinates
{
    Coordinate::GeocentricFixed PCF;

    Coordinate::GeocentricFixed PCF_velocity;

    Coordinate::Geodetic LLA;

    Coordinate::ENU ENU;

    Coordinate::Spherical RTP;
};

struct AeroEnvironment
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

struct ElectroMagneticEnvironment
{
    Eigen::Vector3d magnetic_field;

    Eigen::Vector3d electric_field;

    Eigen::Vector3d radiation_pressure;

    double charge_density;

    double total_radiant_intensity;

    XYTable spectral_intensity;

};

class PlanetaryFrame
{
    PlanetaryBody* _current_planet;

    double _local_to_TT;

    PlanetCoordinates _coordinates;

public:

    void update(const Eigen::Vector3d& pos,
        const Eigen::Vector3d& vel,
        double talo);

    void update(const FrameReference& frame_ref,
        const Eigen::Vector3d& pos,
        const Eigen::Vector3d& vel,
        double talo);

    void set_current_planet(PlanetaryBody* current_planet)
    {
        _current_planet = current_planet;
    }

    const PlanetaryBody& get_current_planet() const
    {
        return *_current_planet;
    }

    const PlanetCoordinates& get_planet_coordinates() const
    {
        return _coordinates;
    }

};

class ICRFFrame
{

    SolarSystem* _icrf;

    Eigen::Vector3d _ICRF_position;

    Eigen::Matrix3d _ICRF_orientation;

public:

    const Eigen::Vector3d& get_ICRF_position() const
    {
        return _ICRF_position;
    }

    const Eigen::Matrix3d& get_ICRF_orientation() const
    {
        return _ICRF_orientation;
    }
    
};


class Environment
{
    Eigen::Vector3d _frame_acceleration; // includes gravity

    AeroEnvironment* _aero_environment = nullptr;

    ElectroMagneticEnvironment* _em_environment = nullptr;

public:

    const Eigen::Vector3d& get_frame_acceleration() const
    {
        return _frame_acceleration;
    }

    const AeroEnvironment* get_aero() const
    {
        return _aero_environment;
    }

    const ElectroMagneticEnvironment* get_em() const
    {
        return _em_environment;
    }

    void update_from_planetary_frame(const PlanetaryFrame& planetaryFrame,
        const Eigen::Vector3d& pos, 
        const Eigen::Vector3d& vel, 
        double time);

    void update_icrf_frame(const ICRFFrame& icrfFrame,
        const Eigen::Vector3d& pos,
        const Eigen::Vector3d& vel,
        double time);
};
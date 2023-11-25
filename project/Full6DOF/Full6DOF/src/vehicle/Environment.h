#pragma once

#include <memory>
#include "../physics/Coordinates.h"
#include "../physics/SolarSystem.h"
#include "../physics/Geometry.h"

enum class SIMULATION_STATE
{
    NONE = -1,
    LAUNCH_LANDING = 0,
    ATMOSPHERIC,
    LOW_ORBIT,
    HIGH_ORBIT,
    INTERPLANETARY,
    COASTING
};

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
    EpochTime time;

    double TALO;

    Coordinate::GeocentricFixed ground;

    void set(EpochTime launch_time,
        Coordinate::Geodetic launch)
    {
        time = launch_time;
        ground = WGS84::LLA2ECEF(launch); // in future make some function for geodetic to ecef for each planet
    }
};

class Environment
{
protected:

    SolarSystemBody* _current_planet;

    GroundAndTimeReference _ref;

    Coordinate::GeocentricInertial _PCI;

    Coordinate::GeocentricFixed _PCF;

    Eigen::Vector3d _frame_acceleration; // includes gravity

public:

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

    const Coordinate::GeocentricInertial& get_PCI() const
    {
        return _PCI;
    }

    const Coordinate::GeocentricFixed& get_PCF() const
    {
        return _PCF;
    }

    const Eigen::Vector3d get_frame_acceleration() const
    {
        return _frame_acceleration;
    }

    void update(const Eigen::Vector3d& position, 
        const Eigen::Vector3d& velocity, 
        double time);

};

class Environment_NearBody 
{
    Coordinate::Geodetic _LLA;

    Coordinate::ENU _ENU;

    Coordinate::Spherical _RTP;

    AeroData _aero_data;

public:

    void update(const Environment& environment);

    const Coordinate::Geodetic& get_LLA() const
    {
        return _LLA;
    }

    const Coordinate::ENU& get_ENU() const
    {
        return _ENU;
    }

    const Coordinate::Spherical& get_RTP() const
    {
        return _RTP;
    }

    const AeroData& get_aero_data() const
    {
        return _aero_data;
    }
};
#pragma once

#include <memory>
#include "../physics/Coordinates.h"
#include "../physics/SolarSystem.h"
#include "../physics/Geometry.h"
#include "../util/Table.h"

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

    double charge_density;

    double total_radiant_intensity;

    XYTable spectral_intensity;

};

struct Environment
{
    Eigen::Vector3d frame_acceleration; // includes gravity

    AeroEnvironment* aero_environment = nullptr;

    ElectroMagneticEnvironment* em_environment = nullptr;
};

class LocalFrame
{
    SolarSystemBody* _current_planet;

    Coordinate::GeocentricFixed _ref;

    double _jd;

    double _local_to_TDB;

    PlanetCoordinates _coordinates;

public:

    void update(const Coordinate::GeocentricFixed& new_ref,
        double julianDate);

    void update_coordinates(double* pos_vel,
        double talo);

    void set_current_planet(SolarSystemBody* current_planet)
    {
        _current_planet = current_planet;
    }

    const SolarSystemBody& get_current_planet() const
    {
        return *_current_planet;
    }

    const PlanetCoordinates& get_planet_coordinates() const
    {
        return _coordinates;
    }

};

class SystemFrame
{

    SolarSystem* _solar_system;

    Eigen::Vector3d _position;

    Eigen::Matrix3d _orientation;

public:

};


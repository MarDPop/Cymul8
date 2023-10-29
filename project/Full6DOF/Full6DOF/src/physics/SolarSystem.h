#pragma once

#include "Time.h"
#include "Gravity.h"
#include "Atmosphere.h"
#include "Ephemeris.h"
#include "Geometry.h"
#include "../../lib/Eigen/Dense"

class EphemerisHistory
{
    std::vector<std::array<double, 6>> _ephemeris;

    std::vector<std::array<double, 6>> _dephemeris;

    std::vector<double> _mjd;

    Ephemeris _current;

    Eigen::Vector3d _position;

public:

    void load(std::string file);

    void set(double mjd);

    const Eigen::Vector3d& get_position() const
    {
        return _position;
    }
};

/**
* @brief orientation in ICRF
*/
class OrientationHistory
{
protected:

    Eigen::Matrix3d _pci2pcef; // ECI 2 ECEF

public:

    virtual void set(double mjd) = 0;

    const Eigen::Matrix3d& get()
    {
        return _pci2pcef;
    }

};

class OrientationHistory_Constant : public virtual OrientationHistory
{
    const Eigen::Matrix3d _axis_ref;

    const double _mjd_ref;

    const double _rotation_rate;

public:

    OrientationHistory_Constant(Eigen::Matrix3d __axis_ref,
        double __mjd_ref,
        double __rotation_rate) :
        _axis_ref(__axis_ref),
        _mjd_ref(__mjd_ref),
        _rotation_rate(__rotation_rate) 
    {
        _pci2pcef = __axis_ref;
    }

    void set(double mjd) override;

};

class OrientationHistory_IERS
{

};

class SolarSystemBody
{
    std::unique_ptr<Gravity> _gravity; // gravity from this body only

    std::unique_ptr<Atmosphere> _atmosphere;

    std::unique_ptr<Geometry> _geometry;

    std::unique_ptr<OrientationHistory> _orientation;

    EphemerisHistory _ephemeris;

    std::vector<SolarSystemBody> _orbiting_bodies;

public:

    enum FRAMENUM
    {
        SOLAR = 0,
        MERCURY = 1,
        VENUS = 2,
        EARTH = 3,
        MARS = 4,
        JUPITER = 5,
        SATURN = 6,
        URANUS = 7,
        NEPTUNE = 8
    };

    enum FRAME_IDENTIFIER
    {
        BARYCENTER = 0,
        MAINBODY = 99
    };

    const SolarSystemBody* const parent;

    const int _id;

    SolarSystemBody(
            std::unique_ptr<Gravity> __gravity,
            std::unique_ptr<Atmosphere> __atmosphere,
            std::unique_ptr<Geometry> __geometry,
            std::unique_ptr<OrientationHistory> __orientation,
            EphemerisHistory __ephemeris,
            int __id,
            const SolarSystemBody* __parent = nullptr) :
        _gravity(std::move(__gravity)),
        _atmosphere(std::move(__atmosphere)),
        _geometry(std::move(__geometry)),
        _orientation(std::move(__orientation)),
        _ephemeris(std::move(__ephemeris)),
        _id(__id),
        parent(__parent) {}

    void set(EpochTime time);

    const Gravity& gravity() { return *_gravity; }

    const Atmosphere& atmosphere() { return *_atmosphere; }

    const Geometry& geometry() { return *_geometry; }

    const OrientationHistory& orientation() { return *_orientation; }

    const EphemerisHistory& ephemeris() { return _ephemeris; }

    const std::vector<SolarSystemBody>& bodies() { return _orbiting_bodies; }


};


class SolarSystem
{

    SolarSystemBody _sol;

    std::vector<SolarSystemBody*> _bodies;

public:

    void add(int identifier);

    void remove(int identifier);

};
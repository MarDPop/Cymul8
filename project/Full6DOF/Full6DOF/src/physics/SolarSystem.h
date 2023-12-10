#pragma once

#include "Time.h"
#include "Gravity.h"
#include "Atmosphere.h"
#include "Ephemeris.h"
#include "Geometry.h"
#include "../util/Table.h"
#include "../../lib/Eigen/Dense"

/**
* @brief provides position relative to parent body, but all in the ICRF (J2000) frame
*/
class EphemerisHistory
{
protected:

    Eigen::Vector3d _position;

public:

    EphemerisHistory() {}

    virtual void load(const std::string& file) {};

    virtual void set(double jd2000) {}

    const Eigen::Vector3d& get_position() const
    {
        return _position;
    }
};

class EphemerisHistory_KeplerElements final : public virtual EphemerisHistory
{
private:
    std::vector<std::array<double, 6>> _ephemeris;

    std::vector<std::array<double, 6>> _dephemeris;

    std::vector<double> _jd2000;

    Ephemeris _current;

    unsigned _tidx = 0;

public:

    EphemerisHistory_KeplerElements() {}

    void load(const std::string& file) override;

    void set(double jd2000) override;
};

/**
* @brief orientation in ICRF
*/
class OrientationHistory
{
protected:

    Eigen::Matrix3d _icrf2fixed; // ECI 2 ECEF

    double _rotation_rate = 0;

public:

    virtual void set(double jd2000) {};

    const Eigen::Matrix3d& get_icrf2fixed() const
    {
        return _icrf2fixed;
    }

    const double get_rotation_rate() const
    {
        return _rotation_rate;
    }
};

class OrientationHistory_Constant final : public virtual OrientationHistory
{
    const Eigen::Matrix3d _axis_ref;

    const double _jd2000_ref;

public:

    OrientationHistory_Constant(Eigen::Matrix3d __axis_ref,
        double __jd2000_ref,
        double __rotation_rate) :
        _axis_ref(__axis_ref),
        _jd2000_ref(__jd2000_ref)
    {
        _icrf2fixed = __axis_ref;
        _rotation_rate = __rotation_rate;
    }

    void set(double jd2000) override;

};

class OrientationHistory_IERS
{
    double xP; // polar x component

    double yP;

    double CIPx;

    double CIPy;

    Eigen::Matrix3d W;

public:

};

class Barycenter;

enum class CLASSIFICATION
{
    NONE = -1,
    PLANET = 0,
    MOON,
    ASTEROID,
    COMET,
    MANMADE
};

class PlanetaryBody
{
    std::unique_ptr<OrientationHistory> _orientation;

    std::unique_ptr<EphemerisHistory> _ephemeris;

    std::unique_ptr<Gravity> _gravity; // gravity from this body only

    std::unique_ptr<Atmosphere> _atmosphere;

    std::unique_ptr<Geometry> _geometry;

    double _TDB2TT;

public:

    const Barycenter& barycenter;

    const int id;

    const CLASSIFICATION classification;

    PlanetaryBody(
        std::unique_ptr<Gravity> __gravity,
        std::unique_ptr<Atmosphere> __atmosphere,
        std::unique_ptr<Geometry> __geometry,
        std::unique_ptr<OrientationHistory> __orientation,
        std::unique_ptr<EphemerisHistory> __ephemeris,
        const Barycenter& __barycenter,
        int __id,
        CLASSIFICATION __classification) :
            _gravity(std::move(__gravity)),
            _atmosphere(std::move(__atmosphere)),
            _geometry(std::move(__geometry)),
            _orientation(std::move(__orientation)),
            _ephemeris(std::move(__ephemeris)),
            barycenter(__barycenter),
            id(__id),
            classification(__classification) {}

    void set_barycentric_dynamical_time(double tdb)
    {
        _ephemeris->set(tdb);
        _orientation->set(tdb + _TDB2TT);
    }

    const OrientationHistory& orientation() { return *_orientation; }

    const EphemerisHistory& ephemeris() { return *_ephemeris; }

    const Gravity& gravity() { return *_gravity; }

    const Atmosphere& atmosphere() { return *_atmosphere; }

    const Geometry& geometry() { return *_geometry; }

    Eigen::Vector3d get_position_in_ICRF(const Eigen::Vector3d& position_around_planet)
    {

    }
};

class Barycenter
{
    std::vector<PlanetaryBody> _bodies;

    std::vector<Barycenter> _barycenters;

    std::unique_ptr<EphemerisHistory> _ephemeris;

public:

    const Barycenter* const parent;

    const double MU;

    const int id;

    Barycenter(const Barycenter* _parent,
        double _MU,
        int _id) : parent(_parent), MU(_MU), id(_id) {}

    void set_ephemeris(std::unique_ptr<EphemerisHistory> ephemeris)
    {
        _ephemeris = std::move(ephemeris);
    }

    void set_barycentric_dynamical_time(double tdb);
};



class SolarSystem
{

    Barycenter _solar_system_barycenter;

public:

    static constexpr double SSB_MU = 1.3289051882026906e+11;

    SolarSystem() : _solar_system_barycenter(nullptr, SSB_MU, 0) {}


};



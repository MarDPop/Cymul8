#pragma once

#include "Time.h"
#include "Gravity.h"
#include "Atmosphere.h"
#include "Ephemeris.h"
#include "Geometry.h"
#include "../util/Table.h"
#include "../../lib/Eigen/Dense"

class EphemerisHistory
{
    std::vector<std::array<double, 6>> _ephemeris;

    std::vector<std::array<double, 6>> _dephemeris;

    std::vector<double> _jd2000;

    Ephemeris _current;

    Eigen::Vector3d _position;

    unsigned _tidx = 0;

public:

    void load(std::string file);

    void set(double jd2000);

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

class OrientationHistory_Constant : public virtual OrientationHistory
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

};

class SolarSystemBody
{
    std::unique_ptr<Gravity> _gravity; // gravity from this body only

    std::unique_ptr<Atmosphere> _atmosphere;

    std::unique_ptr<Geometry> _geometry;

    std::unique_ptr<OrientationHistory> _orientation;

    EphemerisHistory _ephemeris;

    std::vector<SolarSystemBody> _orbiting_bodies;

    double _TDB2TT;

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

    const double MU;

    const SolarSystemBody* const parent;

    const int id;

    SolarSystemBody(
            std::unique_ptr<Gravity> __gravity,
            std::unique_ptr<Atmosphere> __atmosphere,
            std::unique_ptr<Geometry> __geometry,
            std::unique_ptr<OrientationHistory> __orientation,
            EphemerisHistory __ephemeris,
            double __MU,
            int __id,
            const SolarSystemBody* __parent = nullptr) :
        _gravity(std::move(__gravity)),
        _atmosphere(std::move(__atmosphere)),
        _geometry(std::move(__geometry)),
        _orientation(std::move(__orientation)),
        _ephemeris(std::move(__ephemeris)),
        MU(__MU),
        id(__id),
        parent(__parent) {}

    void set_barycentric_dynamical_time(double tdb)
    {
        _ephemeris.set(tdb);
        _orientation->set(tdb + _TDB2TT);
    }

    const Gravity& gravity() { return *_gravity; }

    const Atmosphere& atmosphere() { return *_atmosphere; }

    const Geometry& geometry() { return *_geometry; }

    const OrientationHistory& orientation() { return *_orientation; }

    const EphemerisHistory& ephemeris() { return _ephemeris; }

    const std::vector<SolarSystemBody>& bodies() { return _orbiting_bodies; }
};

class RadiationPressure
{

public:

    virtual double get(double R, double T) { return 0.0; }
};

class Planetary_System
{
public:
    virtual void set_TDB_reference(double jd2000);

    virtual void set_utc_jd2000(double jd2000);
};

class SolarSystem final : public virtual Planetary_System
{

    SolarSystemBody _sol;

    std::unique_ptr<RadiationPressure> _radiation_pressure;

    std::vector<SolarSystemBody*> _bodies;

    double UTC2TDB;

public:

    void add(int identifier);

    void remove(int identifier);

    void set_TDB_reference(double jd2000) override;

    void set_utc_jd2000(double jd2000) override
    {
        double tdb = jd2000 + UTC2TDB; // 
        for (auto* body : _bodies)
        {
            body->set_barycentric_dynamical_time(tdb);
        }
    }

};

enum MODEL_FIDELITY
{
    LOW,
    MEDIUM,
    HIGH
};

template<class G_EARTH, class G_MOON, class ATM_EARTH>
class Earth_Moon_Sun final : public virtual Planetary_System
{
    static_assert(std::is_base_of<Gravity, G_EARTH>::value, "G not derived from Guidance");
    static_assert(std::is_base_of<Gravity, G_MOON>::value, "G not derived from Guidance");
    static_assert(std::is_base_of<Atmosphere, ATM_EARTH>::value, "G not derived from Guidance");

public:

    struct Earth
    {
        G_EARTH gravity;

        ATM_EARTH atmosphere;

        BasicTable<double, 9> ECI2ECF;

        static constexpr double MU = 3.986004418e5; // km/ kg s3
    };

    struct Moon
    {
        BasicTable<double, 3> pos_ECI;

        BasicTable<double, 9> ECI2LCF;

        G_MOON gravity;

        static constexpr double MU = 4.9048695e3; // km/ kg s3
    };

    struct Sun
    {
        BasicTable<double, 3> pos_ECI;

        static constexpr double MU = 1.32712440018e11; // km/ kg s3
    };

    void load(const std::string& ephemeris_earth_icrf,
        const std::string& ephemeris_moon_icrf,
        const std::string& ephemeris_sun_icrf,
        const std::string& moon_orientation,
        const std::string& finals_file);

    void load_eci(const std::string& ephemeris_moon_eci,
        const std::string& ephemeris_sun_eci,
        const std::string& finals_file);

    void set(double time);

    void get_gravity()

private:

    double _JD_reference;

    Earth _earth;

    Moon _moon;

    Sun _sun;

};
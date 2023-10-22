#pragma once

#include <memory>
#include "Coordinates.h"

#include "Gravity.h"
#include "Atmosphere.h"
#include "Geometry.h"

class Environment
{
public:

    enum STATE
    {
        LOCAL = 0,
        FLYING_HIGH,
        LEO,
        HEO,
        CISLUNAR
    };

private:
    const Vehicle& _vehicle;

    Coordinate::ECI _ECI;

    Coordinate::ECEF _ECEF;

    Coordinate::Geodetic _LLA;

    Coordinate::LocalTangentPlane _LTP;

    Coordinate::Spherical _RTP;

    EpochTime _launch_time;

    Coordinate::Geodetic _launch_lla;

    STATE _current_state;

    std::unique_ptr<Gravity> gravity;

    std::unique_ptr<Atmosphere> atmosphere;

    std::unique_ptr<Geometry> geometry;

public:

    Environment(const Vehicle& vehicle) :
        _vehicle(vehicle) {}
    
    void set_launch(EpochTime launch_time,
                    Coordinate::Geodetic launch)
    {
        _launch_time = launch_time;
        _launch_lla = launch;
    }

    void update(double T);

};
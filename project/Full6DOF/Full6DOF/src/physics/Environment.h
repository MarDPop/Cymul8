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
        NONE = -1,
        LAUNCH_LANDING = 0,
        ATMOSPHERIC,
        LOW_ORBIT,
        HIGH_ORBIT,
        INTERPLANETARY
    };

private:
    const Vehicle& _vehicle;

    Planet _current_planet;

    Coordinate::ECI _ECI;

    Coordinate::ECEF _ECEF;

    Coordinate::Geodetic _LLA;

    Coordinate::LocalTangentPlane _LTP;

    Coordinate::Spherical _RTP;

    EpochTime _launch_time;

    Coordinate::Geodetic _launch_lla;

    STATE _current_state;

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
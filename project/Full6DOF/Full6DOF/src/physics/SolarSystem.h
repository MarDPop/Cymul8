#pragma once

#include "Time.h"
#include "Gravity.h"
#include "Atmosphere.h"
#include "Ephemeris.h"
#include "Geometry.h"

struct SolarSystemBody
{
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

    EphemerisHistory ephemeris;

    Gravity& gravity;

    Atmosphere& atmosphere;

    Geometry& geometry;

    std::vector<SolarSystemBody> orbiting_bodies;

    SolarSystemBody(Gravity& __gravity,
           Atmosphere& __atmosphere,
           EphemerisHistory __ephemeris,
           Geometry& __geometry) :
        gravity(__gravity),
        atmosphere(__atmosphere),
        ephemeris(std::move(__ephemeris)),
        geometry(__geometry) {}

};
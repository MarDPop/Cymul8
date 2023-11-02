#pragma once

#include <memory>

#include "../../lib/tinyxml/tinyxml2.h"

#include "../vehicle/Vehicle.h"
#include "../ode/ODE_standalone.h"
#include "../physics/SolarSystem.h"

class SimulationConfiguration
{
    double DISTANCE_NOT_LOCAL = 30;//km
};

template<class V>
class Simulation
{
    SimulationConfiguration _config;
    
    V _vehicle;

    ode<V, double> _ode;

    SolarSystem solar_system;

public:

    Simulation(const tinyxml2::XMLDocument& config);
    virtual ~Simulation() {}

};
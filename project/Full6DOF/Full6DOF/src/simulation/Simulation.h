#pragma once

#include <memory>

#include "../../lib/tinyxml/tinyxml2.h"

#include "../vehicle/Vehicle.h"
#include "../ode/ODE.h"
#include "../physics/SolarSystem.h"

class SimulationConfiguration
{
    double DISTANCE_NOT_LOCAL = 30;//km
};

class Simulation
{
    SimulationConfiguration _config;
    
    std::unique_ptr<Vehicle> _vehicle;

    std::unique_ptr<ODE> _ode;

    SolarSystem solar_system;

public:

    Simulation(const tinyxml2::XMLDocument& config);
    virtual ~Simulation() {}

};
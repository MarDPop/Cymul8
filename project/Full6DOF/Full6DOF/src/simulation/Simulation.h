#pragma once

#include <memory>

#include "../../lib/tinyxml/tinyxml2.h"

#include "../physics/Stage.h"
#include "../vehicle/Vehicle.h"
#include "../ode/ODE.h"

class SimulationConfiguration
{
    double ALTITUDE_RARIFIED_DRAG = 100; // km
};

class Simulation
{
    SimulationConfiguration _config;
    
    std::unique_ptr<Vehicle> _vehicle;

    std::unique_ptr<ODE> _ode;

public:

    Simulation(const tinyxml2::XMLDocument& config);
    virtual ~Simulation() {}

};
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

enum class SIMULATION_STATE
{
    NONE = -1,
    LAUNCH_LANDING = 0,
    ATMOSPHERIC,
    LOW_ORBIT,
    HIGH_ORBIT,
    INTERPLANETARY,
    COASTING
};

template<class V>
class Simulation
{
    SimulationConfiguration _config;
    
    V _vehicle;

    ode<V, double> _ode;

    FrameReference _local_frame;

    SolarSystem _solar_system;

    SIMULATION_STATE _state;

public:

    Simulation(const tinyxml2::XMLDocument& config);
    virtual ~Simulation() {}



};
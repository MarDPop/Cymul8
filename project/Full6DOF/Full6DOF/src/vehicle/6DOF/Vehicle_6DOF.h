#pragma once

#include <memory>

#include "../Vehicle.h"
#include "../Action.h"
#include "Aerodynamics.h"
#include "Propulsion.h"

class Vehicle_6DOF : public virtual Vehicle
{
    std::unique_ptr<Aerodynamics> _aero;

    std::unique_ptr<Propulsion> _propulsion;

    BodyAction _action_sum;

public:

    Vehicle_6DOF(Environment& environment) : Vehicle(environment) {}
};
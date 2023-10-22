#pragma once

#include "../Vehicle.h"

#include <vector>

class Vehicle_3DOF : public virtual Vehicle
{
    //Thruster _thruster;

    //Aerodynamics _aero;

    Eigen::Vector3d _sum_forces;

public:

    Vehicle_3DOF(Environment& environment) : Vehicle(environment) {}

};
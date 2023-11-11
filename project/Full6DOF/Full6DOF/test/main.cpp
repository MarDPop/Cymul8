#include "../src/ode/ODE.h"

#include <cmath>
#include <string>
#include <fstream>
#include <iostream>

#include "../lib/Eigen/Dense"

#include "ODE_standalone.h"
#include "../src/vehicle/Vehicle.h"

void test2()
{
    
    Vehicle_3DOF_Standard rocket;

    ode<Vehicle_3DOF_Standard, double> rocket_ode(rocket);


}

int main(int argc, char** argv) 
{
    test2();

    return 0;
}
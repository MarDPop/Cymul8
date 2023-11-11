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
    ode<Vehicle_3DOF_Standard, double> rocket;
    // Vehicle_3DOF_Standard rocket;
    // std::cout << rocket.get_num_states() << std::endl;
}

int main(int argc, char** argv) 
{
    test2();

    return 0;
}
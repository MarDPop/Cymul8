#pragma once

#include <array>
#include <algorithm>
#include "../util/Vector.h"

/**
 * @brief a class which defines behavior for dynamics which can be implemented by an ODE
 * 
 * @tparam NSTATES 
 */
template<unsigned NSTATES>
struct Fixed_Size_Dynamics 
{
    /**
     * @brief this method describes the rate of change to a state given it's current state and time
     * 
     * @param x current state
     * @param time current time
     * @param dx rate of change of state
     * @return true if state is valid
     * @return false if state is invalid and ODE should stop
     */
    virtual bool set_state(const std::array<double, NSTATES>& x,
        const double& time, 
        std::array<double, NSTATES>& dx) = 0;
};

/**
 * @brief a class which defines behavior for dynamics which can be implemented by an ODE
 * 
 * 
 */
class Dynamics 
{

public:
    /**
     * @brief number of states
     * 
     */
    const unsigned NSTATES;

    /**
     * @brief this method describes the rate of change to a state given it's current state and time
     * 
     * @param x current state
     * @param time current time
     * @param dx rate of change of state
     * @return true if state is valid
     * @return false if state is invalid and ODE should stop
     */
    virtual bool set_state(const double* x, 
        const double& time, 
        double* dx) = 0;

    Dynamics(unsigned _NSTATES) : 
        NSTATES(_NSTATES) {}
};
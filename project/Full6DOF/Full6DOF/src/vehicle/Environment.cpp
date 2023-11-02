#include "Environment.h"

void Environment::update_launch(const double* state, double talo)
{
    _TALO = talo;
}

void Environment::update_atm(const double* state, double talo)
{
    _TALO = talo;
}

void Environment::update_low_orbit(const double* state, double talo)
{

}

void Environment::update_high_orbit(const double* state, double talo)
{

}

void Environment::update_interplanetary(const double* state, double talo)
{

}

void Environment::update_coast(const double* state, double talo)
{

}


void Environment::set_state(STATE state)
{
    _current_state = state;
    switch (state)
    {
    case LAUNCH_LANDING:
        update = &Environment::update_launch;
        break;
    case ATMOSPHERIC:
        update = &Environment::update_atm;
        break;
    case LOW_ORBIT:
        update = &Environment::update_low_orbit;
        break;
    case HIGH_ORBIT:
        update = &Environment::update_high_orbit;
        break;
    case INTERPLANETARY:
        update = &Environment::update_interplanetary;
        break;
    default:
        update = &Environment::update_coast;
        break;
    }
}
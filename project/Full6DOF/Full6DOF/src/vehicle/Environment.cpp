#include "Environment.h"

void Environment::update_launch_landing(const Eigen::Vector3d& position,
    const Eigen::Vector3d& velocity,
    double talo)
{
    _TALO = talo;
}

void Environment::update_atm(const Eigen::Vector3d& position,
    const Eigen::Vector3d& velocity,
    double talo)
{
    _TALO = talo;
}

void Environment::update_low_orbit(const Eigen::Vector3d& position,
    const Eigen::Vector3d& velocity,
    double talo)
{

}

void Environment::update_high_orbit(const Eigen::Vector3d& position,
    const Eigen::Vector3d& velocity,
    double talo)
{

}

void Environment::update_interplanetary(const Eigen::Vector3d& position,
    const Eigen::Vector3d& velocity,
    double talo)
{

}

void Environment::update_coast(const Eigen::Vector3d& position,
    const Eigen::Vector3d& velocity,
    double talo)
{

}


void Environment::set_state(STATE state)
{
    _current_state = state;
    switch (state)
    {
    case STATE::LAUNCH_LANDING:
        update = &Environment::update_launch_landing;
        break;
    case STATE::ATMOSPHERIC:
        update = &Environment::update_atm;
        break;
    case STATE::LOW_ORBIT:
        update = &Environment::update_low_orbit;
        break;
    case STATE::HIGH_ORBIT:
        update = &Environment::update_high_orbit;
        break;
    case STATE::INTERPLANETARY:
        update = &Environment::update_interplanetary;
        break;
    default:
        update = &Environment::update_coast;
        break;
    }
}
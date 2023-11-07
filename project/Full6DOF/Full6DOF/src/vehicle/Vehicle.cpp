#include "Vehicle.h"


void Vehicle_3DOF_Standard::update_accelerations(double t) 
{
    Eigen::Vector3d force;
    force.setZero();

    if (_thruster->is_active())
    {
        _thruster->update(_environment.get_air(), _environment.get_aero_data(), t);
        force += _orientation.col(0) * _thruster->get_thrust();
    }

    if (_environment.in_air())
    {
        _aerodynamics->update(_environment.get_air(), _environment.get_aero_data(), t);
        force += _aerodynamics->get_action();
    }

    _acceleration += force * (1.0 / _state.mass);
}


void Vehicle_6DOF_Standard::update_accelerations(double t)
{
    double dm = _state.mass - _inertia_empty.mass;
    _center_of_mass = _inertia_empty.center_of_mass + _center_of_mass_delta * dm;
    for (auto i = 0u; i < MOMENT_CONSTANTS::FULL; i++)
    {
        _moment_of_inertia.I[i] = _inertia_empty.moment_of_inertia.I[i] + _moment_of_inertia_delta.I[i] * dm;
    }

    if (_propulsion->get_thruster().is_active())
    {
        _propulsion->update_thrust(_environment.get_air(), _environment.get_aero_data(), t);
    }

    if (_environment.in_air())
    {
        _aerodynamics->update(_environment.get_air(), _environment.get_aero_data(), t);
    }

    BodyAction totalAction = _propulsion->get_action() + _aerodynamics->get_action(); // at zero

    _acceleration += totalAction.force * (1.0 / _state.mass);

    _moment_of_inertia.get_angular_acceleration_body(_state.angular_velocity,
        totalAction.get_torque(_center_of_mass),
        _angular_acceleration.data());

    _mass_rate = -_propulsion->get_thruster().get_mass_rate();
}

void Vehicle_6DOF_Standard::set_inertia(const Inertia<MOMENT_CONSTANTS::FULL>& empty, 
                                        const Inertia<MOMENT_CONSTANTS::FULL>& full)
{
    _inertia_empty = empty;
    double dm = full.mass - empty.mass;
    _center_of_mass_delta = (full.center_of_mass - empty.center_of_mass) / dm;
    for (auto i = 0u; i < MOMENT_CONSTANTS::FULL; i++)
    {
        _moment_of_inertia_delta.I[i] = (full.moment_of_inertia.I[i] - empty.moment_of_inertia.I[i]) / dm;
    }
}

void Vehicle_Components::update_accelerations(double t)
{

}
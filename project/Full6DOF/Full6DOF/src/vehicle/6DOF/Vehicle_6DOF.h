#pragma once

#include <memory>

#include "../../physics/Body.h"
#include "../Vehicle.h"
#include "../Action.h"
#include "Aerodynamics.h"
#include "Propulsion.h"

template<MOMENT_CONSTANTS NDEG>
class Vehicle_6DOF : public virtual Vehicle
{
    Eigen::Quaterniond _orientation;

    Eigen::Vector3d _angular_velocity;

    Eigen::Vector3d _angular_acceleration;

    Inertia<NDEG> _inertia;

    Inertia<NDEG> _inertia_rate;

    std::unique_ptr<Aerodynamics> _aerodynamics;

    std::unique_ptr<Propulsion> _propulsion;

    BodyAction _action_sum;

public:

    unsigned get_num_states() const override
    {
        return 17 + DEG + gnc.get_number_control_states();
    }

    void operator()(const double* x, const double time, double* dx) override
    {
        memcpy(_position.data(), x, 3*sizeof(double));
        memcpy(_velocity.data(), x + 3, 3*sizeof(double));
        memcpy(_velocity.data(), x + 3, 3 * sizeof(double));
    }

    void set_aerodynamics(std::unique_ptr<Aerodynamics> aero)
    {
        _aero = std::move(aero);
    }

    void set_propulsion(std::unique_ptr<Aerodynamics> propulsion)
    {
        _propulsion = std::move(propulsion);
    }

    const Eigen::Quaterniond& get_orientation() const { return _orientation; }

    const Eigen::Vector3d& get_angular_velocity() const { return _angular_velocity; }

    const Eigen::Vector3d& get_angular_acceleration() const { return _angular_acceleration; }

    const Inertia<NDEG>& get_inertia() const { return _inertia; }

    const Inertia<NDEG>& get_inertia_rate() const { return _inertia_rate; }

    const Aerodynamics& get_aerodynamics() const { return *_aerodynamics; }

    const Propulsion& get_propulsion() const { return *_propulsion; }

    const BodyAction& get_action_sum() const { return _action_sum; }
};
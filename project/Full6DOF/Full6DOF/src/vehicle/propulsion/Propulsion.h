#pragma once

#include "Thruster.h"
#include "../Environment.h"
#include "../../physics/Action.h"
#include <memory>

template<class T>
class Propulsion_T
{
protected:

    BodyAction _action;

    Eigen::Vector3d _thrust_vector; // in body

    T _thruster;

public:

    Propulsion_T() : _thrust_vector{ 1, 0, 0 }
    {
        _action.zero();
    }

    void update_thrust(const Air& air, const AeroData& aero_data, double t)
    {
        _thruster.update(air, aero_data, t);
        _action.force = _thrust_vector*_thruster.get_thrust();
    }

    void set_thrust_vector(const Eigen::Vector3d& __thrust_vector)
    {
        _thrust_vector = __thrust_vector;
    }

    const T& get_thruster() const
    {
        return _thruster;
    }

    const BodyAction get_action() const
    {
        return _action;
    }
};

class Propulsion 
{
protected:

    BodyAction _action;

    Eigen::Vector3d _thrust_vector; // in body

    std::unique_ptr<Thruster> _thruster;

public:

    Propulsion() : _thrust_vector{ 1, 0, 0 }
    {
        _action.setZero();
    }

    void set_thruster(std::unique_ptr<Thruster> __thruster)
    {
        _thruster = std::move(__thruster);
    }

    void update_thrust(const Air& air, 
        const AeroData& aero_data, 
        double t)
    {
        _thruster->update(air, aero_data, t);
        _action.force = _thrust_vector * _thruster->get_thrust();
    }

    void set_thrust_vector(const Eigen::Vector3d& __thrust_vector)
    {
        _thrust_vector = __thrust_vector;
    }

    const Thruster& get_thruster() const
    {
        return *_thruster;
    }

    const BodyAction get_action() const
    {
        return _action;
    }
};
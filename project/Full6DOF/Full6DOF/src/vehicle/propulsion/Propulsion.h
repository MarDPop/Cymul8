#pragma once

#include "Thruster.h"
#include "../Environment.h"
#include "../../physics/Action.h"

template<class T, class F>
class Propulsion : public virtual BodyAction
{
protected:

    Eigen::Vector3d _thrust_vector; // in body

    T _thruster;

    F _tank;

public:

    Propulsion() : _thrust_vector{ 1, 0, 0 }
    {
        zero();
    }

    void update_thrust(const Air& air, const AeroData& aero_data, double t)
    {
        _thruster.update(air, aero_data, t);
        force = _thrust_vector*_thruster.get_thrust();
    }

    void update_inertia(double m, double t, const Eigen::Vector3d& acceleration_body)
    {
        _tank.update(m, t, acceleration_body);
    }

    void set_thrust_vector(const Eigen::Vector3d& __thrust_vector)
    {
        _thrust_vector = __thrust_vector;
    }

    const T& get_thruster() const
    {
        return _thruster;
    }

    const F& get_tank() const
    {
        return _tank;
    }
};
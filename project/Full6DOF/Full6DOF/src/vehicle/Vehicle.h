#pragma once

#include "GNC.h"
#include "Environment.h"

#include <memory>


template<class G, class N, class C>
class Vehicle
{
    
protected:

    Eigen::Vector3d _position;

    Eigen::Vector3d _velocity;

    double _mass;

    Eigen::Vector3d _acceleration;

    double _mass_rate;

public:

    Environment environment;

    GuidanceNavigationControl<G,N,C> gnc;

    virtual void operator()(const double* x, const double t, double* dx) = 0;

    virtual unsigned get_num_states() const = 0;

    virtual void get_state(double* state)
    {
        memcpy(state, _position.data(), 3 * sizeof(double));
        memcpy(state + 3, _velocity.data(), 3 * sizeof(double));
        state[6] = _mass;
    }

    const Eigen::Vector3d& get_position() const
    {
        return _position;
    }

    const Eigen::Vector3d& get_velocity() const
    {
        return _velocity;
    }

    const Eigen::Vector3d& get_accelaration() const
    {
        return _acceleration;
    }

    double get_mass() const
    {
        return _mass;
    }

    double get_mass_rate() const
    {
        return _mass_rate;
    }
};


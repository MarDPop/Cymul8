#pragma once

#include "GNC.h"
#include "Environment.h"

#include <memory>

class Vehicle
{
protected:

    Eigen::Vector3d _position;

    Eigen::Vector3d _velocity;

    Eigen::Vector3d _acceleration;

    double _mass;

    double _mass_rate;

public:

    Environment environment;

    std::unique_ptr<GNC> gnc;

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


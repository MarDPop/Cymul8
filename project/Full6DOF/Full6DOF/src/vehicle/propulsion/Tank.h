#pragma once

#include "../../physics/Body.h"

template<MOMENT_CONSTANTS NDEG>
class Tank
{
    Eigen::Quaterniond _orientation;

    Inertia<NDEG> _inertia;

public:

    virtual void update(double m, double t, const Eigen::Vector3d& acceleration_body) = 0;

    const Inertia<NDEG>& get_inertia() const
    {
        return _inertia;
    }

};
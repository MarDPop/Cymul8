#pragma once

#include "../../physics/Body.h"

template<MOMENT_CONSTANTS NDEG>
class Tank
{
    Eigen::Quaterniond _orientation;

    Inertia<NDEG> _inertia;

public:

    const Inertia<NDEG>& get_inertia() const
    {
        return _inertia;
    }

    virtual void update(double m, double t) = 0;

    virtual Inertia<NDEG> get_inertia_rate(double mass_rate, const Eigen::Vector3d& acceleration_body) const = 0;

};
#pragma once

#include "../Vehicle.h"
#include "Propulsion.h"
#include "Aerodynamics.h"


template<class GNC>
class Vehicle_3DOF : public virtual Vehicle<BODY_POINT_MASS<double>, GNC>
{
protected:

    std::unique_ptr<Propulsion> _propulsion;

    std::unique_ptr<Aerodynamics> _aero;

    Eigen::Vector3d _sum_forces;

public:

    void set_aerodynamics(std::unique_ptr<Aerodynamics> aero)
    {
        _aero = std::move(aero);
    }

    void set_propulsion(std::unique_ptr<Aerodynamics> propulsion)
    {
        _propulsion = std::move(propulsion);
    }

    void operator()(const Float* x, const Float t, Float* dx) override
    {
        memcpy(_body.state.data(), x, sizeof(_body.state));
    }
};

template<class T, class A, class GNC>
class Vehicle_3DOF_T : public virtual Vehicle<BODY_POINT_MASS<double>, GNC>
{
    static_assert(std::is_base_of<Aerodynamics>, A > ::value, "A not derived from Aerodynamics");
    static_assert(std::is_base_of<Propulsion>, P > ::value, "P not derived from Propulsion");

protected:

    T _thruster;

    A _aero;

    Eigen::Vector3d _sum_forces;

public:

    void operator()(const Float* x, const Float t, Float* dx) override
    {
        memcpy(_body.state.data(), x, sizeof(_body.state));
    }
};
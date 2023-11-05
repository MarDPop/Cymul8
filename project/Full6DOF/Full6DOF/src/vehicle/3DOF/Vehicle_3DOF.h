#pragma once

#include "../Vehicle.h"
#include "Aero.h"
#include "Thruster.h"


template<class GNC>
class Vehicle_3DOF : public virtual Vehicle<BODY_POINT_MASS<double>, GNC>
{
protected:

    std::unique_ptr<Thruster> _thruster;

    std::unique_ptr<Aero> _aero;

public:

    void set_aerodynamics(std::unique_ptr<Aero> aero)
    {
        _aero = std::move(aero);
    }

    void set_thruster(std::unique_ptr<Thruster> thruster)
    {
        _thruster = std::move(thruster);
    }

    void operator()(const Float* x, const Float t, Float* dx) override
    {
        memcpy(_body.state.data(), x, sizeof(_body.state));
    }
};

template<class T, class A, class GNC>
class Vehicle_3DOF_T : public virtual Vehicle<BODY_POINT_MASS<double>, GNC>
{
    static_assert(std::is_base_of<Thruster>, T > ::value, "P not derived from Propulsion");
    static_assert(std::is_base_of<Aero>, A > ::value, "A not derived from Aerodynamics");

protected:

    T _thruster;

    A _aero;

public:

    void operator()(const Float* x, const Float t, Float* dx) override
    {
        memcpy(_body.state.data(), x, sizeof(_body.state));
    }
};
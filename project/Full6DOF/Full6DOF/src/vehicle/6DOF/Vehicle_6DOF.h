#pragma once

#include <memory>

#include "../../physics/Body.h"
#include "../../physics/Action.h"
#include "../Vehicle.h"
#include "Aerodynamics.h"
#include "Propulsion.h"

template<MOMENT_CONSTANTS NDEG, class GNC>
class Vehicle_6DOF : public virtual Vehicle<BODY<NDEG>, GNC>
{
protected:

    std::unique_ptr<Aerodynamics> _aerodynamics;

    std::unique_ptr<Propulsion> _propulsion;

public:

    void set_aerodynamics(std::unique_ptr<Aerodynamics> aero)
    {
        _aero = std::move(aero);
    }

    void set_propulsion(std::unique_ptr<Aerodynamics> propulsion)
    {
        _propulsion = std::move(propulsion);
    }

    void operator()(const double* x, const double t, double* dx)
    {
        memcpy(_body.state.data(), x, sizeof(_body.state));
    }
};

template<class A, class P, MOMENT_CONSTANTS NDEG, class GNC>
class Vehicle_6DOF_T : public virtual Vehicle<BODY<NDEG>, GNC>
{
    static_assert(std::is_base_of<Aerodynamics>, A > ::value, "A not derived from Aerodynamics");
    static_assert(std::is_base_of<Propulsion>, P > ::value, "P not derived from Propulsion");

protected:

    A _aerodynamics;

    P _propulsion;

public:

    void operator()(const double* x, const double t, double* dx)
    {
        memcpy(_body.state.data(), x, sizeof(_body.state));
        _environment.update(_body.position, _body.velocity, t);

        _aerodynamics.update(_environment, t);
        _propulsion.update(_environment, t);

        BodyAction sumActions = _aerodynamics.get_action() + 
                                _propulsion.get_action();


    }
};
#pragma once

#include "../Environment.h"
#include "../../physics/Action.h"

template<class T>
class Aerodynamics
{
    T _action;

public:

    virtual void update(const Air& air, const AeroData& aeroData, double t) = 0;

    const T& get_action() const
    {
        return _action;
    }
};


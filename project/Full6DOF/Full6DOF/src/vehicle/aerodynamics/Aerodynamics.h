#pragma once

#include "../Environment.h"
#include "../../physics/Action.h"

template<class T>
class Aerodynamics
{
protected:

    T _action;

public:

    Aerodynamics()
    {
        _action.setZero();
    }

    virtual void update(const AeroEnvironment& aeroData, 
        double t) {}

    const T& get_action() const
    {
        return _action;
    }
};


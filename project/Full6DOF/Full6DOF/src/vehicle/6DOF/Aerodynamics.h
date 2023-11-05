#pragma once

#include "../Environment.h"
#include "../../physics/Action.h"

class Aerodynamics
{
protected:

    BodyAction _action;

public:

    const BodyAction& get_action()
    {
        return _action;
    }

    virtual void update(const Environment& environment, double time) = 0;
};
#pragma once

#include "../Environment.h"

class Thruster
{

    double _thrust;

    double _mass_rate;

    bool _active;

public: 

    virtual void update(const Air& air, const AeroData& aero_data, double t) = 0;

    virtual void start()
    {
        _active = true;
    }

    virtual void shutdown()
    {
        _active = false;
    }

    bool is_active() const
    {
        return _active;
    }

    double get_thrust() const
    {
        return _thrust;
    }

    double get_mass_rate() const
    {
        return _mass_rate;
    }
};
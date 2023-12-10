#pragma once

#include "../Environment.h"

class Thruster
{
protected:

    double _thrust = 0;

    double _mass_rate = 0;

    bool _active = false;

public: 

    Thruster() {}

    Thruster(double thrust, double mass_rate) :
        _thrust(thrust),
        _mass_rate(mass_rate)
    {}

    virtual ~Thruster() {}

    virtual void update(const AeroEnvironment* aero, double t) {}

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
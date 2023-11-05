#pragma once

#include "GNC.h"
#include "Environment.h"
#include "../physics/Body.h"

template<class B, class G>
class Vehicle
{
    static_assert(std::is_base_of<GNC<B>>, G>::value, "G not derived from Guidance");

protected:

    friend class Vehicle_ODE;

    B _body;

    G _gnc;

    Environment _environment;

public:

    virtual void operator()(const double* x, const double t, double* dx) = 0;

    static unsigned get_num_states() const
    {
        B::N_STATES + gnc.get_control_states();
    }

    const B& get_body() const
    {
        return _body;
    }

    const G& get_GNC() const
    {
        return _gnc;
    }

    const Environment& get_environment() const
    {
        return _environment;
    }
};

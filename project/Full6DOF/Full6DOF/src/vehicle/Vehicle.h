#pragma once

#include "GNC.h"
#include "../physics/Environment.h"

#include <memory>

class Vehicle
{
protected:
    friend class Environment;

    std::unique_ptr<GNC> _gnc;

    Environment _environment;

public:

    void set_gnc(std::unique_ptr<GNC> gnc)
    {
        _gnc = std::move(gnc);
    }

    GNC& get_gnc()
    {
        return *_gnc;
    }

    const GNC& get_gnc() const
    {
        return *_gnc;
    }
};


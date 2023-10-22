#pragma once

#include "../../lib/Eigen/Dense"

#include "Time.h"
#include "Coordinates.h"

#include <vector>
#include <memory>

/**
 * @brief 
 * 
 */
struct Gravity : public virtual Eigen::Vector3d
{

    virtual void compute(   const double* position,
                            const double R,
                            const double T) = 0;
};

class Multi_Gravity : public virtual Gravity
{

    std::vector<std::unique_ptr<Gravity>> _gravities;

public:

    Multi_Gravity() {}

    void add_gravity(std::unique_ptr<Gravity>& gravity)
    {
        this->_gravities.push_back(std::move(gravity));
    }

    
};
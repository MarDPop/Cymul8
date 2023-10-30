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
class Gravity
{

public:

    virtual void set_time(EpochTime time) = 0;

    virtual void compute(   const Coordinate::GeocentricFixed& position,
                            const double R,
                            Eigen::Vector3d& acceleration) = 0;
};
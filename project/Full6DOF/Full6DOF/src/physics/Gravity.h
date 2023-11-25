#pragma once

#include "../../lib/Eigen/Dense"

#include "EpochTime.h"
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

    virtual void compute(   const Coordinate::GeocentricFixed& position,
                            const double R,
                            Eigen::Vector3d& acceleration) const = 0;

    //remember all position and velocity referenced in body

    static Eigen::Vector3d fictional_forces_Z_rotation(const Eigen::Vector3d& position,
        const Eigen::Vector3d& velocity,
        const double& rotation);

    static Eigen::Vector3d fictional_forces_const_rotation(const Eigen::Vector3d& position,
        const Eigen::Vector3d& velocity,
        const Eigen::Vector3d& rotation);

    static Eigen::Vector3d fictional_forces(const Eigen::Vector3d& position,
        const Eigen::Vector3d& velocity,
        const Eigen::Vector3d& rotation,
        const Eigen::Vector3d& rotation_rate);
};

class Gravity_Local : virtual Gravity
{

    double _g0;

    double _R0;

    Eigen::Vector3d _earth_axis_orientation_in_frame;

public:

    Gravity_Local(double g0,
        double R0,
        Eigen::Vector3d earth_axis_orientation_in_frame) :
        _g0(-g0),
        _R0(R0),
        _earth_axis_orientation_in_frame(earth_axis_orientation_in_frame) {}

    void compute(const Coordinate::GeocentricFixed& position,
        const double R,
        Eigen::Vector3d& acceleration) const override;

};
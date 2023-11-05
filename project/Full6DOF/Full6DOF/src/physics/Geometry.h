#pragma once

#include "../../lib/Eigen/Dense"
#include "Coordinates.h"

namespace WGS84
{

}

class Geometry
{

public:

    virtual bool inside(const Eigen::Vector3d& point,
                        const double R) const
    {
        return false;
    }

    virtual bool intersects(const Eigen::Vector3d& start,
                            const Eigen::Vector3d& LOS) const
    {
        return false;
    }

    virtual Eigen::Matrix3d LTP(const Eigen::Vector3d& point);
};

class Geometry_Sphere : public virtual Geometry
{

    const double _R;

public:

    Geometry_Sphere(double __R) :
        _R(__R) {}

    bool inside(const Eigen::Vector3d& point,
        const double R) const override
    {
        return R > _R;
    }

    bool intersects(const Eigen::Vector3d& start,
        const Eigen::Vector3d& LOS) const
    {
        return false;
    }
};

class Geometry_Ellipsoid : public virtual Geometry
{
    const double _polar_radius;

    const double _equatorial_radius;

public:

    Geometry_Ellipsoid(double __A,
                    double __B) :
        _polar_radius(__A),
        _equatorial_radius(__B) {}

    bool inside(const Eigen::Vector3d& point,
        const double R) const override
    {
        return R > _polar_radius;
    }

    bool intersects(const Eigen::Vector3d& start,
        const Eigen::Vector3d& LOS) const
    {
        return false;
    }
};
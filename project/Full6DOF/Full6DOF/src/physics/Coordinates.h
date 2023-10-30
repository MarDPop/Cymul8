#pragma once

#include "../../lib/Eigen/Dense"

namespace Coordinate
{
    struct Frame
    {
        Frame* reference_frame = nullptr;

        Eigen::Matrix3d orientation;

        Eigen::Vector3d location;

        Frame() : orientation(Eigen::Matrix3d::Identity()),
            location(Eigen::Vector3d::Zero()) {}

        Frame(const Eigen::Matrix3d& _orientation, const Eigen::Vector3d& _location) : 
            orientation(_orientation), location(_location) {}
    };

    struct FrameWrapper
    {
        Frame* reference_frame = nullptr;

        Eigen::Matrix3d& orientation;

        Eigen::Vector3d& location;

        FrameWrapper(Eigen::Matrix3d& _orientation, Eigen::Vector3d& _location) :
            orientation(_orientation), location(_location) {}
    };

    typedef Eigen::Vector3d _3D;

    class GeocentricFixed : public virtual Eigen::Vector3d
    {
    public:

    };

    class GeocentricInertial : public virtual Eigen::Vector3d
    {
    public:

    };

    class Geodetic : public virtual Eigen::Vector3d
    {
    public:

    };

    class ENU : public virtual Eigen::Vector3d
    {
    public:

    };

    class Spherical : public virtual Eigen::Vector3d
    {
    public:

    };

}

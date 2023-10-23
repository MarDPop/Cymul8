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

    class ECEF : public virtual _3D
    {
    public:

    };

    class ECI : public virtual _3D
    {
    public:

    };

    class Geodetic : public virtual _3D
    {
    public:

    };

    class LocalTangentPlane : public virtual _3D
    {
    public:

    };

    class Spherical : public virtual _3D
    {
    public:

    };

}

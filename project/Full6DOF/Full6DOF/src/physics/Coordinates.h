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

    union GeocentricFixed
    {
        Eigen::Vector3d vec;
        struct 
        {
            double x;
            double y;
            double z;
        };

        GeocentricFixed() {}
        GeocentricFixed(const GeocentricFixed& x) : 
            vec(x.vec){}
        GeocentricFixed(GeocentricFixed&& x) :
            vec(std::move(x.vec)) {}

        void operator=(const GeocentricFixed& x)
        {
            vec = x.vec;
        }
    };

    union GeocentricInertial
    {
        Eigen::Vector3d vec;
        struct
        {
            double x;
            double y;
            double z;
        };

        GeocentricInertial() {}
        GeocentricInertial(const GeocentricInertial& x) :
            vec(x.vec) {}
        GeocentricInertial(GeocentricInertial&& x) noexcept :
            vec(std::move(x.vec)) {}

        void operator=(const GeocentricInertial& x)
        {
            vec = x.vec;
        }
    };

    struct Geodetic
    {
        double latitude;
        double longitude;
        double altitude;
    };

    union ENU
    {
        Eigen::Vector3d vec;
        struct
        {
            double east;
            double north;
            double up;
        };

        ENU() {}
        ENU(const ENU& x) :
            vec(x.vec) {}
        ENU(ENU&& x) noexcept :
            vec(std::move(x.vec)) {}

        void operator=(const ENU& x)
        {
            vec = x.vec;
        }
    };

    struct Spherical
    {
        double radius;
        double polar;
        double azimuth;
    };

}

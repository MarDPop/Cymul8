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
    };

    class Geodetic
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
    };

    struct Spherical
    {
        double radius;
        double polar;
        double azimuth;
    };

}

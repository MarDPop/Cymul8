#pragma once

#include "../../lib/Eigen/Dense"

#include "../physics/Coordinates.h"

/**
* struct to contain an action to a body
* ie, force and moment.
* all units are in SI
*/
struct BodyAction
{
    /**
    * Force Vector in body frame
    */
    Eigen::Vector3d force;

    /**
    * Force Vector in body frame
    */
    Eigen::Vector3d moment;

    /**
    * Force Vector reference point in body
    */
    Eigen::Vector3d location;

    void setZero()
    {
        this->force.setZero();
        this->moment.setZero();
        this->location.setZero();
    }

    Eigen::Vector3d get_torque(const Eigen::Vector3d& location) const
    {
        Eigen::Vector3d arm = location - this->location;
        return this->moment + arm.cross(this->force);
    }

    BodyAction operator+(const BodyAction& otherAction) const;

    void operator+=(const BodyAction& otherAction)
    {
        this->force += otherAction.force;
        this->moment += otherAction.get_torque(this->location);
    }
};


/**
* struct to contain an action to a body
* ie, force and moment.
* all units are in SI
*/
struct GeneralAction
{
    /**
    * Force Vector in frame
    */
    Eigen::Vector3d force;

    /**
    * Moment Vector in frame
    */
    Eigen::Vector3d moment;

    /**
    * Location at which the force and moment are action from the origin of reference frame
    */
    Eigen::Vector3d arm;

    /**
    * if nullptr, inertial frame
    */
    Coordinate::Frame* reference_frame = nullptr;

    void setZero()
    {
        this->force.setZero();
        this->moment.setZero();
        this->arm.setZero();
    }

    Eigen::Vector3d get_torque() const
    {
        return this->moment + this->arm.cross(this->force);
    }

    GeneralAction operator+(const GeneralAction& otherAction) const;

    void operator+=(const GeneralAction& otherAction);
};

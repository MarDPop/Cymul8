#include "Action.h"

BodyAction BodyAction::operator+(const BodyAction& otherAction) const
{
    BodyAction output;
    output.location.setZero();

    Eigen::Vector3d thisTorque = this->location.cross(this->force);
    Eigen::Vector3d otherTorque = otherAction.location.cross(otherAction.force);

    for (int i = 0; i < 3; i++)
    {
        output.force[i] = this->force[i] + otherAction.force[i];
        output.moment[i] = this->moment[i] + otherAction.moment[i];
    }
    for (int i = 0; i < 3; i++)
    {
        output.moment[i] += thisTorque[i] + otherTorque[i];
    }
    return output;
}


GeneralAction GeneralAction::operator+(const GeneralAction& otherAction) const
{
    GeneralAction inertialAction;
    inertialAction.arm.setZero();
    if (this->reference_frame)
    {
        inertialAction.force = this->reference_frame->orientation.transpose() * this->force;
        inertialAction.moment = this->reference_frame->orientation.transpose() * this->get_torque();
    }
    else
    {
        inertialAction.force = this->force;
        inertialAction.moment = this->get_torque();
    }
    if (otherAction.reference_frame)
    {
        inertialAction.force += otherAction.reference_frame->orientation.transpose() * otherAction.force;
        inertialAction.moment += otherAction.reference_frame->orientation.transpose() * otherAction.get_torque();
    }
    else
    {
        inertialAction.force += otherAction.force;
        inertialAction.moment += otherAction.get_torque();
    }
    return inertialAction;
}

void GeneralAction::operator+=(const GeneralAction& otherAction)
{
    if (otherAction.reference_frame == this->reference_frame)
    {
        this->force += otherAction.force;
        this->moment += otherAction.get_torque();
    }
    else
    {
        if (!this->reference_frame)
        {
            this->force += otherAction.reference_frame->orientation.transpose() * otherAction.force;
            this->moment += otherAction.reference_frame->orientation.transpose() * otherAction.get_torque();
        }
        else
        {
            // Don't do this... lol
            // Also find a better way to do this
            Eigen::Vector3d forceInertial = otherAction.reference_frame->orientation.transpose() * otherAction.force;
            Eigen::Vector3d momentInertial = otherAction.reference_frame->orientation.transpose() * otherAction.get_torque();
            Eigen::Vector3d armInertial = otherAction.reference_frame->orientation.transpose() * otherAction.get_torque();
        }
    }
}
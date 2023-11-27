#include "Gravity.h"

Eigen::Vector3d Gravity::fictional_forces_Z_rotation(const Eigen::Vector3d& position,
    const Eigen::Vector3d& velocity,
    const double& rotation)
{
    double x = -rotation*(rotation*position[0] - 2.0*velocity[1]);
    double y = -rotation*(rotation*position[1] + 2.0*velocity[0]);

    return Eigen::Vector3d(x,y,0);
}

Eigen::Vector3d Gravity::fictional_forces_const_rotation(const Eigen::Vector3d& position,
    const Eigen::Vector3d& velocity,
    const Eigen::Vector3d& rotation)
{
    return (rotation.cross(position) + 2*velocity).cross(rotation); // negates by reverse cross product
}

Eigen::Vector3d Gravity::fictional_forces(const Eigen::Vector3d& position,
    const Eigen::Vector3d& velocity,
    const Eigen::Vector3d& rotation,
    const Eigen::Vector3d& rotation_rate)
{
    Eigen::Vector3d coriolis = -2.0 * rotation.cross(velocity);
    Eigen::Vector3d centripetal = rotation.cross(rotation.cross(position));
    Eigen::Vector3d euler = rotation_rate.cross(position);
    return coriolis - centripetal - euler;
}

void Gravity_Local::compute(const Coordinate::GeocentricFixed& position,
    const double R,
    Eigen::Vector3d& acceleration) const
{
    double ratio = _R0/(_R0 + position.z);
    acceleration[2] = _g0*ratio*ratio;
}

void Gravity_Kepler::compute(const Coordinate::GeocentricFixed& position,
    const double R,
    Eigen::Vector3d& acceleration) const
{
    acceleration = position.vec * (-_MU / (R*R*R));
}
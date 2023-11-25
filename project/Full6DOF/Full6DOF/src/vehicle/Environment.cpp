#include "Environment.h"

void AeroData::compute(const Air& air,
    const Eigen::Vector3d& air_ecef_velocity)
{
    airspeed = air_ecef_velocity.norm();
    air_velocity_ecef_unit = air_ecef_velocity*(1.0/airspeed);
    mach = airspeed*air.inv_sound_speed;
    dynamic_pressure = 0.5*air.density*airspeed*airspeed;
    double beta = 1.0 + 0.2*mach*mach;
    impact_pressure = air.pressure*(1.0 - beta*beta*beta*sqrt(beta));
}

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


void Environment::update(const Eigen::Vector3d& position,
    const Eigen::Vector3d& velocity,
    double time,
    bool near_body)
{
    _ref.TALO = time;

    _current_planet->set_ut1_jd2000(_ref.jd2000_launch + );

    if (near_body)
    {
        _near_body->update(*this);
        _current_planet->gravity().compute(_near_body->get_coordinates().PCF, time, _frame_acceleration);
    }
    else
    {

    }
}
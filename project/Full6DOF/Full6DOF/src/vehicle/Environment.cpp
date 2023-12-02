#include "Environment.h"

void AeroEnvironment::compute(const Air& air,
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

    double jd2000_utc = _ref.jd2000_utc_launch + time*Time::SEC_2_DAY_FRACTION;

    _planetary_system->set_utc_jd2000(jd2000_utc);

    _current_planet->gravity().compute(_near_body.get_coordinates().PCF, 
        _near_body.get_coordinates().RTP.radius,
        _frame_acceleration);

    
}
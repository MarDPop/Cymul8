#pragma once

#include "GNC.h"
#include "Environment.h"
#include "../physics/Body.h"
#include "../physics/Action.h"

#include "aerodynamics/Aerodynamics.h"
#include "propulsion/Thruster.h"
#include "propulsion/Propulsion.h"

template<class Body, class G>
class Vehicle : public virtual Body
{
    static_assert(std::is_base_of<Body>, Body_Base > ::value, ""G"Body not derived from Body_Base");
    static_assert(std::is_base_of<GNC<Body>>, G>::value, "G not derived from Guidance");

protected:

    friend class Vehicle_ODE;

    G _gnc;

    Environment _environment;

    virtual void update_accelerations(double t) = 0;

public:

    virtual unsigned get_num_states() const
    {
        _state_vector.size() + gnc.get_control_states(); // component states
    }

    const G& get_GNC() const
    {
        return _gnc;
    }

    const Environment& get_environment() const
    {
        return _environment;
    }

    void operator()(const Float* x, const Float t, Float* dx)
    {
        set_state(x);

        _environment.update(_state.position, _state.velocity, t);

        _gnc.update(x + _state_vector.size(), t, dx + _state_vector.size());

        _state.acceleration = _environment.get_frame_acceleration();

        this->update_accelerations(t);

        get_state_rate(dx);
    }

};


template<class T, class A, class GNC>
class Vehicle_3DOF_T : public virtual Vehicle<Body_Point_Mass<double>, GNC>
{
    static_assert(std::is_base_of<Thruster>, T > ::value, "P not derived from Propulsion");
    static_assert(std::is_base_of<Aerodynamics<Eigen::Vector3d>>, A > ::value, "A not derived from Aerodynamics");

protected:

    Eigen::Matrix3d _orientation;

    T _thruster;

    A _aero;

    void update_accelerations(double t)
    {
        Eigen::Vector3d force;
        force.setZero();

        if (_thruster.is_active())
        {
            _thruster.update(_environment.get_air(), _environment.get_aero_data(), t);
            force += _orientation.col(0)*_thruster.get_thrust();
        }

        if (_environment.in_air())
        {
            force += _aero.update(_environment.get_air(), _environment.get_aero_data(), t);
        }
          
        _acceleration += force*(1.0/_state.mass);
    }

public:

    void set_orientation(const Eigen::Matrix3d& __orientation)
    {
        _orientation = __orientation;
    }
};


template<MOMENT_CONSTANTS NDEG, class GNC>
class Vehicle_6DOF_Mass_Dependent : public virtual Vehicle< Body_Mass_Dependent_Inertia<NDEG>, GNC>
{
protected:

    Inertia<NDEG> _inertia_empty;

    std::unique_ptr<Aerodynamics<BodyAction>> _aerodynamics;

    std::unique_ptr<Propulsion> _propulsion;

    void update_accelerations(double t)
    {
        if (_propulsion->get_thruster().is_active())
        {
            _propulsion->update_thrust(_environment.get_air(), _environment.get_aero_data(), t);
        }

        if (_environment.in_air())
        {
            _aerodynamics->update(_environment.get_air(), _environment.get_aero_data(), t);
        }

        BodyAction totalAction = _propulsion->get_action() + _aerodynamics->get_action(); // at zero

        _acceleration += totalAction.force * (1.0 / _state.mass);

        _inertia.moment_of_inertia.get_angular_acceleration_body(_state.angular_velocity,
            totalAction.get_torque(_inertia.center_of_mass),
            _angular_acceleration);
        
        _mass_rate = -_propulsion->get_thruster().get_mass_rate();

        _propulsion->update_inertia(_state.mass, t, _state.acceleration);
    }


public:

    void set_aerodynamics(std::unique_ptr<Aerodynamics> __aero)
    {
        _aero = std::move(__aero);
    }

    void set_propulsion(std::unique_ptr<Aerodynamics> __propulsion)
    {
        _propulsion = std::move(__propulsion);
    }

};

template<class A, class P, MOMENT_CONSTANTS NDEG, class GNC>
class Vehicle_6DOF_Mass_Dependent_T : public virtual Vehicle<Body_Mass_Dependent_Inertia<NDEG>, GNC>
{
    static_assert(std::is_base_of<Aerodynamics<BodyAction>>, A > ::value, "A not derived from Aerodynamics");
    static_assert(std::is_base_of<Propulsion>, P > ::value, "P not derived from Propulsion");

protected:

    A _aerodynamics;

    P _propulsion;

    void update_accelerations(double t)
    {
        _acceleration.setZero();

        if (_thruster.is_active())
        {
            _thruster.update(_environment.get_air(), _environment.get_aero_data(), t);
            _acceleration += _orientation.col(0) * _thruster.get_thrust();
        }

        if (_environment.in_air())
        {
            _acceleration += _aero.update(_environment.get_air(), _environment.get_aero_data(), t);
        }

        _acceleration *= (1.0 / _state.mass);
    }

};

template<class A, class P, MOMENT_CONSTANTS NDEG, class GNC>
class Vehicle_6DOF_T : public virtual Vehicle< Body<NDEG>, GNC>
{
    static_assert(std::is_base_of<Aerodynamics<BodyAction>>, A > ::value, "A not derived from Aerodynamics");
    static_assert(std::is_base_of<Propulsion>, P > ::value, "P not derived from Propulsion");

protected:

    A _aerodynamics;

    P _propulsion;

    void update_accelerations(double t)
    {
        _acceleration.setZero();

        if (_thruster.is_active())
        {
            _thruster.update(_environment.get_air(), _environment.get_aero_data(), t);
            _acceleration += _orientation.col(0) * _thruster.get_thrust();
        }

        if (_environment.in_air())
        {
            _acceleration += _aero.update(_environment.get_air(), _environment.get_aero_data(), t);
        }

        _acceleration *= (1.0 / _state.mass);
    }

};
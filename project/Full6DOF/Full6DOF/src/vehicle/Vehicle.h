#pragma once

#include "GNC.h"
#include "Environment.h"
#include "../physics/Body.h"
#include "../physics/Action.h"

#include "aerodynamics/Aerodynamics.h"
#include "propulsion/Thruster.h"
#include "propulsion/Propulsion.h"
#include "component/Component.h"

// Intent is for this to be a vehicle factory... so you "build" vehicles in 
// static libraries and theses are just tools to help you build them

template<class Body, class G>
class Vehicle : public virtual Body
{
    static_assert(std::is_base_of<GNC<Body>, G>::value, "G not derived from Guidance");

protected:

    friend class Vehicle_ODE;

    G _gnc;

    Environment _environment;

    virtual void update_accelerations(double t) = 0;

public:

    virtual unsigned get_num_states() const
    {
        _state_vector.size() + gnc.get_control().N_CONTROL_STATES; // component states
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

    void update_accelerations(double t) override
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


template<class A, class P, MOMENT_CONSTANTS NDEG, class GNC>
class Vehicle_6DOF_T : public virtual Vehicle< Body_Mass_Dependent_Inertia<NDEG>, GNC>
{
    static_assert(std::is_base_of<Aerodynamics<BodyAction>>, A > ::value, "A not derived from Aerodynamics");
    static_assert(std::is_base_of<Propulsion>, P > ::value, "P not derived from Propulsion");

protected:

    A _aerodynamics;

    P _propulsion;

    Inertia<NDEG> _inertia_empty;

    /**
     * current center of mass location from nose (m) [should be negative]
     */
    Eigen::Vector3d _center_of_mass_delta;

    /**
    * current principal moment of inertia  (kg m2)
    */
    MomentOfInertia<NDEG> _moment_of_inertia_delta;


    void update_accelerations(double t) override
    {
        double dm = _state.mass - _inertia_empty.mass;
        _inertia.center_of_mass = _inertia_empty.center_of_mass + _center_of_mass_delta*dm;
        for (auto i = 0u; i < NDEG; i++)
        {
            _inertia.moment_of_inertia.I[i] = _inertia_empty.moment_of_inertia.I[i] + _moment_of_inertia_delta.I[i]*dm;
        }

        if (_propulsion.get_thruster().is_active())
        {
            _propulsion.update_thrust(_environment.get_air(), _environment.get_aero_data(), t);
        }

        if (_environment.in_air())
        {
            _aerodynamics.update(_environment.get_air(), _environment.get_aero_data(), t);
        }

        BodyAction totalAction = _propulsion.get_action() + _aerodynamics.get_action(); // at zero

        _acceleration += totalAction.force * (1.0 / _state.mass);

        _inertia.moment_of_inertia.get_angular_acceleration_body(_state.angular_velocity,
            totalAction.get_torque(_inertia.center_of_mass),
            _angular_acceleration);
        
        _mass_rate = -_propulsion.get_thruster().get_mass_rate();        
    }

public:

    void set_inertia(const Inertia<NDEG>& empty, const Inertia<NDEG>& full)
    {
        _inertia_empty = empty;
        double dm = full.mass - empty.mass;
        _center_of_mass_delta = (full.center_of_mass - empty.center_of_mass) / dm;
        for (auto i = 0u; i < NDEG; i++)
        {
            _moment_of_inertia_delta.I[i] = (full.moment_of_inertia.I[i] - empty.moment_of_inertia.I[i]) / dm;
        }
    }
};

template<class A, class P, class T, class GNC>
class Vehicle_6DOF_Full_T : public virtual Vehicle< Body<MOMENT_CONSTANTS::FULL>, GNC>
{
    static_assert(std::is_base_of<Aerodynamics<BodyAction>>, A > ::value, "A not derived from Aerodynamics");
    static_assert(std::is_base_of<Propulsion>, P > ::value, "P not derived from Propulsion");
    // could put static assert for NDEG == Tank::get_inertia()::size() instead of assuming full... but the templates are getting ridiculous

protected:

    A _aerodynamics;

    P _propulsion;

    T _tank;

    Inertia<MOMENT_CONSTANTS::FULL> _inertia_empty;

    Inertia<MOMENT_CONSTANTS::FULL> _inertia_vehicle;

    void update_accelerations(double t) override
    {
        if (_propulsion.get_thruster().is_active())
        {
            _propulsion.update_thrust(_environment.get_air(), _environment.get_aero_data(), t);
        }

        if (_environment.in_air())
        {
            _aerodynamics.update(_environment.get_air(), _environment.get_aero_data(), t);
        }

        BodyAction totalAction = _propulsion.get_action() + _aerodynamics.get_action(); // at zero

        _acceleration += totalAction.force * (1.0 / _state.mass);

        _state.inertia.moment_of_inertia.get_angular_acceleration_body(_state.angular_velocity,
            totalAction.get_torque(_state.inertia.center_of_mass),
            _angular_acceleration.data());

        _mass_rate = -_propulsion.get_thruster().get_mass_rate();

        _tank.update_inertia(_state.mass, t, _acceleration, _inertia_rate);
    }

};

class Vehicle_3DOF_Standard : public virtual Vehicle<Body_Point_Mass<double>, 
    GuidanceNavigationControl<Body_Point_Mass<double>>>
{
protected:

    Eigen::Matrix3d _orientation;

    std::unique_ptr<Thruster> _thruster = std::make_unique<Thruster>();

    std::unique_ptr<Aerodynamics<Eigen::Vector3d>> _aerodynamics = std::make_unique<Aerodynamics<Eigen::Vector3d>>();

    void update_accelerations(double t) override;

public:

    void set_orientation(const Eigen::Matrix3d& __orientation)
    {
        _orientation = __orientation;
    }

    void set_aerodynamics(std::unique_ptr<Aerodynamics<Eigen::Vector3d>> __aerodynamics)
    {
        _aerodynamics = std::move(__aerodynamics);
    }

    void set_thruster(std::unique_ptr<Thruster> __thruster)
    {
        _thruster = std::move(__thruster);
    }

    void set_guidance(std::unique_ptr<Guidance<Body_Point_Mass<double>>> __guidance)
    {
        _gnc.set_guidance(std::move(__guidance));
    }

    void set_navigation(std::unique_ptr<Navigation<Body_Point_Mass<double>>> __navigation)
    {
        _gnc.set_navigation(std::move(__navigation));
    }

    void set_control(std::unique_ptr<Control<Body_Point_Mass<double>>> __control)
    {
        _gnc.set_control(std::move(__control));
    }
};

class Vehicle_6DOF_Standard : 
    public virtual Vehicle<Body_Mass_Dependent_Inertia<MOMENT_CONSTANTS::FULL>, 
    GuidanceNavigationControl<Body_Mass_Dependent_Inertia<MOMENT_CONSTANTS::FULL>>>
{
    // LMAO at that template
protected:

    std::unique_ptr<Aerodynamics<BodyAction>> _aerodynamics;

    std::unique_ptr<Propulsion> _propulsion;

    Inertia<MOMENT_CONSTANTS::FULL> _inertia_empty;

    /**
     * current center of mass location from nose (m) [should be negative]
     */
    Eigen::Vector3d _center_of_mass_delta;

    /**
    * current principal moment of inertia  (kg m2)
    */
    MomentOfInertia<MOMENT_CONSTANTS::FULL> _moment_of_inertia_delta;


    void update_accelerations(double t) override;

public:

    void set_inertia(   const Inertia<MOMENT_CONSTANTS::FULL>& empty, 
                        const Inertia<MOMENT_CONSTANTS::FULL>& full);

    void set_aerodynamics(std::unique_ptr<Aerodynamics<BodyAction>> __aerodynamics)
    {
        _aerodynamics = std::move(__aerodynamics);
    }

    void set_propulsion(std::unique_ptr<Propulsion> __propulsion)
    {
        _propulsion = std::move(__propulsion);
    }

    void set_guidance(std::unique_ptr<Guidance<Body_Mass_Dependent_Inertia<MOMENT_CONSTANTS::FULL>>> __guidance)
    {
        _gnc.set_guidance(std::move(__guidance));
    }

    void set_navigation(std::unique_ptr<Navigation<Body_Mass_Dependent_Inertia<MOMENT_CONSTANTS::FULL>>> __navigation)
    {
        _gnc.set_navigation(std::move(__navigation));
    }

    void set_control(std::unique_ptr<Control<Body_Mass_Dependent_Inertia<MOMENT_CONSTANTS::FULL>>> __control)
    {
        _gnc.set_control(std::move(__control));
    }

};

class Vehicle_Components : public virtual Vehicle<Body<MOMENT_CONSTANTS::FULL>, GuidanceNavigationControl<Body<MOMENT_CONSTANTS::FULL>>>
{

protected:

    std::vector<Component> _components;

    void update_accelerations(double t) override;

public:

    void set_guidance(std::unique_ptr<Guidance<Body<MOMENT_CONSTANTS::FULL>>> __guidance)
    {
        _gnc.set_guidance(std::move(__guidance));
    }

    void set_navigation(std::unique_ptr<Navigation<Body<MOMENT_CONSTANTS::FULL>>> __navigation)
    {
        _gnc.set_navigation(std::move(__navigation));
    }

    void set_control(std::unique_ptr<Control<Body<MOMENT_CONSTANTS::FULL>>> __control)
    {
        _gnc.set_control(std::move(__control));
    }
};

class Staged_6DOF
{
    std::vector<Vehicle_6DOF_Standard> _stages;

public:


};
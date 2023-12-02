#pragma once

#include "GNC.h"
#include "Environment.h"
#include "../physics/Body.h"
#include "../physics/Action.h"
#include "../simulation/Simulation.h"

#include "aerodynamics/Aerodynamics.h"
#include "propulsion/Thruster.h"
#include "propulsion/Propulsion.h"
#include "component/Component.h"

// Intent is for this to be a vehicle factory... so you "build" vehicles in 
// static libraries and theses are just tools to help you build them

template<class B, class G>
class Vehicle : public virtual B
{
    static_assert(std::is_base_of<GNC, G>::value, "G not derived from Guidance");

protected:

    friend class Vehicle_ODE;

    G _gnc;

    Environment _environment;

    Simulation* _sim;

    virtual void update_accelerations(double t) = 0;

public:

    Vehicle() {}
    virtual ~Vehicle() {}

    void set_simulation(Simulation* sim)
    {
        _sim = sim;
    }

    virtual unsigned get_num_states() const
    {
        return static_cast<unsigned>(B::_state_vector.size())
            + _gnc.get_control().N_CONTROL_STATES; // component states
    }

    const G& get_GNC() const
    {
        return _gnc;
    }

    const Environment& get_environment() const
    {
        return _environment;
    }

    void operator()(const double* x, const double t, double* dx)
    {
        B::set_state(x);

        _environment.update(B::_state.position, B::_state.velocity, t);

        _gnc.update(x + B::_state_vector.size(), t, dx + B::_state_vector.size());

        B::_state.acceleration = _environment.get_frame_acceleration();

        this->update_accelerations(t);

        B::get_state_rate(dx);
    }

};

template<class T, class A, class GNC>
class Vehicle_3DOF_T final : public virtual Vehicle<Body_Point_Mass, GNC>
{
    static_assert(std::is_base_of<Thruster, T>::value, "P not derived from Propulsion");
    static_assert(std::is_base_of<Aerodynamics<Eigen::Vector3d>, A>::value, "A not derived from Aerodynamics");

protected:

    Eigen::Matrix3d _orientation;

    T _thruster;

    A _aero;

    void update_accelerations(double t) final override
    {
        Eigen::Vector3d force;
        force.setZero();

        if (_thruster.is_active())
        {
            _thruster.update(t);
            force += _orientation.col(0)*_thruster.get_thrust();
        }

        if (Vehicle<Body_Point_Mass, GNC>::_environment.in_air())
        {
            force += _aero.update(Vehicle<Body_Point_Mass, GNC>::_environment.get_aero_data(), t);
        }
          
        Body_Point_Mass::_acceleration += force*(1.0/ Body_Point_Mass::_state.mass);
    }

public:

    Vehicle_3DOF_T() {}

    void set_orientation(const Eigen::Matrix3d& __orientation)
    {
        _orientation = __orientation;
    }
};


template<class A, class P, MOMENT_CONSTANTS NDEG, class G>
class Vehicle_6DOF_T final : public virtual Vehicle< Body_Mass_Dependent_Inertia<NDEG>, G>
{
    static_assert(std::is_base_of<Aerodynamics<BodyAction>, A > ::value, "A not derived from Aerodynamics");
    static_assert(std::is_base_of<Propulsion, P >::value, "P not derived from Propulsion");

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

    void update_inertia() final override
    {
        double dm = Body_Mass_Dependent_Inertia<NDEG>::_state.mass - _inertia_empty.mass;
        Body_Mass_Dependent_Inertia<NDEG>::_inertia.center_of_mass = _inertia_empty.center_of_mass + _center_of_mass_delta * dm;
        for (auto i = 0u; i < NDEG; i++)
        {
            Body_Mass_Dependent_Inertia<NDEG>::_inertia.moment_of_inertia.I[i] =
                _inertia_empty.moment_of_inertia.I[i] + _moment_of_inertia_delta.I[i] * dm;
        }
    }

    void update_accelerations(double t) final override
    {
        if (_propulsion.get_thruster().is_active())
        {
            _propulsion.update_thrust(t);
        }

        if (Vehicle< Body_Mass_Dependent_Inertia<NDEG>, G>::_environment.in_air())
        {
            _aerodynamics.update(Vehicle< Body_Mass_Dependent_Inertia<NDEG>, G>::_environment.get_air(),
                Vehicle< Body_Mass_Dependent_Inertia<NDEG>, G>::_environment.get_aero_data(), t);
        }

        BodyAction totalAction = _propulsion.get_action() + _aerodynamics.get_action(); // at zero

        Body_Mass_Dependent_Inertia<NDEG>::_acceleration += totalAction.force * (1.0 / Body_Mass_Dependent_Inertia<NDEG>::_state.mass);

        Body_Mass_Dependent_Inertia<NDEG>::get_angualar_acceleration(
            totalAction.get_torque(Body_Mass_Dependent_Inertia<NDEG>::_inertia.center_of_mass));
        
        Body_Mass_Dependent_Inertia<NDEG>::_mass_rate = -_propulsion.get_thruster().get_mass_rate();
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

template<class A, class P, class T, class G>
class Vehicle_6DOF_Full_T final : public virtual Vehicle< Body<MOMENT_CONSTANTS::FULL>, G>
{
    static_assert(std::is_base_of<Aerodynamics<BodyAction>, A>::value, "A not derived from Aerodynamics");
    static_assert(std::is_base_of<Propulsion, P>::value, "P not derived from Propulsion");
    // could put static assert for NDEG == Tank::get_inertia()::size() instead of assuming full... but the templates are getting ridiculous

protected:

    A _aerodynamics;

    P _propulsion;

    T _tank;

    Inertia<MOMENT_CONSTANTS::FULL> _inertia_empty;

    Inertia<MOMENT_CONSTANTS::FULL> _inertia_vehicle;

    void update_accelerations(double t) final override
    {
        if (_propulsion.get_thruster().is_active())
        {
            _propulsion.update_thrust(Vehicle< Body<MOMENT_CONSTANTS::FULL>, G>::_environment.get_air(), 
                Vehicle< Body<MOMENT_CONSTANTS::FULL>, G>::_environment.get_aero_data(), t);
        }

        if (Vehicle< Body<MOMENT_CONSTANTS::FULL>, G>::_environment.in_air())
        {
            _aerodynamics.update(Vehicle< Body<MOMENT_CONSTANTS::FULL>, G>::_environment.get_air(), 
                Vehicle< Body<MOMENT_CONSTANTS::FULL>, G>::_environment.get_aero_data(), t);
        }

        BodyAction totalAction = _propulsion.get_action() + _aerodynamics.get_action(); // at zero

        Body<MOMENT_CONSTANTS::FULL>::_acceleration += totalAction.force * (1.0 / Body<MOMENT_CONSTANTS::FULL>::_state.inertia.mass);

        Body<MOMENT_CONSTANTS::FULL>::get_angualar_acceleration(
            totalAction.get_torque(Body<MOMENT_CONSTANTS::FULL>::_inertia.center_of_mass));

        Body<MOMENT_CONSTANTS::FULL>::_mass_rate = -_propulsion.get_thruster().get_mass_rate();

        _tank.update_inertia(Body<MOMENT_CONSTANTS::FULL>::_state.inertia.mass, t, 
            Body<MOMENT_CONSTANTS::FULL>::_acceleration, Body<MOMENT_CONSTANTS::FULL>::_inertia_rate);
    }

};

class Vehicle_3DOF_Standard final : public virtual Vehicle<Body_Point_Mass, 
    GuidanceNavigationControl<State_Point>>
{
protected:

    Eigen::Matrix3d _orientation;

    std::unique_ptr<Thruster> _thruster = std::make_unique<Thruster>();

    std::unique_ptr<Aerodynamics<Eigen::Vector3d>> _aerodynamics = std::make_unique<Aerodynamics<Eigen::Vector3d>>();

    void update_accelerations(double t) override;

public:

    Vehicle_3DOF_Standard() {}

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

    void set_guidance(std::unique_ptr<Guidance<State_Point>> __guidance)
    {
        _gnc.set_guidance(std::move(__guidance));
    }

    void set_navigation(std::unique_ptr<Navigation<State_Point>> __navigation)
    {
        _gnc.set_navigation(std::move(__navigation));
    }

    void set_control(std::unique_ptr<Control<State_Point>> __control)
    {
        _gnc.set_control(std::move(__control));
    }
};

class Vehicle_6DOF_Standard final : 
    public virtual Vehicle<Body_Mass_Dependent_Inertia<MOMENT_CONSTANTS::FULL>, 
    GuidanceNavigationControl<State_Rigid_Body_<MOMENT_CONSTANTS::FULL>>>
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


    void update_accelerations(double t) final override;

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

    void set_guidance(std::unique_ptr<Guidance<State_Rigid_Body_<MOMENT_CONSTANTS::FULL>>> __guidance)
    {
        _gnc.set_guidance(std::move(__guidance));
    }

    void set_navigation(std::unique_ptr<Navigation<State_Rigid_Body_<MOMENT_CONSTANTS::FULL>>> __navigation)
    {
        _gnc.set_navigation(std::move(__navigation));
    }

    void set_control(std::unique_ptr<Control<State_Rigid_Body_<MOMENT_CONSTANTS::FULL>>> __control)
    {
        _gnc.set_control(std::move(__control));
    }

};

class Vehicle_Components : public virtual Vehicle<Body<MOMENT_CONSTANTS::FULL>, GuidanceNavigationControl<State_Rigid_Body_<MOMENT_CONSTANTS::FULL>>>
{

protected:

    std::vector<Component> _components;

    void update_accelerations(double t) final override;

public:

    void set_guidance(std::unique_ptr<Guidance<State_Rigid_Body_<MOMENT_CONSTANTS::FULL>>> __guidance)
    {
        _gnc.set_guidance(std::move(__guidance));
    }

    void set_navigation(std::unique_ptr<Navigation<State_Rigid_Body_<MOMENT_CONSTANTS::FULL>>> __navigation)
    {
        _gnc.set_navigation(std::move(__navigation));
    }

    void set_control(std::unique_ptr<Control<State_Rigid_Body_<MOMENT_CONSTANTS::FULL>>> __control)
    {
        _gnc.set_control(std::move(__control));
    }
};

class Staged_6DOF
{
    std::vector<Vehicle_6DOF_Standard> _stages;

public:


};
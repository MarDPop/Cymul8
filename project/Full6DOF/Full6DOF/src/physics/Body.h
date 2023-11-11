#pragma once

#include "../ode/Dynamics.h"
#include "../../lib/Eigen/Dense"
#include <cstring>
#include <array>

/**
 * @brief
 *
 */
 /*
 Note on constants:
     - for Plane Symmetry the plane is assumed to be in the XZ axis
     - for Axisymmetric the axis is assumed to be on the X axis
 */
enum MOMENT_CONSTANTS : unsigned
{
    EQUAL = 1u,
    AXISYMMETRIC = 2u,
    PRINCIPAL_AXIS = 3u,
    PLANE_SYMMETRY = 4u,
    FULL = 6u
};


template<MOMENT_CONSTANTS NDEG>
struct MomentOfInertia
{
    std::array<double, NDEG> I; // Order is Ixx, Iyy, Izz, Ixy, Ixz, Iyz
    // EXCEPT FOR PLANE SYMMETRY: Order is Ixx, Iyy, Izz, Ixz,

    MomentOfInertia operator+(const MomentOfInertia& moi) const
    {
        MomentOfInertia output;
        for (unsigned idx = 0; idx < NDEG; idx++)
        {
            output.I[idx] = I[idx] + moi.I[idx];
        }
        return output;
    }

    void operator=(const MomentOfInertia& moi)
    {
        memset(I.data(), moi.I.data(), NDEG * sizeof(double));
    }

    Eigen::Matrix3d get_inertia_matrix() const;

    void get_angular_acceleration_inertial(const Eigen::Quaterniond& orientation,
        const Eigen::Vector3d& angular_velocity_inertial,
        const Eigen::Vector3d& torque,
        double* angular_acceleration_inertial) const;

    void get_angular_acceleration_body(const Eigen::Vector3d& angular_velocity_body,
        const Eigen::Vector3d& torque_body,
        double* angular_acceleration_body) const;
};

template<MOMENT_CONSTANTS NDEG>
struct Inertia
{
    /**
    * current mass  (kg)
    */
    double mass;

    /**
    * current center of mass location from nose (m) [should be negative]
    */
    Eigen::Vector3d center_of_mass;

    /**
    * current principal moment of inertia  (kg m2)
    */
    MomentOfInertia<NDEG> moment_of_inertia;

    void operator=(const Inertia<NDEG>& inertia)
    {
        this->mass = inertia.mass;
        this->center_of_mass = inertia.center_of_mass;
        memcpy(this->moment_of_inertia.I.data(), inertia.moment_of_inertia.I.data(), NDEG * sizeof(double));
    }

    Inertia<MOMENT_CONSTANTS::FULL> operator+(const Inertia<NDEG>& inertia) const;

};

struct State_Point
{
    static constexpr unsigned N_STATES = 7;

    Eigen::Vector3d position;

    Eigen::Vector3d velocity;

    double mass;
};

struct State_Rigid_Body
{
    static constexpr unsigned N_STATES = 14;

    Eigen::Vector3d position;

    Eigen::Vector3d velocity;

    Eigen::Quaterniond orientation;

    Eigen::Vector3d angular_velocity;

    double mass;
};

template<MOMENT_CONSTANTS NDEG>
struct State_Rigid_Body_
{
    static constexpr unsigned N_STATES = 17 + NDEG;

    Eigen::Vector3d position;

    Eigen::Vector3d velocity;

    Eigen::Quaterniond orientation;

    Eigen::Vector3d angular_velocity;

    Inertia<NDEG> inertia;
};

template<typename State>
class Body_Base
{
protected:

    union
    {
        alignas(32) std::array<double, State::N_STATES> _state_vector;
        alignas(32) State _state;
    };

public:

    virtual void get_state_rate(double* dx) = 0;

    virtual void set_state(const double* x)
    {
        memcpy(_state_vector.data(), x, sizeof(_state_vector));
    }

    const double* get_state_vector() const
    {
        return _state_vector.data();
    }

    const State& get_state() const
    {
        return _state;
    }
};

namespace body
{
    void get_orientation_rate(const Eigen::Vector3d& angular_velocity,
        const Eigen::Quaterniond& orientation,
        double* q_dot);

    void get_angular_acceleration(const Eigen::Vector3d& angular_velocity,
        const Eigen::Matrix3d& I,
        const Eigen::Vector3d& torque,
        double* angular_acceleration);
}

/**
 * @brief 
 * 
 */
class Body_Point_Mass : public virtual Body_Base<State_Point>
{
protected:
    
    Eigen::Vector3d _acceleration;

    double _mass_rate;

public:

    void get_state_rate(double* dx) override
    {
        memcpy(dx, _state.velocity.data(), 3*sizeof(double));
        memcpy(dx + 3, _acceleration.data(), 3*sizeof(double));
        dx[6] = _mass_rate;
    }
};

template<MOMENT_CONSTANTS NDEG>
class Body_Mass_Dependent_Inertia : public virtual Body_Base<State_Rigid_Body>
{
protected:

    /**
     * current center of mass location from nose (m) [should be negative]
     */
    Eigen::Vector3d _center_of_mass;

    /**
    * current principal moment of inertia  (kg m2)
    */
    MomentOfInertia<NDEG> _moment_of_inertia;

    Eigen::Vector3d _acceleration;

    Eigen::Vector3d _angular_acceleration;

    double _mass_rate;

    virtual void update_inertia() = 0;

    void get_angualar_acceleration(const Eigen::Vector3d& torque)
    {
        _moment_of_inertia.get_angular_acceleration_body(
            _state.angular_velocity,
            torque,
            _angular_acceleration);
    }

public:

    void get_state_rate(double* dx) override
    {
        memcpy(dx, _state.velocity.data(), 3*sizeof(double));
        memcpy(dx + 3, _acceleration.data(), 3* sizeof(double));
        body::get_orientation_rate(_state.angular_velocity, _state.orientation, dx + 6);
        memcpy(dx + 10, _angular_acceleration.data(), 3*sizeof(double));
        dx[13] = _mass_rate;
    }

    void set_state(const double* x) override
    {
        memcpy(_state_vector.data(), x, sizeof(_state_vector));
        this->update_inertia();
    }
};

/**
 * @brief 
 * 
 */
template<MOMENT_CONSTANTS NDEG>
class Body : public virtual Body_Base<State_Rigid_Body_<NDEG>>
{
protected:

    Eigen::Vector3d _acceleration;

    Eigen::Vector3d _angular_acceleration;

    Inertia<NDEG> _inertia_rate;

    void get_angualar_acceleration(const Eigen::Vector3d& torque)
    {
        Body_Base<State_Rigid_Body_<NDEG>>::_inertia.moment_of_inertia.get_angular_acceleration_body(
            Body_Base<State_Rigid_Body_<NDEG>>::_state.angular_velocity,
            torque,
            _angular_acceleration);
    }

public:

    void get_state_rate(double* dx) override
    {
        memcpy(dx, Body_Base<State_Rigid_Body_<NDEG>>::_state.velocity.data(),3*sizeof(double));
        memcpy(dx + 3,_acceleration.data(),3*sizeof(double));
        body::get_orientation_rate(Body_Base<State_Rigid_Body_<NDEG>>::_state.angular_velocity,
            Body_Base<State_Rigid_Body_<NDEG>>::_state.orientation, dx + 6);
        memcpy(dx + 10,_angular_acceleration.data(),3*sizeof(double));
        dx[13] = _inertia_rate.mass;
        memcpy(dx + 14,_inertia_rate.center_of_mass.data(),3*sizeof(double));
        memcpy(dx + 17,_inertia_rate.moment_of_inertia.I.data(), NDEG*sizeof(double));
    }

};
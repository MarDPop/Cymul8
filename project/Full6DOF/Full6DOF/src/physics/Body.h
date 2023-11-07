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

    void operator=(const Inertia<NDEG>& inertia) const
    {
        this->mass = inertia.mass;
        this->center_of_mass = inertia.center_of_mass;
        memcpy(this->moment_of_inertia.I.data(), inertia.moment_of_inertia.I.data(), NDEG * sizeof(double));
    }

    Inertia<MOMENT_CONSTANTS::FULL> operator+(const Inertia<NDEG>& inertia) const;

};

template<typename Float, class S>
class Body_Base
{
protected:

    union
    {
        alignas(32) std::array<Float, S::N_STATES> _state_vector;
        alignas(32) S _state;
    };

public:

    virtual void get_state_rate(Float* dx) = 0;

    virtual void set_state(const Float* x)
    {
        memcpy(_state_vector.data(), x, sizeof(_state_vector));
    }

    const Float* get_state_vector() const
    {
        return _state_vector.data();
    }

    const S& get_state() const
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

template<typename Float>
struct State_Point
{
    static_assert(std::is_arithmetic<Float>::value, "Must use a arithmetic type");

    static constexpr unsigned N_STATES = 7;

    Eigen::Matrix<Float, 3, 1> position;

    Eigen::Matrix<Float, 3, 1> velocity;

    Float mass;
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

/**
 * @brief 
 * 
 */
template<typename Float>
class Body_Point_Mass : public virtual Body_Base<Float, State_Point<Float>>
{
protected:
    
    Eigen::Matrix<Float, 3, 1> _acceleration;

    Float _mass_rate;

public:

    void get_state_rate(Float* dx) override
    {
        memcpy(dx, velocity.data(), 3*sizeof(Float));
        memcpy(dx + 3, acceleration.data(), 3*sizeof(Float));
        dx[6] = mass_rate;
    }
};

template<MOMENT_CONSTANTS NDEG>
class Body_Mass_Dependent_Inertia : public virtual Body_Base<double, State_Rigid_Body>
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

public:

    void get_state_rate(double* dx) override
    {
        memcpy(dx, this->velocity.data(), 3 * sizeof(double));
        memcpy(dx + 3, this->acceleration.data(), 3 * sizeof(double));
        body::get_orientation_rate(angular_velocity, orientation, dx + 6);
        memcpy(dx + 10, this->angular_acceleration.data(), 3 * sizeof(double));
        dx[13] = this->_mass_rate;
    }

    void set_state(const Float* x) override
    {
        memcpy(state.data(), x, sizeof(state));
        this->update_inertia();
    }
};

/**
 * @brief 
 * 
 */
template<MOMENT_CONSTANTS NDEG>
class Body : public virtual Body_Base<double, State_Rigid_Body_<NDEG>>
{
protected:

    Eigen::Vector3d _acceleration;

    Eigen::Vector3d _angular_acceleration;

    Inertia<NDEG> _inertia_rate;

public:

    void get_state_rate(double* dx) override
    {
        memcpy(dx,this->_state.velocity.data(),3*sizeof(double));
        memcpy(dx + 3,this->_acceleration.data(),3*sizeof(double));
        body::get_orientation_rate(angular_velocity, orientation, dx + 6);
        memcpy(dx + 10,this->angular_acceleration.data(),3*sizeof(double));
        dx[13] = this->inertia_rate.mass;
        memcpy(dx + 14,this->inertia_rate.center_of_mass.data(),3*sizeof(double));
        memcpy(dx + 17,this->inertia_rate.moment_of_inertia.I.data(), NDEG*sizeof(double));
    }

};
#pragma once

#include "../ode/Dynamics.h"
#include "../../lib/Eigen/Dense"
#include <cstring>
#include <array>

template<typename Float>
struct Body_Base
{
    virtual void get_state_rate(Float* dx) = 0;
};

/**
 * @brief 
 * 
 */
template<typename Float>
struct Body_Point_Mass : Body_Base<Float>
{
    static_assert(std::is_arithmetic<Float>::value, "Must use a arithmetic type");

    static constexpr unsigned N_STATES = 7;

    union
    {
        alignas(32) std::array<Float, N_STATES> state;
        struct alignas(32)
        {
            Eigen::Matrix<Float, 3, 1> position;

            Eigen::Matrix<Float, 3, 1> velocity;

            Float mass;
        };
    };
    
    Eigen::Matrix<Float, 3, 1> acceleration;

    Float mass_rate;

    void get_state_rate(Float* dx) override
    {
        memcpy(dx, velocity.data(), 3*sizeof(Float));
        memcpy(dx + 3, acceleration.data(), 3*sizeof(Float));
        dx[6] = mass_rate;
    }
};

/**
 * @brief
 *
 */
class Body_Point_Mass_Dynamics : virtual protected Body_Point_Mass<double>,
                                    virtual Fixed_Size_Dynamics<7>
{
protected:

    double time;

    virtual void compute_acceleration() {}

    virtual void compute_mass_rate() {}

    virtual bool stop_conditions()
    {
        return true;
    }

public:

    bool set_state(const std::array<double, 7>& x, const double& time, std::array<double, 7>& dx);

};

/**
 * @brief 
 * 
 */
/*
Note on constants:
    - for Plane Symmetry the plane is assumed to be in the XZ axis
    - for Axisymmetric the axis is assumed to be on the Z axis
*/
enum class MOMENT_CONSTANTS : unsigned 
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

    MomentOfInertia operator+(const MomentOfInertia& moi) const
    {
        MomentOfInertia output;
        for(unsigned idx = 0; idx < NDEG; idx++)
        {
            output.I[idx] = I[idx] + moi.I[idx];
        }
        return output;
    }

    void operator=(const MomentOfInertia& moi)
    {
        memset(I.data(), moi.I.data(), NDEG*sizeof(double));
    }

    Eigen::Matrix3d get_inertia_matrix() const;
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
        memcpy(this->moment_of_inertia.I.data(),inertia.moment_of_inertia.I.data(),NDEG*sizeof(double));
    }

    Inertia<MOMENT_CONSTANTS::FULL> operator+(const Inertia<NDEG>& inertia) const;

};


namespace body
{
    void get_orientation_rate(const Eigen::Vector3d& angular_velocity,
        const Eigen::Quaterniond& orientation,
        double* q_dot);


}

/**
 * @brief 
 * 
 */
template<MOMENT_CONSTANTS NDEG>
struct Body : virtual Body_Base<double>
{
    static constexpr unsigned N_STATES = 17 + NDEG;

    union 
    {
        alignas(32) std::array<double, 17 + NDEG> state;
        struct alignas(32)
        {
            Eigen::Vector3d position;

            Eigen::Vector3d velocity;

            Eigen::Quaterniond orientation;

            Eigen::Vector3d angular_velocity;

            Inertia<NDEG> inertia;
        };
    };

    Eigen::Vector3d acceleration;

    Eigen::Vector3d angular_acceleration;

    Inertia<NDEG> inertia_rate;

    void get_state_rate(double* dx) override
    {
        memcpy(dx,this->velocity.data(),3*sizeof(double));
        memcpy(dx + 3,this->acceleration.data(),3*sizeof(double));
        body::get_orientation_rate(angular_velocity, orientation, dx + 6);
        memcpy(dx + 10,this->angular_acceleration.data(),3*sizeof(double));
        dx[13] = this->inertia_rate.mass;
        memcpy(dx + 14,this->inertia_rate.center_of_mass.data(),3*sizeof(double));
        memcpy(dx + 17,this->inertia_rate.moment_of_inertia.I.data(), NDEG*sizeof(double));
    }
};


/**
 * @brief
 *
 */
template<MOMENT_CONSTANTS NDEG>
class Body_Dynamics : protected virtual Body<NDEG>,  Fixed_Size_Dynamics<17 + NDEG>
{
protected:

    double _time;

    virtual void compute_state_rate() = 0;

    virtual bool stop_conditions()
    {
        return true;
    }

public:

    bool set_state(const std::array<double, 17 + NDEG>& x, const double& time, std::array<double, 17 + NDEG>& dx) override
    {
        this->_state = x;
        this->_time = time;
        this->compute_state_rate();
        this->get_state_rate(dx.data());
        return this->stop_conditions();
    }

};
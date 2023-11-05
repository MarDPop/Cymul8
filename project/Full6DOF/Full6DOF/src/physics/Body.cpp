#include "Body.h"

bool Body_Point_Mass_Dynamics::set_state(const std::array<double, 7>& x, const double& time, std::array<double, 7>& dx)
{
    memcpy(this->position.data(), &x[0], 3 * sizeof(double));
    memcpy(this->velocity.data(), &x[3], 3 * sizeof(double));
    this->mass = x[6];
    this->time = time;

    this->compute_acceleration();
    this->compute_mass_rate();

    dx[0] = x[3];
    dx[1] = x[4];
    dx[2] = x[5];
    dx[3] = _acceleration[0];
    dx[4] = _acceleration[1];
    dx[5] = _acceleration[2];
    dx[6] = _mass_rate;

    return this->stop_conditions();
}


template<>
inline Eigen::Matrix3d MomentOfInertia<MOMENT_CONSTANTS::FULL>::get_inertia_matrix() const
{
    Eigen::Matrix3d inertia;
    double* data = inertia.data();
    data[0] = this->I[0];
    data[4] = this->I[1];
    data[8] = this->I[2];
    data[1] = data[3] = -this->I[3];
    data[2] = data[6] = -this->I[4];
    data[5] = data[7] = -this->I[5];
    return inertia;
}

template<>
inline Eigen::Matrix3d MomentOfInertia<MOMENT_CONSTANTS::PLANE_SYMMETRY>::get_inertia_matrix() const
{
    Eigen::Matrix3d inertia;
    double* data = inertia.data();
    data[0] = this->I[0];
    data[4] = this->I[1];
    data[8] = this->I[2];
    data[1] = data[3] = 0.0;
    data[2] = data[6] = -this->I[3];
    data[5] = data[7] = 0.0;
    return inertia;
}

template<>
inline Eigen::Matrix3d MomentOfInertia<MOMENT_CONSTANTS::PRINCIPAL_AXIS>::get_inertia_matrix() const
{
    Eigen::Matrix3d inertia;
    double* data = inertia.data();
    data[0] = this->I[0];
    data[4] = this->I[1];
    data[8] = this->I[2];
    data[1] = data[2] = data[3] = 0.0;
    data[5] = data[6] = data[7] = 0.0;
    return inertia;
}

template<>
inline Eigen::Matrix3d MomentOfInertia<MOMENT_CONSTANTS::AXISYMMETRIC>::get_inertia_matrix() const
{
    Eigen::Matrix3d inertia;
    double* data = inertia.data();
    data[0] = data[4] = this->I[0];
    data[8] = this->I[1];
    data[1] = data[2] = data[3] = 0.0;
    data[5] = data[6] = data[7] = 0.0;
    return inertia;
}

template<>
inline Eigen::Matrix3d MomentOfInertia<MOMENT_CONSTANTS::EQUAL>::get_inertia_matrix() const
{
    Eigen::Matrix3d inertia;
    double* data = inertia.data();
    data[0] = data[4] = data[8] = this->I[0];
    data[1] = data[2] = data[3] = 0.0;
    data[5] = data[6] = data[7] = 0.0;
    return inertia;
}

template<MOMENT_CONSTANTS NDEG>
Inertia<MOMENT_CONSTANTS::FULL> Inertia<NDEG>::operator+(const Inertia<NDEG>& inertia) const
{
    Inertia<MOMENT_CONSTANTS::FULL> output;

    // Get mass first
    output.mass = this->mass + inertia.mass;

    // Compute center of mass
    output.center_of_mass = this->center_of_mass * this->mass + inertia.center_of_mass * inertia.mass;
    output.center_of_mass *= (1.0 / output.mass);

    output.moment_of_inertia = this->moment_of_inertia + inertia.moment_of_inertia;

    Eigen::Vector3d r = this->center_of_mass - output.center_of_mass;
    Eigen::Vector3d mr = r * this->mass;
    double mr2 = mr.dot(r);
    output.moment_of_inertia.I[0] += mr2;
    output.moment_of_inertia.I[1] += mr2;
    output.moment_of_inertia.I[2] += mr2;

    output.moment_of_inertia.I[0] -= mr.x() * r.x();
    output.moment_of_inertia.I[1] -= mr.y() * r.y();
    output.moment_of_inertia.I[2] -= mr.z() * r.z();
    output.moment_of_inertia.I[3] -= mr.x() * r.y();
    output.moment_of_inertia.I[4] -= mr.x() * r.z();
    output.moment_of_inertia.I[5] -= mr.y() * r.z();

    r = inertia.center_of_mass - output.center_of_mass;
    mr = r * inertia.mass;
    mr2 = mr.dot(r);
    output.moment_of_inertia.I[0] += mr2;
    output.moment_of_inertia.I[1] += mr2;
    output.moment_of_inertia.I[2] += mr2;

    output.moment_of_inertia.I[0] -= mr.x() * r.x();
    output.moment_of_inertia.I[1] -= mr.y() * r.y();
    output.moment_of_inertia.I[2] -= mr.z() * r.z();
    output.moment_of_inertia.I[3] -= mr.x() * r.y();
    output.moment_of_inertia.I[4] -= mr.x() * r.z();
    output.moment_of_inertia.I[5] -= mr.y() * r.z();

    return output;
}
template<>
Inertia<MOMENT_CONSTANTS::FULL> Inertia<MOMENT_CONSTANTS::FULL>::operator+(const Inertia<MOMENT_CONSTANTS::FULL>& inertia) const;

template<>
Inertia<MOMENT_CONSTANTS::FULL> Inertia<MOMENT_CONSTANTS::PLANE_SYMMETRY>::operator+(const Inertia<MOMENT_CONSTANTS::PLANE_SYMMETRY>& inertia) const;

template<>
Inertia<MOMENT_CONSTANTS::FULL> Inertia<MOMENT_CONSTANTS::PRINCIPAL_AXIS>::operator+(const Inertia<MOMENT_CONSTANTS::PRINCIPAL_AXIS>& inertia) const;

template<>
Inertia<MOMENT_CONSTANTS::FULL> Inertia<MOMENT_CONSTANTS::AXISYMMETRIC>::operator+(const Inertia<MOMENT_CONSTANTS::AXISYMMETRIC>& inertia) const;

template<>
Inertia<MOMENT_CONSTANTS::FULL> Inertia<MOMENT_CONSTANTS::EQUAL>::operator+(const Inertia<MOMENT_CONSTANTS::EQUAL>& inertia) const;


#define METHOD 1
#if METHOD == 0
void body::get_orientation_rate(const Eigen::Vector3d& angular_velocity,
    const Eigen::Quaterniond& orientation,
    double* q_dot)
{

    Eigen::Matrix4d angular_matrix{
        {-angular_velocity.x(), -angular_velocity.y(), -angular_velocity.z(), 0.0},
        {0.0, angular_velocity.z(), -angular_velocity.y(), angular_velocity.x()},
        {-angular_velocity.z(), 0.0, angular_velocity.x(), angular_velocity.y()},
        {angular_velocity.y(), -angular_velocity.x(), 0.0, angular_velocity.z()}
    }; // outer product w x q

    Eigen::Map<Eigen::Vector4d> quat(orientation.coeffs().data());

    Eigen::Vector4d mult = angular_matrix * quat;

    auto q = mult.data();
    q_dot[0] = q[0] * 0.5;
    q_dot[1] = q[1] * 0.5;
    q_dot[2] = q[2] * 0.5;
    q_dot[3] = q[3] * 0.5;
}
#elif METHOD == 1 && defined(__AVX__)

#include <immintrin.h>

void body::get_orientation_rate(const Eigen::Vector3d& angular_velocity,
    const Eigen::Quaterniond& orientation,
    double* q_dot)
{
    static const __m256d half = _mm256_set1_pd(0.5);

    __m256d col1 = _mm256_set_pd(-angular_velocity.x(), 0.0, -angular_velocity.z(), angular_velocity.y());
    __m256d col2 = _mm256_set_pd(-angular_velocity.y(), angular_velocity.z(), 0.0, -angular_velocity.x());
    __m256d col3 = _mm256_set_pd(-angular_velocity.z(), -angular_velocity.y(), angular_velocity.x(), 0.0);
    __m256d col4 = _mm256_set_pd(0.0, angular_velocity.x(), angular_velocity.y(), angular_velocity.z());

    col1 = _mm256_mul_pd(col1, _mm256_set1_pd(orientation.x()));
    col2 = _mm256_mul_pd(col2, _mm256_set1_pd(orientation.y()));
    col3 = _mm256_mul_pd(col3, _mm256_set1_pd(orientation.z()));
    col4 = _mm256_mul_pd(col4, _mm256_set1_pd(orientation.w()));

    __m256d xy = _mm256_add_pd(col1, col2);
    __m256d zw = _mm256_add_pd(col3, col4);

    _mm256_storeu_pd(q_dot, _mm256_mul_pd(_mm256_add_pd(xy, zw), half)); // hopefull uses multadd
}

#else 
void body::get_orientation_rate(const Eigen::Vector3d& angular_velocity,
    const Eigen::Quaterniond& orientation,
    double* q_dot)
{
    const double* q = orientation.coeffs().data();

    q_dot[0] = (-angular_velocity.x() * q[0], -angular_velocity.y() * q[1], -angular_velocity.z() * q[2]) * 0.5;
    q_dot[1] = (angular_velocity.z() * q[1], -angular_velocity.y() * q[2], angular_velocity.x() * q[3]) * 0.5;
    q_dot[2] = (-angular_velocity.z() * q[0], angular_velocity.x() * q[2], angular_velocity.y() * q[3]) * 0.5;
    q_dot[3] = (angular_velocity.y() * q[0], -angular_velocity.x() * q[1], angular_velocity.z() * q[3]) * 0.5;
}

#endif
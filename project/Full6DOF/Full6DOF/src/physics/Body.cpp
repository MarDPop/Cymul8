#include "Body.h"

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

template<>
 Eigen::Matrix3d MomentOfInertia<MOMENT_CONSTANTS::FULL>::get_inertia_matrix() const
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
 Eigen::Matrix3d MomentOfInertia<MOMENT_CONSTANTS::PLANE_SYMMETRY>::get_inertia_matrix() const
{
     // symmetry about 
    Eigen::Matrix3d inertia;
    double* data = inertia.data();
    data[0] = this->I[0];
    data[4] = this->I[1];
    data[8] = this->I[2];
    data[1] = data[3] = 0.0;
    data[2] = data[6] = -this->I[3]; // 
    data[5] = data[7] = 0.0;
    return inertia;
}

template<>
 Eigen::Matrix3d MomentOfInertia<MOMENT_CONSTANTS::PRINCIPAL_AXIS>::get_inertia_matrix() const
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
 Eigen::Matrix3d MomentOfInertia<MOMENT_CONSTANTS::AXISYMMETRIC>::get_inertia_matrix() const
{
    Eigen::Matrix3d inertia;
    double* data = inertia.data();
    data[0] = this->I[0];
    data[4] = data[8] = this->I[1];
    data[1] = data[2] = data[3] = 0.0;
    data[5] = data[6] = data[7] = 0.0;
    return inertia;
}

template<>
 Eigen::Matrix3d MomentOfInertia<MOMENT_CONSTANTS::EQUAL>::get_inertia_matrix() const
{
    Eigen::Matrix3d inertia;
    double* data = inertia.data();
    data[0] = data[4] = data[8] = this->I[0];
    data[1] = data[2] = data[3] = 0.0;
    data[5] = data[6] = data[7] = 0.0;
    return inertia;
}

 template<>
 void MomentOfInertia<MOMENT_CONSTANTS::FULL>::get_angular_acceleration_inertial(const Eigen::Quaterniond& orientation,
     const Eigen::Vector3d& angular_velocity_inertial,
     const Eigen::Vector3d& torque_inertial,
     double* angular_acceleration_inertial) const
 {

 }

 template<>
 void MomentOfInertia<MOMENT_CONSTANTS::PLANE_SYMMETRY>::get_angular_acceleration_inertial(const Eigen::Quaterniond& orientation,
     const Eigen::Vector3d& angular_velocity_inertial,
     const Eigen::Vector3d& torque_inertial,
     double* angular_acceleration_inertial) const
 {

 }

 template<>
 void MomentOfInertia<MOMENT_CONSTANTS::PRINCIPAL_AXIS>::get_angular_acceleration_inertial(const Eigen::Quaterniond& orientation,
     const Eigen::Vector3d& angular_velocity_inertial,
     const Eigen::Vector3d& torque_inertial,
     double* angular_acceleration_inertial) const
 {

 }

 template<>
 void MomentOfInertia<MOMENT_CONSTANTS::AXISYMMETRIC>::get_angular_acceleration_inertial(const Eigen::Quaterniond& orientation,
     const Eigen::Vector3d& angular_velocity_inertial,
     const Eigen::Vector3d& torque_inertial,
     double* angular_acceleration_inertial) const
 {

 }

 template<>
 void MomentOfInertia<MOMENT_CONSTANTS::EQUAL>::get_angular_acceleration_inertial(const Eigen::Quaterniond& orientation,
     const Eigen::Vector3d& angular_velocity_inertial,
     const Eigen::Vector3d& torque_inertial,
     double* angular_acceleration_inertial) const
 {
     
 }

template<>
void MomentOfInertia<MOMENT_CONSTANTS::FULL>::get_angular_acceleration_body(const Eigen::Vector3d& angular_velocity,
     const Eigen::Vector3d& torque_body,
     double* angular_acceleration) const
{
    Eigen::Map< Eigen::Vector3d> omega_dot(angular_acceleration);
    Eigen::Matrix3d I = this->get_inertia_matrix();
    Eigen::Vector3d wxIw = angular_velocity.cross(I*angular_velocity);
    omega_dot = I.lu().solve(torque_body - wxIw); // see if LDLT could be used instead
}

template<>
void MomentOfInertia<MOMENT_CONSTANTS::PLANE_SYMMETRY>::get_angular_acceleration_body(const Eigen::Vector3d& angular_velocity,
     const Eigen::Vector3d& torque_body,
     double* angular_acceleration) const
{
    Eigen::Map< Eigen::Vector3d> omega_dot(angular_acceleration);
    Eigen::Matrix3d I = this->get_inertia_matrix();
    Eigen::Vector3d wxIw = angular_velocity.cross(I * angular_velocity);
    omega_dot = I.lu().solve(torque_body - wxIw); // see if LDLT could be used instead
    // TODO: replace with more optimum version
}

template<>
void MomentOfInertia<MOMENT_CONSTANTS::PRINCIPAL_AXIS>::get_angular_acceleration_body(const Eigen::Vector3d& angular_velocity,
    const Eigen::Vector3d& torque_body,
    double* angular_acceleration) const
{
    Eigen::Vector3d M_minus_wxIw = torque_body;
    M_minus_wxIw[0] -= (this->I[2] - this->I[1])*angular_velocity[2]*angular_velocity[1];
    M_minus_wxIw[1] -= (this->I[0] - this->I[2])*angular_velocity[2]*angular_velocity[0];
    M_minus_wxIw[2] -= (this->I[1] - this->I[0])*angular_velocity[1]*angular_velocity[0];

    angular_acceleration[0] = M_minus_wxIw[0] / this->I[0];
    angular_acceleration[1] = M_minus_wxIw[1] / this->I[1];
    angular_acceleration[2] = M_minus_wxIw[2] / this->I[2];
}

template<>
void MomentOfInertia<MOMENT_CONSTANTS::AXISYMMETRIC>::get_angular_acceleration_body(const Eigen::Vector3d& angular_velocity,
     const Eigen::Vector3d& torque_body,
     double* angular_acceleration) const
{
    double dIw = (this->I[0] - this->I[1])*angular_velocity[0];

    double inv_I = 1.0/this->I[1];
    angular_acceleration[0] = torque_body[0]/this->I[0];
    angular_acceleration[1] = (torque_body[1] + dIw*angular_velocity[2])*inv_I;
    angular_acceleration[2] = (torque_body[2] - dIw*angular_velocity[1])*inv_I;
}

template<>
void MomentOfInertia<MOMENT_CONSTANTS::EQUAL>::get_angular_acceleration_body(const Eigen::Vector3d& angular_velocity,
     const Eigen::Vector3d& torque_body,
     double* angular_acceleration) const
{
    double inv_I = 1.0/this->I[0];
    angular_acceleration[0] = torque_body[0]*inv_I;
    angular_acceleration[1] = torque_body[1]*inv_I;
    angular_acceleration[2] = torque_body[2]*inv_I;
}


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


void body::get_angular_acceleration(const Eigen::Vector3d& angular_velocity,
    const Eigen::Matrix3d& I,
    const Eigen::Vector3d& torque,
    double* angular_acceleration)
{

}
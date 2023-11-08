#include "Ephemeris.h"

#include "../util/functions.h"

constexpr double TWOPI = 6.283185307179586476925286766559;

double Ephemeris::trueAnomalyFromMeanAnomalyApprox(const double MA, const double eccentricity)
{
    double e2 = eccentricity * eccentricity;
    return MA + eccentricity * (2.0 - 0.25 * e2) * sin(MA) + 1.2 * e2 * sin(MA + MA) + 1.0833333333333333333 * e2 * eccentricity * sin(3.0 * MA);
}

double Ephemeris::meanAnomalyFromTrueAnomalyApprox(const double TA, const double eccentricity)
{
    double e2 = eccentricity * eccentricity;
    return TA - 2.0 * eccentricity * sin(TA) + e2 * (
        0.125 * (6.0 + e2) * sin(TA + TA) - eccentricity * 0.33333333333333333 * sin(3.0 * TA) + 0.15625 * e2 * sin(4.0 * TA));
}

double Ephemeris::eccentricAnomalyFromMeanAnomaly(const double MA, const double eccentricity)
{
    double EA = MA;
    constexpr int MAX_ITER = 20;
    constexpr double ACCURACY = 1e-15;
    // Halley's method
    for (int iter = 0; iter < MAX_ITER; iter++)
    {
        double sA = sin(EA);
        double cA = cos(EA); // hopefully optimizes to use fsincos
        double f = EA - eccentricity*sA - MA;
        double f_prime = 1.0 + eccentricity*cA;
        double dEA = f*f_prime/(f_prime*f_prime + 0.5*f*(eccentricity*sA)); // f'' = -e*sin(EA)
        EA -= dEA;
        if (fabs(dEA) < ACCURACY)
        {
            break;
        }
    }
    return EA;
}

double Ephemeris::meanAnomalyFromEccentricAnomaly(const double EA, const double eccentricity)
{
    return EA + eccentricity * sin(EA);
}

double Ephemeris::eccentricAnomalyFromTrueAnomaly(const double TA, const double eccentricity)
{
    return atan(sqrt(1.0 - eccentricity*eccentricity)*sin(TA) / (1.0 + eccentricity * cos(TA)));
}

double Ephemeris::meanAnomalyFromTrueAnomaly(const double TA, const double eccentricity)
{
    double num = sqrt(1.0 - eccentricity * eccentricity) * sin(TA);
    double cE = cos(TA);
    return atan2(num, (eccentricity + cE)) - eccentricity * num / (1.0 + eccentricity * cE);
}

double Ephemeris::trueAnomalyFromEccentricAnomaly(const double EA, const double eccentricity)
{
    double beta = eccentricity / (1.0 + sqrt(1.0 - eccentricity * eccentricity));
    return EA - 2.0 * atan(beta * sin(EA) / (1.0 + beta * cos(EA)));
}

double Ephemeris::trueAnomalyFromMeanAnomaly(const double MA, const double eccentricity)
{
    return trueAnomalyFromEccentricAnomaly(eccentricAnomalyFromMeanAnomaly(MA, eccentricity), eccentricity);
}

std::array<double, 6> Ephemeris::elements2cartesian(const std::array<double, 6>& elements, const double MU)
{
    double tmp = 1.0 - elements[1]*elements[1];
    double st = sin(elements[5]);
    double ct = cos(elements[5]);
    double tmp2 = 1.0 + elements[1] * ct;
    double radius = elements[0] * tmp / tmp2;
    double x = radius*ct;
    double y = radius*st;
    tmp = sqrt(MU*elements[0]*tmp)/(radius * tmp2);
    double v_x = -st * tmp;
    double v_y = (elements[1] + ct) * tmp;

    if (fabs(elements[2]) < 1e-8)
    {
        std::array<double, 6> out = { x,y,0.0,v_x,v_y,0.0 };
        return out;
    }

    double sw = sin(elements[4]);
    double cw = cos(elements[4]);
    double so = sin(elements[3]);
    double co = cos(elements[3]);

    st = sin(elements[2]);
    ct = sqrt(1.0 - st*st);
    double Rxx = cw * co - sw * ct * so;
    double Rxy = -(sw * co + cw * ct * so);
    double Ryx = cw * so + sw * ct * co;
    double Ryy = cw * ct * co - sw * so;
    double Rzx = sw * st;
    double Rzy = cw * st;

    std::array<double, 6> out;
    out[0] = Rxx*x + Rxy*y;
    out[1] = Ryx*x + Ryy*y;
    out[2] = Rzx*x + Rzy*y;
    out[3] = Rxx*v_x + Rxy*v_y;
    out[4] = Ryx*v_x + Ryy*v_y;
    out[5] = Rzx*v_x + Rzy*v_y;
    return out;
}

std::array<double, 6> Ephemeris::cartesian2elements(const std::array<double, 6>& cartesian, const double MU)
{
    std::array<double, 6> oe{};
    double h[3] = {cartesian[1]*cartesian[5] - cartesian[2]*cartesian[4],
                   cartesian[2]*cartesian[3] - cartesian[0]*cartesian[5], 
                   cartesian[0]*cartesian[4] - cartesian[1]*cartesian[3] };

    // std::array<double,2> n = {h[1],-h[0]}; // z is implicit 0
    double v2 = func::dot3(&cartesian[3], &cartesian[3]);
    double r_inv = 1.0 / func::mag3(cartesian.data());
    double rv = func::dot3(cartesian.data(), &cartesian[3]);
    double e[3];
    const double invMU = 1.0 / MU;
    double tmp1 = v2*invMU - r_inv;
    double tmp2 = rv*invMU;
    e[0] = cartesian[0]*tmp1 + cartesian[3]*tmp2;
    e[1] = cartesian[1]*tmp1 + cartesian[4]*tmp2;
    e[2] = cartesian[2]*tmp1 + cartesian[5]*tmp2;
    double egy = v2*0.5 - MU*r_inv;

    oe[0] = -MU / (2*egy);
    oe[1] = func::mag3(e);
    double inv_e = 1 / oe[1];
    double nmag = h[0]*h[0] + h[1]*h[1];
    oe[2] = acos(h[2] / sqrt(nmag + h[2]*h[2]));

    if (fabs(oe[2]) > 1e-9) {
        nmag = 1.0 / sqrt(nmag);
        oe[3] = acos(h[1]*nmag);
        if (h[0] > 0) 
        {
            oe[3] = TWOPI - oe[3];
        }

        oe[4] = acos((h[1]*e[0] - h[0]*e[1])*nmag*inv_e);
        if (e[2] < 0) 
        {
            oe[4] = TWOPI - oe[4];
        }
    }

    oe[5] = acos(func::dot3(e,cartesian.data())*r_inv*inv_e);
    if (rv < 0) 
    {
        oe[5] = TWOPI - oe[5];
    }
    return oe;
}

void Ephemeris::kepler2position(const double* oe, double* pos)
{
    double st = sin(oe[5]);
    double ct = cos(oe[5]);
    double radius = oe[0]*(1.0 - oe[1]*oe[1]) / (1.0 + oe[1]*ct);
    double x = radius*ct;
    double y = radius*st;

    if (fabs(oe[2]) < 1e-8) 
    {
        pos[0] = x;
        pos[1] = y;
        pos[0] = 0.0;
        return;
    }

    double sw = sin(oe[4]);
    double cw = cos(oe[4]);
    double so = sin(oe[3]);
    double co = cos(oe[3]);

    st = sin(oe[2]);
    ct = sqrt(1.0 - st*st);
    double Rxx = cw * co - sw * ct * so;
    double Rxy = -(sw * co + cw * ct * so);
    double Ryx = cw * so + sw * ct * co;
    double Ryy = cw * ct * co - sw * so;
    double Rzx = sw * st;
    double Rzy = cw * st;

    pos[0] = Rxx*x + Rxy*y;
    pos[1] = Ryx*x + Ryy*y;
    pos[2] = Rzx*x + Rzy*y;
}


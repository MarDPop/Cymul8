#include "Geometry.h"

constexpr double a1 = 42697.67270715754; //a1 = a*e2
constexpr double a2 = 1.8230912546075456E9; //a2 = a1*a1
constexpr double a3 = 142.91722289812412; //a3 = a1*e2/2
constexpr double a4 = 4.557728136518864E9; //a4 = 2.5*a2
constexpr double a5 = 42840.589930055656; //a5 = a1+a3
constexpr double a6 = 0.9933056200098622; //a6 = 1-e2
constexpr double e2 = 0.006694380066765; //WGS-84 first eccentricity squared
constexpr double EARTH_EQUATOR_R = 6378137; // m
constexpr double EARTH_POLAR_R = 6356752.314; // m

Coordinate::Geodetic WGS84::ECEF2LLA(const Coordinate::GeocentricFixed& ecef)
{
    Coordinate::Geodetic lla;   //Results go here (Lat, Lon, Altitude)
    double zp = fabs(ecef.z);
    double z2 = zp * zp;
    double w2 = ecef.x * ecef.x + ecef.y * ecef.y;
    double w = sqrt(w2);
    double r2 = 1.0 / (w2 + z2);
    double r_inv = sqrt(r2);
    lla.longitude = atan2(ecef.y, ecef.x);       //Lon (final)
    double s2 = z2 * r2;
    double c2 = w2 * r2;
    double u = a2 * r_inv;
    double v = a3 - a4 * r_inv;
    double s, c, ss;
    if (c2 > 0.3) {
        s = (zp * r_inv) * (1.0 + c2 * (a1 + u + s2 * v) * r_inv);
        lla.latitude = asin(s);
        ss = s * s;
        c = sqrt(1.0 - ss);
    }
    else {
        c = (w * r_inv) * (1.0 - s2 * (a5 - u - c2 * v) * r_inv);
        lla.latitude = acos(c);
        ss = 1.0 - c * c;
        s = sqrt(ss);
    }
    double g = 1.0 - e2 * ss;
    double rg = EARTH_EQUATOR_R / sqrt(g);
    double rf = a6 * rg;
    u = w - rg * c;
    v = zp - rf * s;
    double f = c * u + s * v;
    double m = c * v - s * u;
    double p = m / (rf / g + f);
    lla.latitude += p;
    lla.altitude = f + m * p * 0.5;
    if (ecef.z < 0) {
        lla.latitude = -lla.latitude;
    }
    return lla; 
}

Coordinate::GeocentricFixed WGS84::LLA2ECEF(const Coordinate::Geodetic& lla)
{

    double s = sin(lla.latitude);
    double n = EARTH_EQUATOR_R / sqrt(1.0 - e2 * s * s);
    double tmp = (n + lla.altitude) * sqrt(1 - s * s);
    Coordinate::GeocentricFixed ecef;
    ecef.x = tmp * cos(lla.longitude);
    ecef.y = tmp * sin(lla.longitude);
    ecef.z = s * (n * (1.0 - e2) + lla.altitude);
    return ecef;
}

double WGS84::vincentyFormulae(double long1, double lat1, double long2, double lat2) {
    constexpr double EARTH_FLATTENING = 0.003352810664747;
    constexpr double EARTH_FLATTENING_1 = 1.0 - EARTH_FLATTENING;
    constexpr double EARTH_EQUATOR_R = 6378137; // m

    double U1 = atan(EARTH_FLATTENING_1 * tan(lat1));
    double U2 = atan(EARTH_FLATTENING_1 * tan(lat2));
    double L = long2 - long1;
    double l = L;
    double calpha;
    double st;
    double ct;
    double d;
    double a;
    double b;
    calpha = st = ct = d = a = b = 0;
    for (int i = 0; i < 20; i++) {
        double su2 = sin(U2);
        double su1 = sin(U1);
        double cu2 = cos(U2);
        double cu1 = cos(U1);
        a = cu2 * sin(l);
        b = cu1 * su2 - su1 * cu2 * cos(l);
        st = sqrt(a * a + b * b);
        ct = su2 * su1 + cu1 * cu2 * cos(l);
        a = atan2(st, ct);
        b = cu2 * cu1 * sin(l) / st;
        calpha = 1 - b * b;
        d = ct - 2 * su1 * su2 / calpha;
        double C = EARTH_FLATTENING / 16 * calpha * (4 + EARTH_FLATTENING * (4 - 3 * calpha));
        l = L + (1 - C) * EARTH_FLATTENING * b * (a + C * st * (d + C * ct * (2 * d - 1)));
    }
    double u2 = calpha * (EARTH_EQUATOR_R * EARTH_EQUATOR_R / (EARTH_POLAR_R * EARTH_POLAR_R) - 1);
    double A = 1 - u2 / 16384 * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)));
    double B = u2 / 1024 * (256 + u2 * (-128 + u2 * (74 - 47 * u2)));
    L = B * st * (d + 0.25 * B * (ct * (2 * d - 1)) - 0.1666 * B * d * (-3 + 4 * b * b) * (-3 + 4 * d));
    return EARTH_POLAR_R * A * (a - L);
}

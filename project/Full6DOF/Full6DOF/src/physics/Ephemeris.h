#pragma once

#pragma warning(push)
#pragma warning(disable : 26495)

#include "EpochTime.h"

#include <array>
#include <vector>
#include <cmath>

struct Ephemeris
{
    union 
    {
        std::array<double,6> elements;
        struct 
        {
            double semi_major_axis;
            double eccentricity;
            double inclination;
            double longitude_of_ascending_node;
            double argument_of_periapse;
            double true_anomaly;
        };
    };

    double MU; // km3 / s2

    Ephemeris() {}

    Ephemeris(  const std::array<double,6>& _elements, 
                double _MU) : 
                    elements(_elements), 
                    MU(_MU)
                    {}

    static double trueAnomalyFromMeanAnomalyApprox(const double MA, const double eccentricity);

    static double meanAnomalyFromTrueAnomalyApprox(const double TA, const double eccentricity);

    static double eccentricAnomalyFromMeanAnomaly(const double MA, const double eccentricity);

    static double meanAnomalyFromEccentricAnomaly(const double EA, const double eccentricity);

    static double eccentricAnomalyFromTrueAnomaly(const double TA, const double eccentricity);

    static double meanAnomalyFromTrueAnomaly(const double TA, const double eccentricity);

    static double trueAnomalyFromEccentricAnomaly(const double EA, const double eccentricity);

    static double trueAnomalyFromMeanAnomaly(const double MA, const double eccentricity);

    static std::array<double, 6> elements2cartesian(const std::array<double, 6>& elements, const double MU);

    static std::array<double, 6> cartesian2elements(const std::array<double, 6>& cartesian, const double MU);

    static void kepler2position(const double* oe, double* pos);

};
#pragma warning(pop)
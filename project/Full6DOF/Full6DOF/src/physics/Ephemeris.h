#pragma once

#include "Time.h"
#include "../util/fast_math.h"
#include "../util/Table.h"

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

    const double MU;

    Ephemeris(  const std::array<double,6>& _elements, 
                double _MU) : 
                    elements(_elements), 
                    MU(_MU)
                    {}

    Ephemeris(double _MU) : MU(_MU) {}

    static double trueAnomalyFromEccentricAnomaly(const double EA, const double eccentricity) 
    {
        double beta = eccentricity/(1.0 + sqrt(1.0 - eccentricity*eccentricity));
        return EA - 2.0*atan(beta * sin(EA)/(1.0 + beta*cos(EA));
    }

    static double trueAnomalyFromMeanAnomaly(const double MA, const double eccentricity) 
    {
        return trueAnomalyFromEccentricAnomaly(eccentricAnomalyFromMeanAnomaly(MA, eccentricity), eccentricity);
    }

    static double trueAnomalyFromMeanAnomalyApprox(const double MA, const double eccentricity) 
    {
        
    }

    static double meanAnomalyFromTrueAnomalyApprox(const double TA, const double eccentricity) 
    {
        
    }

    static double eccentricAnomalyFromMeanAnomaly(const double MA, const double eccentricity) 
    {
           
    }

    static double meanAnomalyFromEccentricAnomaly(const double EA, const double eccentricity) 
    {
        return EA + eccentricity*sin(EA);
    }

    static double eccentricAnomalyFromTrueAnomaly(const double TA, const double eccentricity) 
    {
        return atan((sqrt(1.0 - eccentricity*eccentricity)*sin(TA)/(1.0 + eccentricity*cos(TA));
    }

    static double meanAnomalyFromTrueAnomaly(const double TA, const double eccentricity) 
    {
        double num = sqrt(1.0 - eccentricity*eccentricity)*sin(TA);
        double cE = cos(TA);
        return atan2(num, (eccentricity + cE)) - eccentricity*num/(1.0 + eccentricity*cE);
    }

    static std::array<double,6> elements2cartesian(const std::array<double,6>& elements) {}

    static std::array<double,6> cartesian2elements(const std::array<double,6>& cartesian) {}

};


class EphemerisHistory
{
    double _MU = 1.0;

    std::vector<std::array<double, 6>> _ephemeris;

    std::vector<std::array<double, 6>> _dephemeris;

    std::vector<double> _mjd;

public:

    void load(std::string file);

    Ephemeris get_at_time(double mjd);
};
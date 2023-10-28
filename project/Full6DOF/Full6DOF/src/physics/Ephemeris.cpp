#include "Ephemeris.h"
#include <fstream>
#include "../util/StringHelper.h"

bool file_is_horizons(std::string filename)
{
    int idx_ext = filename.size() - 1;
    int idx = 3;
    const char* horizons_ext = ".JPL";
    while (idx > 0)
    {
        if (filename[idx_ext] != horizons_ext[idx])
        {
            return false;
        }
        idx_ext--;
        idx--;
    }
    return true;
}

void EphemerisHistory::load(std::string filename)
{
    std::ifstream file(filename);
    if (file.is_open())
    {
        if (file_is_horizons(filename))
        {
            // use CSV output
        }
        else
        {
            int lineNo = 0;
            for (std::string line; std::getline(file, line); )
            {
                if (lineNo == 0)
                {
                    _MU = std::stod(line);
                    lineNo = 1;
                    continue;
                }
                auto data = strings::split(line);

            }
        }
    }
}

Ephemeris EphemerisHistory::get_at_time(double mjd)
{
    Ephemeris interpolated(_MU);

    bool before = mjd < _mjd[0];
    bool after = mjd > _mjd.back();
    if (before || after)
    {
        auto idx = after * (_mjd.size() - 1);
        interpolated.elements = _ephemeris[idx];
        double dt = (mjd - _mjd[idx]) * Time::JULIAN_DAY_SEC_D;
        double mean_motion = sqrt(_MU / fmath::CB(interpolated.semi_major_axis));
        double MA = interpolated.elements.back() + mean_motion * dt;
        interpolated.true_anomaly = Ephemeris::trueAnomalyFromMeanAnomaly(MA, interpolated.eccentricity);
    }
    else
    {
        auto it = std::lower_bound(_mjd.begin(), _mjd.end(), mjd);
        auto idx = std::distance(_mjd.begin(), it);
        double delta = mjd - *it;

        interpolated.elements = _ephemeris[idx];
        const auto& d = _dephemeris[idx];
        for (unsigned i = 0; i < 6; i++)
        {
            interpolated.elements[i] += d[i] * delta;
        }

        interpolated.true_anomaly = Ephemeris::trueAnomalyFromMeanAnomaly(interpolated.true_anomaly, interpolated.eccentricity);
    }

    return interpolated;
}
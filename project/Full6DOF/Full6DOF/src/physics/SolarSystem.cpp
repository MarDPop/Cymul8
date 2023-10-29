#include "SolarSystem.h"

#include <fstream>
#include "Time.h"
#include "../util/StringHelper.h"
#include "../util/functions.h"

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
        constexpr double DEG = 1.74532925199432957692369e-2;
        if (file_is_horizons(filename))
        {
            // use CSV output
            bool startEphem = false;
            for (std::string line; std::getline(file, line); )
            {
                if (line[0] == 'K' && line[1] == 'e')
                {
                    int idx_start = 18;
                    int idx_end = 19;
                    while (line[idx_end] != ' ' && idx_end < line.size())
                    {
                        idx_end++;
                    }
                    _current.MU = std::stod(line.substr(idx_start, idx_end)); // km3 / s2
                }
                if (!startEphem)
                {
                    startEphem = line[0] == '$';
                    continue;
                }
                if (line[0] == '$')
                {
                    break;
                }
                auto data = strings::split(line, ',');
                auto mjd = std::stod(data[0]) - EpochTime::MJD_JULIAN_DATE;
                _mjd.push_back(mjd);

                std::array<double, 6> e;
                e[1] = std::stod(data[2]);
                e[2] = std::stod(data[4]) * DEG;
                e[3] = std::stod(data[5]) * DEG;
                e[4] = std::stod(data[6]) * DEG;
                e[5] = std::stod(data[9]) * DEG; // mean anomaly
                e[0] = std::stod(data[11]); // km
            }
        }
        else
        {
            int lineNo = 0;
            for (std::string line; std::getline(file, line); )
            {
                if (lineNo == 0)
                {
                    _current.MU = std::stod(line);
                    lineNo = 1;
                    continue;
                }
                auto data = strings::split(line);
                _mjd.push_back(std::stod(data[0]));
                std::array<double, 6> e;
                for (auto i = 0; i < 6; i++)
                {
                    e[i] = std::stod(data[i + 1]);
                }
                _ephemeris.push_back(e);
            }
        }

        for (auto nE = 1u; nE < _ephemeris.size(); nE++)
        {
            std::array<double, 6> de;
            const auto& lo = _ephemeris[nE - 1];
            const auto& hi = _ephemeris[nE];
            auto dmjd = 1.0 / (_mjd[nE] - _mjd[nE - 1]);
            for (auto i = 0; i < 6; i++)
            {
                de[i] = (hi[i] - lo[i]) * dmjd;
            }
            _dephemeris.push_back(de);
        }
    }
}

void EphemerisHistory::set(double mjd)
{
    bool before = mjd < _mjd[0];
    bool after = mjd > _mjd.back();
    if (before || after)
    {
        auto idx = after * (_mjd.size() - 1);
        _current.elements = _ephemeris[idx];
        double dt = (mjd - _mjd[idx]) * EpochTime::JULIAN_DAY_SEC_D;
        double mean_motion = sqrt(_current.MU / func::CB(_current.semi_major_axis));
        double MA = _current.elements.back() + mean_motion * dt;
        _current.true_anomaly = Ephemeris::trueAnomalyFromMeanAnomaly(MA, _current.eccentricity);
    }
    else
    {
        auto it = std::lower_bound(_mjd.begin(), _mjd.end(), mjd);
        auto idx = std::distance(_mjd.begin(), it);
        double delta = mjd - *it;

        _current.elements = _ephemeris[idx];
        const auto& d = _dephemeris[idx];
        for (unsigned i = 0; i < 6; i++)
        {
            _current.elements[i] += d[i] * delta;
        }

        _current.true_anomaly = Ephemeris::trueAnomalyFromMeanAnomaly(_current.true_anomaly, _current.eccentricity);
    }
}

void OrientationHistory_Constant::set(double mjd)
{
    const double angle = _rotation_rate*(mjd - _mjd_ref);
    
    _pci2pcef = Eigen::AngleAxisd(angle, _axis_ref.row(2));
}
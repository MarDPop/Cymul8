#pragma once

#include <vector>
#include <array>

typedef long UNIX_TIMESTAMP;

typedef long GPS_TIMESTAMP;

class EpochTime
{
    static constexpr std::array<int, 28> LEAP_SECONDS_UTC = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27};

    // days in which leap second is added at 23:59:60
    static constexpr std::array<int,28> LEAP_SECOND_MJD = 
    {
        40587,
        41498,
        41682,
        42047,
        42412,
        42777,
        43143,
        43508,
        43873,
        44238,
        44785,
        45150,
        45515,
        46246,
        47160,
        47891,
        48256,
        48803,
        49168,
        49533,
        50082,
        50629,
        51178,
        53735,
        54831,
        56108,
        57203,
        57753
    };

    double _day_sec; 

    int _mjdn;

public:

    static constexpr long SECONDS2NANOSECONDS = 1000000000L;

    static constexpr double JULIAN_DATE_J2000 = 2451545.0;

    static constexpr unsigned JDN_J2000 = 2451545;

    static constexpr unsigned JULIAN_DAY_SEC = 86400;

    static constexpr unsigned JULIAN_DAY_SEC_D = 86400.0;

    static constexpr long JULIAN_DAY_NANOSEC = JULIAN_DAY_SEC*SECONDS2NANOSECONDS;

    static constexpr double SEC_2_DAY_FRACTION = 1.0/86400.0;

    static constexpr double NANOSEC_2_DAY_FRACTION = 1.0/JULIAN_DAY_NANOSEC;

    static constexpr double MJD_JULIAN_DATE = 2400000.5;

    static constexpr double JULIAN_DATE_UNIX_EPOCH = 2440587.5;

    static constexpr int MJD_UNIX_EPOCH = 40587;

    static constexpr int UNIX_TAI_LEAPSECONDS = 10;

    static constexpr double HALF_DAY_SEC = 43200.0;

    EpochTime() : _mjdn(JDN_J2000), _day_sec(HALF_DAY_SEC) {}

    EpochTime(int mjdn, double day_sec) : _mjdn(mjdn), _day_sec(day_sec) {}

    EpochTime(UNIX_TIMESTAMP timestamp_ns);

    int get_MJDN() const
    {
        return this->_mjdn;
    }

    double get_seconds_past_midnight() const
    {
        return this->_day_sec;
    }

    void operator+=(double sec)
    {
        this->_day_sec += sec;
        if(this->_day_sec > JULIAN_DAY_SEC_D) 
        {
            double days = this->_day_sec*SEC_2_DAY_FRACTION;
            int nDays = static_cast<int>(days);

            this->_mjdn += nDays;
            this->_day_sec = (days - nDays)*JULIAN_DAY_SEC_D;
        }
    }

    static int get_UTC_leap_seconds(int mjdn); 

    double to_modified_julian_date() const
    {
        return this->_mjdn + this->_day_sec * SEC_2_DAY_FRACTION;
    }

    double to_julian_date() const
    {
        return static_cast<double>(this->_mjdn) + (this->_day_sec*SEC_2_DAY_FRACTION + MJD_JULIAN_DATE);
    }

    double days_past_J2000() const
    {
        return static_cast<double>(this->_mjdn - JDN_J2000) + (this->_day_sec * SEC_2_DAY_FRACTION - 0.5);
    }

    static double get_besselian_years(double jd)
    {
        return 1900.0 + (jd - 2415020.31352) / 365.242198781;
    }

};
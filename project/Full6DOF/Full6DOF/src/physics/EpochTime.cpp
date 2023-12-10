#include "EpochTime.h"

namespace Time
{

    EpochDate::EpochDate(UNIX_TIMESTAMP timestamp_ns) :
        _epoch(EPOCH::UNIX)
    {
        long days = timestamp_ns / JULIAN_DAY_NANOSEC;
        long day_nanos = timestamp_ns - (days * JULIAN_DAY_NANOSEC);

        _day_number = static_cast<int>(days);
        int leap_seconds_after_unix = get_UTC_leap_seconds(MJD_UNIX_EPOCH + _day_number); // TODO: see if need to subtract UNIX_TAI

        day_nanos -= leap_seconds_after_unix;
        if (day_nanos < 0)
        {
            day_nanos += JULIAN_DAY_NANOSEC;
            _day_number--;
        }

        this->_day_fraction = day_nanos * NANOSEC_2_DAY_FRACTION;
    }

    EpochDate EpochDate::from_unix_timestamp(UNIX_TIMESTAMP timestamp_ns, EPOCH epoch)
    {
        EpochDate date;

        return date;
    }

    int LeapSeconds::get_UTC_leap_seconds_MJD(const int mjdn)
    {
        auto idx = LEAP_SECOND_MJD.size();
        while (--idx > 0)
        {
            if (mjdn > LEAP_SECOND_MJD[idx])
            {
                break;
            }
        }
        return LEAP_SECONDS_UTC[idx];
    }

    uint8_t LeapSeconds::get_UTC_leap_seconds(const unsigned short unix_day)
    {
        auto idx = LEAP_SECOND_UNIX.size();
        while (--idx > 0)
        {
            if (unix_day > LEAP_SECOND_UNIX[idx])
            {
                break;
            }
        }
        return LEAP_SECONDS_UTC[idx];
    }

    bool LeapSeconds::is_leap_day_MJD(int mjdn)
    {
        auto idx = LEAP_SECOND_MJD.size();
        while (--idx > 0)
        {
            if (mjdn < LEAP_SECOND_MJD[idx])
            {
                break;
            }
            if (mjdn == LEAP_SECOND_MJD[idx])
            {
                return true;
            }
        }
        return false;
    }

    bool LeapSeconds::is_leap_day(const unsigned short unix_day)
    {
        auto idx = LEAP_SECOND_UNIX.size() - 1;
        while (unix_day > LEAP_SECOND_UNIX[idx] && --idx > 0){}

        return LEAP_SECOND_UNIX[idx] == unix_day;
    }

    double getGMST(double julianUT1) 
    {
        //http://aa.usno.navy.mil/faq/docs/GAST.php
        int JD0 = (int)julianUT1;
        double H = (julianUT1 - JD0) * 24;
        double T = (julianUT1 - JDN_J2000) / 36525.0;
        return 6.697374558 + 0.06570982441908 * (JD0 - JDN_J2000) + 1.00273790935 * H + 0.000026 * T * T;
    }

    double getGAST(double julianUT1) 
    {
        constexpr double DEG2RAD = 0.01745329251994329576923690768;
        double D = julianUT1 - JDN_J2000;
        double o = 125.04 - 0.052954 * D;
        double L = 280.47 + 0.98565 * D;
        double e = 23.4393 - 0.0000004 * D;
        double d = -0.000319 * sin(o* DEG2RAD) - 0.000024 * sin(DEG2RAD*(2 * L));
        double eqeq = d * cos(DEG2RAD*e);
        return getGMST(julianUT1) + eqeq;
    }

    std::array<int, 3> jdn2ymd(int jdn)
    {
        // All of these must be integer divisions
        std::array<int, 3> ymd;
        int f = jdn + 1401 + (((4 * jdn + 274277) / 146097) * 3) / 4 - 38;
        int e = 4 * f + 3;
        int g = (e % 1461) / 4;
        int h = 5 * g + 2;
        ymd[2] = (h % 153) / 5 + 1;
        ymd[1] = (h / 153 + 2) % 12 + 1;
        ymd[0] = (e / 1461) - 4716 + (12 + 2 - ymd[1]) / 12;
        return ymd;
    }

    int ymd2jdn(int y, int m, int d) 
    {
        // All of these must be integer divisions
        int M2 = (m - 14) / 12;
        return (1461 * (y + 4800 + M2)) / 4 + (367 * (m - 2 - 12 * M2)) / 4 - (3 * (y + 4900 + M2 / 100)) / 4 + d - 32075;
    }

    UNIX_TIMESTAMP to_unix_timestamp(Gregorian date)
    {
        std::tm timeinfo = {};
        timeinfo.tm_year = date.year; // Years since 1900
        timeinfo.tm_mon = date.month;    // Months since January (0-based)
        timeinfo.tm_mday = date.day;         // Day of the month
        timeinfo.tm_hour = date.hour;
        timeinfo.tm_min = date.min;
        timeinfo.tm_sec = date.sec;

        std::time_t time = std::mktime(&timeinfo);
        auto time_point = std::chrono::system_clock::from_time_t(time);

        // Get the duration since the epoch
        auto ms_since_epoch = std::chrono::duration_cast<std::chrono::milliseconds>(time_point.time_since_epoch()).count();

        ms_since_epoch += date.milliseconds;

        return ms_since_epoch * 1000000L;
    }

    EpochDate to_epoch_date_j2000(UNIX_TIMESTAMP ts)
    {
        int days_unix = static_cast<int>(ts / JULIAN_DAY_NANOSEC);

        ts -= LeapSeconds::get_UTC_leap_seconds(days_unix) * Time::SECONDS2NANOSECONDS;

        days_unix += ts < 0;

        long nanos_past_mid = ts - static_cast<long>(days_unix) * JULIAN_DAY_NANOSEC;

        return EpochDate(days_unix, nanos_past_mid * NANOSEC_2_DAY_FRACTION, EPOCH::J2000_UTC);
    }
}
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
}
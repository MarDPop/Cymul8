#include "EpochTime.h"

EpochTime::EpochTime(UNIX_TIMESTAMP timestamp_ns)
{
    long days = timestamp_ns / JULIAN_DAY_NANOSEC;
    long day_nanos = timestamp_ns - (days * JULIAN_DAY_NANOSEC);

    int leap_seconds_after_unix = get_UTC_leap_seconds(MJD_UNIX_EPOCH + static_cast<int>(days)); // TODO: see if need to subtract UNIX_TAI

    day_nanos -= leap_seconds_after_unix;
    if (day_nanos < 0)
    {
        day_nanos += JULIAN_DAY_NANOSEC;
        days--;
    }

    this->_mjdn = days + MJD_UNIX_EPOCH;
    this->_day_sec = day_nanos * NANOSEC_2_DAY_FRACTION;
}

int EpochTime::get_UTC_leap_seconds(int mjdn)
{
    unsigned _leap_idx = 0;
    if (mjdn > LEAP_SECOND_MJD.back())
    {
        return LEAP_SECONDS_UTC.back();
    }

    while (_leap_idx < (LEAP_SECOND_MJD.size() - 1) && mjdn > LEAP_SECOND_MJD[_leap_idx + 1])
    {
        _leap_idx++;
    }

    while (_leap_idx > 0 && mjdn < LEAP_SECOND_MJD[_leap_idx])
    {
        _leap_idx--;
    }

    return LEAP_SECONDS_UTC[_leap_idx];
}
#pragma once

#include <vector>
#include <array>
#include <string>
#include <ctime>
#include <chrono>

#if __cplusplus >= 202002L


#endif


namespace Time
{

    typedef long UNIX_TIMESTAMP;

    typedef long GPS_TIMESTAMP;

    constexpr long SECONDS2NANOSECONDS = 1'000'000'000L;

    constexpr double JULIAN_DATE_J2000_TT = 2451545.0;

    // January 1, 2000, 11:58:55.816 UTC. = 2451545.0 TT
    constexpr UNIX_TIMESTAMP J2000_UNIX_NS = 946'727'935'816'000'000L;

    constexpr int TAI_UTC = 10;

    constexpr double TAI_TT = 32.184;

    constexpr long TAI_TT_NS = 32'184'000'000L;

    constexpr int LEAP_SECONDS_UTC_J2000 = 22;

    constexpr int JDN_J2000 = 2451545;

    constexpr int JULIAN_DAY_SEC = 86400;

    constexpr double JULIAN_DAY_SEC_D = 86400.0;

    constexpr long JULIAN_DAY_NANOSEC = static_cast<long>(86400L*SECONDS2NANOSECONDS);

    constexpr double SEC_2_DAY_FRACTION = 1.0 / 86400.0;

    constexpr double NANOSEC_2_DAY_FRACTION = 1.0 / JULIAN_DAY_NANOSEC;

    constexpr double MJD_JULIAN_DATE = 2400000.5;

    constexpr double JULIAN_DATE_UNIX_EPOCH = 2440587.5;

    constexpr int MJD_UNIX_EPOCH = 40587;

    constexpr int MJD_J2000_MIDNIGHT = 51544;

    constexpr int UNIX_TAI_LEAPSECONDS = 10;

    constexpr double HALF_DAY_SEC = 43200.0;

    struct Gregorian
    {
        unsigned short milliseconds = 0u;
        unsigned char year; // after 1900
        unsigned char month;
        unsigned char day;
        unsigned char hour;
        unsigned char min;
        unsigned char sec;
    };

    class LeapSeconds
    {

    public:

        static constexpr std::array<uint8_t, 28> LEAP_SECONDS_UTC = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27 };

        // days in which leap second is added at 23:59:60
        static constexpr std::array<int, 28> LEAP_SECOND_MJD =
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

        static constexpr std::array<unsigned short, 28> LEAP_SECOND_UNIX =
        {
            0,
            911,
            1095,
            1460,
            1825,
            2190,
            2556,
            2921,
            3286,
            3651,
            4198,
            4563,
            4928,
            5659,
            6573,
            7304,
            7669,
            8216,
            8581,
            8946,
            9495,
            10042,
            10591,
            13148,
            14244,
            15521,
            16616,
            17166
        };

        static int get_UTC_leap_seconds_MJD(int mjdn);

        static uint8_t get_UTC_leap_seconds(unsigned short unix_day);

        static bool is_leap_day_MJD(int mjdn);

        static bool is_leap_day(unsigned short unix_day);
    };

    enum class EPOCH : int
    {
        JULIAN,
        MODIFIED_JULIAN,
        J2000_UTC,
        GPS,
        UNIX,
        TAI,
        TT
    };

    class EpochDate
    {
        double _day_fraction;

        int _day_number;

        EPOCH _epoch;

    public:

        EpochDate() :
            _day_fraction(0.5),
            _day_number(JDN_J2000),
            _epoch(EPOCH::JULIAN) {}

        EpochDate(int dn,
            double day_fraction,
            EPOCH epoch = EPOCH::JULIAN) :
            _day_fraction(day_fraction),
            _day_number(dn),
            _epoch(epoch) {}

        EpochDate(UNIX_TIMESTAMP timestamp_ns);

        static EpochDate from_unix_timestamp(UNIX_TIMESTAMP timestamp_ns, EPOCH epoch);

        int get_day_number() const
        {
            return this->_day_number;
        }

        double get_day_fraction() const
        {
            return this->_day_fraction;
        }

        double get_seconds_past_midnight(double length_of_day = JULIAN_DAY_SEC_D) const
        {
            return this->_day_fraction * length_of_day;
        }

        void operator+=(EpochDate date)
        {
            this->_day_fraction += date._day_fraction;
            this->_day_number += date._day_number;
            if (this->_day_fraction > 1.0)
            {
                this->_day_number++;
                this->_day_fraction--;
            }
        }

        static int get_UTC_leap_seconds(int mjdn);

        double days_past_j2000() const;

        UNIX_TIMESTAMP to_unix_timestamp() const;
    };

    class TimeTable
    {
        // first value at UNIX epoch (Jan 1, 1970 00:00 GMT), 1 per day (first 2.5 years are junk)

        std::vector<float> _dut1;

        std::vector<float> _LOD;

        std::vector<long> _UNIX_midnight;

        TimeTable(const char* filename);

    public:

        static constexpr const char * TABLE_FILENAME = "d";

        static TimeTable& instance()
        {
            static TimeTable table(TABLE_FILENAME);
            return table;
        }

        double get_dut1(double jd2000) const
        {
            unsigned idx = static_cast<unsigned>(jd2000);
            double delta = jd2000 - static_cast<double>(idx);
            return _dut1[idx] + (_dut1[idx + 1] - _dut1[idx]) * delta;
        }

        double get_dut1(const EpochDate& jd2000) const
        {
            unsigned idx = jd2000.get_day_number();
            return _dut1[idx] + (_dut1[idx + 1] - _dut1[idx])*jd2000.get_day_fraction();
        }

        float get_dut1(unsigned jd2000)
        {
            return _dut1[jd2000];
        }

        float get_LOD(unsigned jd2000)
        {
            return _LOD[jd2000];
        }

        double get_day_fraction_ratio(unsigned jd2000)
        {
            return 1.0/(JULIAN_DAY_SEC_D + _LOD[jd2000]);
        }

        long get_UNIX_midnight(unsigned jd2000)
        {
            return _UNIX_midnight[jd2000];
        }
    };

    inline double get_besselian_years(double jd)
    {
        return 1900.0 + (jd - 2415020.31352) / 365.242198781;
    }

    double getGMST(double julianUT1);

    double getGAST(double julianUT1);

    UNIX_TIMESTAMP to_unix_timestamp(Gregorian date);

    std::array<int, 3> jdn2ymd(int jdn);

    int ymd2jdn(int y, int m, int d);

    EpochDate to_epoch_date_j2000(UNIX_TIMESTAMP ts);


}
#include <cstdint>

namespace fmath
{
    float inv_fast(float x) {
        union { float f; int i; } v;
        float w, sx;
        int m;

        sx = (x < 0) ? -1 : 1;
        x = sx * x;

        v.i = (int)(0x7EF127EA - *(uint32_t*)&x);
        w = x * v.f;

        // Efficient Iterative Approximation Improvement in horner polynomial form.
        v.f = v.f * (2 - w);     // Single iteration, Err = -3.36e-3 * 2^(-flr(log2(x)))
        // v.f = v.f * ( 4 + w * (-6 + w * (4 - w)));  // Second iteration, Err = -1.13e-5 * 2^(-flr(log2(x)))
        // v.f = v.f * (8 + w * (-28 + w * (56 + w * (-70 + w *(56 + w * (-28 + w * (8 - w)))))));  // Third Iteration, Err = +-6.8e-8 *  2^(-flr(log2(x)))

        return v.f * sx;
    }

    inline double constexpr divide(double y, double x) 
    {
        // calculates y/x
        union {
            double dbl;
            unsigned long long ull;
        } u;
        u.dbl = x;                      // x = x
        u.ull = (0xbfcdd6a18f6a6f52ULL - u.ull) >> (unsigned char)1;
        // pow( x, -0.5 )
        u.dbl *= u.dbl;                 // pow( pow(x,-0.5), 2 ) = pow( x, -1 ) = 1.0/x
        return u.dbl * y;               // (1.0/x) * y = y/x
    }
}
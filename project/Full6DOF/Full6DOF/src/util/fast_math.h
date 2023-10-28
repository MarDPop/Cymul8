#pragma once

#include <cstdint>

namespace fmath
{
    template<typename T>
    T SQ(T x)
    {
        return x*x;
    }

    template<typename T>
    T CB(T x)
    {
        return x*x*x;
    }

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

    template<typename T>
    inline void inverse3x3(const T* A, T* inv)
    {

    }

    template<typename T>
    inline T constexpr dot3(const T* u, const T* v)
    {
        return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
    }

    template<typename T>
    inline void inverse3x3(const T* A, T* inv)
    {
        inv[0] = A[4]*A[8] - A[5]*A[7];
        inv[3] = -(A[3] * A[8] - A[5] * A[6]);
        inv[6] = A[3] * A[7] - A[4] * A[6];
        T det = 1.0 / (A[0]*inv[0] + A[1]*inv[3] + A[2]*inv[6]);
        inv[0] *= det;
        inv[3] *= det;
        inv[6] *= det;

        inv[1] = -det * (A[1]*A[8] - A[2]*A[7]);
        inv[4] = det * (A[0]*A[8] - A[2]*A[6]);
        inv[7] = -det * (A[0]*A[7] - A[1]*A[6]);

        inv[2] = det * (A[1]*A[5] - A[2]*A[4]);
        inv[5] = -det * (A[0]*A[5] - A[2]*A[3]);
        inv[8] = det * (A[0]*A[4] - A[1]*A[3]);

        return out;
    }

    template<typename T>
    inline void mult3x3(const T* A, const T* x, T* y)
    {
        y[0] = dot3(A, x);
        y[1] = dot3(A + 3, x);
        y[2] = dot3(A + 6, x);
    }
}
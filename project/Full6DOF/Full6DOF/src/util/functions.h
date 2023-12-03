#include <cmath>
#include <array>

namespace func
{

    template<typename T>
    constexpr T SQ(T x)
    {
        return x*x;
    }

    template<typename T>
    constexpr T CB(T x)
    {
        return x*x*x;
    }

    template<typename T>
    inline T constexpr dot2(const T* u, const T* v)
    {
        return u[0]*v[0] + u[1]*v[1];
    }

    template<typename T>
    inline void inverse2x2(const T* A, T* inv)
    {
        T det = 1.0/(A[0]*A[3] - A[1]*A[2]);
        inv[0] = A[3]*det;
        inv[1] = -A[1]*det;
        inv[2] = -A[2]*det;
        inv[3] = A[0]*det;
    }

    template<typename T>
    inline void mult2x2(const T* A, const T* x, T* y)
    {
        y[0] = dot2(A, x);
        y[1] = dot2(A + 2, x);
    }

    template<typename T>
    inline T constexpr dot3(const T* u, const T* v)
    {
        return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
    }

    template<typename T>
    inline T constexpr mag3(const T* u)
    {
        return sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
    }

    template<typename T>
    inline void cross(const T* u, const T* v, T* w)
    {
        w[0] = u[1]*v[2] - u[2]*v[1];
        w[1] = u[2]*v[0] - u[0]*v[2];
        w[2] = u[0]*v[1] - u[1]*v[0];
    }

    template<typename T>
    void inverse3x3(const T* A, T* inv);

    template<typename T>
    inline void mult3x3(const T* A, const T* x, T* y)
    {
        y[0] = dot3(A, x);
        y[1] = dot3(A + 3, x);
        y[2] = dot3(A + 6, x);
    }

    template<typename T>
    inline T dot(const T* u, const T* v, const std::size_t& N);

    template<typename T, uint_fast16_t N>
    inline T dot(const T* u, const T* v)
    {
        T sum = 0;
        for (auto i = 0u; i < N; i++)
        {
            sum += (*u++)*(*v++);
        }
        return sum;
    }

    template<typename T>
    void mult(const T* A, const T* B, T* C, const uint_fast16_t n);

    template<typename T, uint_fast16_t N>
    inline void mult(const T* A, const T* B, T* C)
    {
        T* A_i = A;
        T* C_ik = C;
        for (uint_fast16_t i = 0; i < N; i++)
        {
            for (uint_fast16_t j = 0; j < N; j++)
            {
                T* B_j = B + j;
                for (uint_fast16_t k = 0; k < N; k++)
                {
                    *C_ik += A_i[k] * (*B_j);
                    B_j += N;
                }
                C_ik++;
            }
            A_i += N;
        }
    }

    template<typename T>
    inline void mult(const T* A, const T b, T* C, const uint_fast16_t n)
    {
        const uint_fast32_t n_sq = n*n;
        for (auto i = 0u; i < n_sq; i++)
        {
            *C++ = b*(*A++);
        }
    }

    template<typename T, uint_fast16_t N>
    inline void mult(const T* A, const T b, T* C)
    {
        constexpr uint_fast32_t n_sq = N*N;
        for (auto i = 0u; i < n_sq; i++)
        {
            *C++ = b * (*A++);
        }
    }

    template<typename T>
    inline void add(const T* A, const T* B, T* C, const uint_fast16_t n)
    {
        const uint_fast32_t n_sq = n*n;
        for (auto i = 0u; i < n; i++) 
        {
            *C++ = (*A++) + (*B++);
        }
    }

    template<typename T, uint_fast16_t N>
    inline void add(const T* A, const T* B, T* C)
    {
        constexpr uint_fast32_t n_sq = N*N;
        for (std::size_t i = 0; i < N; i++)
        {
            *C++ = (*A++) + (*B++);
        }
    }

    template<typename T>
    void rotateZ(const T angle, const T* u, T* v)
    {
        T sA = sin(angle);
        T cA = cos(angle);
        v[0] = u[0]*cA - u[1]*sA;
        v[1] = u[1]*cA + u[0]*sA;
        v[2] = u[2];
    }

    template<typename T>
    void LUPSolve(T* A_, T* b, const uint_fast16_t n);

    template<typename T>
    void triDiagonalSolve(T a[], T b[], T c[], T x[], const uint_fast16_t nDiag);

    template<typename T>
    inline void LTS(const T* L, T* b, const uint_fast16_t n) 
    {
        //lower triangular solver by forward substitution
        for (uint_fast16_t i = 0; i < n; i++) 
        {
            for (uint_fast16_t j = 0; j < i; j++)
            {
                b[i] -= L[j] * b[j];
            }
            b[i] /= L[i];
            L += n;
        }
    }

    template<typename T>
    inline void UTS(const T* U, T* b, const uint_fast16_t n)
    {
        //lower triangular solver by backward substitution
        for (uint_fast16_t i = n - 1; i != 0; i--)
        {
            for (uint_fast16_t j = i + 1; j < n; j++) 
            {
                b[i] -= U[j] * b[j];
            }
            b[i] /= U[i];
            U += n;
        }
    }

    template<typename T>
    void getLU(const T* A, T* L, T* U, const uint_fast16_t& n);

    template<typename T>
    T* cholesky(T* A, const uint_fast16_t n);

    constexpr unsigned long long factorial[] = { 1,1,2,6,24,120,720,5040,40320,362880,3628800,39916800,479001600,6227020800,87178291200,1307674368000,20922789888000,355687428096000,6402373705728000,121645100408832000,2432902008176640000 };

    constexpr int arithmetics_primes[] = { 2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,0 };

    constexpr int C[11][11] = { {1},{1,1},{1,2,1},{1,3,3,1},{1,4,6,4,1},{1,5,10,10,5,1},{1,6,15,20,15,6,1},{1,7,21,35,35,21,7,1},{1,8,28,56,70,56,28,8,1},{1,9,36,84,126,126,36,9,1},{1,10,45,120,210,252,210,120,45,10,1} };

    inline double generalBinomial(double alpha, int k) 
    {
        // this can be further optimized for half values required by legendre
        double res = 1;
        for (int i = 1; i <= k; ++i)
            res *= (alpha - static_cast<double>(k + i)) / static_cast<double>(i);
        return res;
    }

    constexpr inline unsigned long binomial(unsigned long& n, unsigned long& k) 
    {
        unsigned long c = 1, i;

        if (k > n - k) // take advantage of symmetry
            k = n - k;

        for (i = 1; i <= k; i++, n--) {
            if (c / i > 4294967295 / n) // return 0 on overflow
                return 0;

            c = c/i*n + c % i*n/i;  // split c * n / i into (c / i * i + c % i) * n / i
        }

        return c;
    }

    unsigned long combination(const unsigned short& n, const unsigned short& k);

    double legendrePoly(const int n, const double x);

    double assocLegendrePoly(int l, int m, double x);

    inline double evalPoly(const double coef[], const int& n, const double& x) {
        int idx = n - 1;
        double y = coef[idx];
        while (idx-- > 0) 
        {
            y = y*x + coef[idx];
        }
        return y;
    }

    template<typename T>
    inline void polyFit(T* x, T* y, const uint_fast16_t& n, const uint_fast16_t& deg, T* coef);

    template<typename T>
    void weightedPolyFit(const T x[], const T y[], const T w[], const uint_fast16_t& n, const uint_fast16_t& deg, T* coef);

    void XYZRotation(double angleX, double angleY, double angleZ, double* mat)
    {
        double cX = cos(angleX);
        double cY = cos(angleY);
        double cZ = cos(angleZ);
        double sX = sin(angleX);
        double sY = sin(angleY);
        double sZ = sin(angleZ);


    }

    void ZYXRotation(double angleX, double angleY, double angleZ, double* mat)
    {
        double cX = cos(angleX);
        double cY = cos(angleY);
        double cZ = cos(angleZ);
        double sX = sin(angleX);
        double sY = sin(angleY);
        double sZ = sin(angleZ);
        mat[0] = cZ*cY;
        mat[1] = cZ*sY*sX - sZ*cX;
        mat[2] = cZ*sY*cX + sZ*sX;
        mat[3] = sZ*cY;
        mat[4] = sZ*sY*sX + cZ*cX;
        mat[5] = sZ*sY*cX - cZ*sX;
        mat[6] = -sY;
        mat[7] = cY*sX;
        mat[8] = cY*cX;
    }

}
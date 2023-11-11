#include "functions.h"

namespace func
{
    template<typename T>
    void inverse3x3(const T* A, T* inv)
    {
        inv[0] = A[4] * A[8] - A[5] * A[7];
        inv[3] = -(A[3] * A[8] - A[5] * A[6]);
        inv[6] = A[3] * A[7] - A[4] * A[6];
        T det = 1.0 / (A[0] * inv[0] + A[1] * inv[3] + A[2] * inv[6]);
        inv[0] *= det;
        inv[3] *= det;
        inv[6] *= det;

        inv[1] = -det * (A[1] * A[8] - A[2] * A[7]);
        inv[4] = det * (A[0] * A[8] - A[2] * A[6]);
        inv[7] = -det * (A[0] * A[7] - A[1] * A[6]);

        inv[2] = det * (A[1] * A[5] - A[2] * A[4]);
        inv[5] = -det * (A[0] * A[5] - A[2] * A[3]);
        inv[8] = det * (A[0] * A[4] - A[1] * A[3]);
    }

    
    template<> void inverse3x3<double>(const double* A, double* inv);
    template<> void inverse3x3<float>(const float* A, float* inv);

    template<typename T>
    T dot(const T* u, const T* v, const std::size_t& N)
    {
        T sum = 0;
        for (auto i = 0u; i < N; i++)
        {
            sum += u[i] * v[i];
        }
        return sum;
    }

    template<> double dot(const double* u, const double* v, const std::size_t& N);
    template<> float dot(const float* u, const float* v, const std::size_t& N);
    template<> int dot(const int* u, const int* v, const std::size_t& N);

    template<typename T>
    void mult(const T* A, const T* B, T* C, const uint_fast16_t n)
    {
        T* A_i = A;
        T* C_ik = C;
        for (uint_fast16_t i = 0; i < n; i++)
        {
            for (uint_fast16_t j = 0; j < n; j++)
            {
                T* B_j = B + j;
                for (uint_fast16_t k = 0; k < n; k++)
                {
                    *C_ik += A_i[k] * (*B_j);
                    B_j += n;
                }
                C_ik++;
            }
            A_i += n;
        }
    }

    template<> void mult<double>(const double* A, const double* B, double* C, const uint_fast16_t n);
    template<> void mult<float>(const float* A, const float* B, float* C, const uint_fast16_t n);

    template<typename T>
    void LUPSolve(T* A_, T* b, const uint_fast16_t n)
    {
        uint_fast16_t i, j, k, i_max;
        T max, absA;

        T** A = new T*[n];
        for (int i = 0; i < n; i++)
        {
            A[i] = A_;
            A_ += n;
        }

        for (i = 0; i < n; i++)
        {
            max = 0;
            i_max = i;

            for (k = i; k < n; k++)
            {
                if ((absA = fabs(A[k][i])) > max)
                {
                    max = absA;
                    i_max = k;
                }
            }

            if (max < 1e-300)
            {
                throw "degenerate Matrix";
            }

            if (i_max != i)
            {
                //pivoting A
                double* rowPtr = std::move(A[i]);
                A[i] = std::move(A[i_max]);
                A[i_max] = std::move(rowPtr);

                //pivoting rows of b
                double tmp = std::move(b[i]);
                b[i] = std::move(b[i_max]);
                b[i_max] = std::move(tmp);
            }

            for (j = i + 1; j < n; j++)
            {
                A[j][i] /= A[i][i];
                for (k = i + 1; k < n; k++)
                {
                    A[j][k] -= A[j][i] * A[i][k];
                }
            }
        }

        i = 0;
        while (i < n)
        {
            for (j = 0; j < i; j++)
            {
                b[i] -= A[i][j] * b[j];
            }
            i++;
        }

        while (i > 0)
        {
            i--;
            for (j = i + 1; j < n; j++)
            {
                b[i] -= A[i][j] * b[j];
            }
            b[i] /= A[i][i];
        }
    }

    template<> void LUPSolve<double>(double* A_, double* b, const uint_fast16_t n);
    template<> void LUPSolve<float>(float* A_, float* b, const uint_fast16_t n);

    template<typename T>
    void triDiagonalSolve(T a[], T b[], T c[], T x[], const uint_fast16_t nDiag)
    {
        uint_fast16_t i1 = 0;
        uint_fast16_t i = 1;
        // forward replace
        while (i < nDiag)
        {
            T w = a[i] / b[i1];
            b[i] -= w * c[i1];
            x[i] -= w * x[i1];
            i1 = i++;
        }
        // backwards solve
        x[i] /= b[i];
        while (i > 0) {
            i--;
            x[i] = (x[i] - c[i] * x[i + 1]) / b[i];
        }
    }

    template<>
    void triDiagonalSolve<double>(double a[], double b[], double c[], double x[], const uint_fast16_t nDiag);
    template<>
    void triDiagonalSolve<float>(float a[], float b[], float c[], float x[], const uint_fast16_t nDiag);

    template<typename T>
    void getLU(const T* A, T* L, T* U, const uint_fast16_t& n)
    {
        const uint_fast32_t n_sq = n*n;
        const uint_fast16_t n1 = n + 1;
        // returns L and U
        memset(U, A, n_sq * sizeof(T));

        uint_fast16_t i, j, k;

        T* Lptr = L;
        for (k = 0; k < n; k++)
        {
            *Lptr = 1.0;
            Lptr += n1;
        }

        for (k = 0; k < n - 1; k++)
        {
            T uInv = 1.0 / U[k*n + k];
            for (i = k + 1; i < n; i++)
            {
                L[i*n + k] = U[i*n + k] * uInv;
            }
            for (i = k + 1; i < n; i++)
            {
                for (j = k + 1; j < n; j++)
                {
                    U[i*n + j] -= L[i*n + k] * U[k*n + j];
                }
            }
        }
    }

    template<>
    void getLU<double>(const double* A, double* L, double* U, const uint_fast16_t& n);
    template<>
    void getLU<float>(const float* A, float* L, float* U, const uint_fast16_t& n);

    template<typename T>
    T* cholesky(T* A, const uint_fast16_t n) 
    {
        T* L = (T*)calloc(n * n, sizeof(T));

        uint_fast32_t i_row = 0;
        for (auto i = 0; i < n; i++)
        {
            uint_fast32_t j_row = 0;
            for (auto j = 0; j < (i + 1); j++)
            {
                T s = 0;
                for (auto k = 0; k < j; k++)
                {
                    s += L[i_row + k] * L[j_row + k];
                }
                L[i_row + j] = (i == j) ? sqrt(A[i_row + i] - s) : (1.0 / L[j_row + j] * (A[i_row + j] - s));

                j_row += n;
            }
            i_row += n;
        }
        return L;
    }

    template<> double* cholesky<double>(double* A, const uint_fast16_t n);
    template<> float* cholesky<float>(float* A, const uint_fast16_t n);

    unsigned long combination(const unsigned short& n, const unsigned short& k)
    {
        if (n <= 10) 
        {
            return C[n][k];
        }
        if (k > n / 2) 
        {
            return combination(n, n - k);
        }
        unsigned long num = n; // might need long here
        unsigned long den = k;
        //vectorizable
        for (auto i = 1; i < k; i++) 
        {
            den *= i;
            num *= (n - i);
        }

        return num / den;
    }

    double legendrePoly(const int n, const double x) {
        if (n == 0)
            return 1;
        if (n == 1)
            return x;

        double sums = 0;

        for (int k = 0; k < n; k++) {
            if (k > 3) {
                sums += pow(x, k) * (combination(n, k) * generalBinomial((n + k - 1) * 0.5, n));
            }
            else {
                if (k == 0) {
                    sums += generalBinomial((n + k - 1) * 0.5, n);
                }
                else {
                    if (k == 1) {
                        sums += x * n * generalBinomial((n + k - 1) * 0.5, n);
                    }
                    else {
                        sums += x * n * generalBinomial((n + k - 1) * 0.5, n);
                    }
                }
            }
        }
        return static_cast<double>(1 << n) * sums;
    }

    double assocLegendrePoly(int l, int m, double x) 
    {
        double sums = 0;
        for (int k = m; k <= l; k++) {
            int prod = k;
            for (int j = m; m < k; m++)
                prod *= j;
            sums += prod*pow(x, k - m)*combination(l, k)*generalBinomial((l + k - 1)*0.5, l);
        }
        return static_cast<double>((1 - 2*(m % 2))*(1 << l))*pow((1.0 - x*x), static_cast<double>(m >> 1))*sums;
    }

    template<typename T>
    void polyFit(const T* x, const T* y, const uint_fast16_t& n, const uint_fast16_t& deg, T* coef)
    {
        const uint_fast16_t m = deg + 1;
        const uint_fast16_t mm = m + m;
        const uint_fast32_t m_sq = m * m;
        uint_fast16_t i, j;

        T* A = new T[m_sq];
        T* sums = new T[mm];

        std::fill_n(coef, m, 0.0);
        std::fill_n(sums, mm, 0.0);
        sums[0] = n;
        for (i = 0; i < n; i++)
        {
            coef[0] += y[i];
            T xn = x[i];
            for (j = 1; j < m; j++)
            {
                sums[j] += xn;
                coef[j] += xn * y[i];
                xn *= x[i];
            }

            for (; j < mm; j++)
            {
                sums[j] += xn;
                xn *= x[i];
            }
        }

        T* Aptr = A;
        for (i = 0; i < m; i++)
        {
            for (j = 0; j < m; j++)
            {
                *Aptr++ = sums[i + j];
            }
        }

        LUPSolve(A, coef, m);

        delete[] A;
        delete[] sums;
    }

    template<>
    void polyFit<double>(const double* x, const double* y, const uint_fast16_t& n, const uint_fast16_t& deg, double* coef);
    template<>
    void polyFit<float>(const float* x, const float* y, const uint_fast16_t& n, const uint_fast16_t& deg, float* coef);


    template<typename T>
    void weightedPolyFit(const T x[], const T y[], const T w[], const uint_fast16_t& n, const uint_fast16_t& deg, T* coef) {
        const uint_fast16_t m = deg + 1;
        const uint_fast16_t mm = m + m;
        uint_fast16_t i, j;

        T* A = new T[m*m];
        T* sums = new T[mm];

        std::fill_n(coef, m, 0.0);
        std::fill_n(sums, mm, 0.0);

        for (i = 0; i < n; i++)
        {
            sums[0] += w[i];
            coef[0] += w[i] * y[i];
            T xn = x[i];
            for (j = 1; j < m; j++)
            {
                T tmp = xn * w[i];
                sums[j] += tmp;
                coef[j] += tmp * y[i];
                xn *= x[i];
            }

            for (; j < mm; j++)
            {
                sums[j] += xn * w[j];
                xn *= x[i];
            }
        }

        T* Aptr = A;
        for (i = 0; i < m; i++)
        {
            for (j = 0; j < m; j++)
            {
                *Aptr++ = sums[i + j];
            }
        }

        LUPSolve(A, coef, m);

        delete[] A;
        delete[] sums;
    }

    template<>
    void weightedPolyFit<double>(const double* x, const double* y, const double* w, const uint_fast16_t& n, const uint_fast16_t& deg, double* coef);
    template<>
    void weightedPolyFit<float>(const float* x, const float* y, const float* w, const uint_fast16_t& n, const uint_fast16_t& deg, float* coef);
}



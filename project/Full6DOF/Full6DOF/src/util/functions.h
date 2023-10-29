#include <cmath>
#include <array>

namespace func
{

    template<typename T>
    constexpr T SQ(T x)
    {
        return x * x;
    }

    template<typename T>
    constexpr T CB(T x)
    {
        return x * x * x;
    }

    template<typename T>
    inline T constexpr dot2(const T* u, const T* v)
    {
        return u[0] * v[0] + u[1] * v[1];
    }

    template<typename T>
    inline void inverse2x2(const T* A, T* inv)
    {
        T det = 1.0 / (A[0] * A[3] - A[1] * A[2]);
        inv[0] = A[3] * det;
        inv[1] = -A[1] * det;
        inv[2] = -A[2] * det;
        inv[3] = A[0] * det;
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
        return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
    }

    template<typename T>
    inline void inverse3x3(const T* A, T* inv)
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

        return out;
    }

    template<typename T>
    inline void mult3x3(const T* A, const T* x, T* y)
    {
        y[0] = dot3(A, x);
        y[1] = dot3(A + 3, x);
        y[2] = dot3(A + 6, x);
    }

    template<typename T>
    inline T dot(const T* u, const T* v, const std::size_t& N) 
    {
        T sum = 0;
        for (auto i = 0u; i < N; i++) 
        {
            sum += u[i] * v[i];
        }
        return sum;
    }

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
    inline void mult(const T* A, const T* B, T* C, const uint_fast16_t n)
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
                    *C_ik += A_i[k]*(*B_j);
                    B_j += n;
                }
                C_ik++;
            }
            A_i += n;
        }
    }

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
                    B_j += n;
                }
                C_ik++;
            }
            A_i += n;
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
        const uint_fast32_t n_sq = n * n;
        for (auto i = 0u; i < n; i++) 
        {
            *C++ = (*A++) + (*B++);
        }
    }

    template<typename T, uint_fast16_t N>
    inline void add(const T* A, const T* B, T* C)
    {
        constexpr uint_fast32_t n_sq = N*N;
        for (std::size_t i = 0; i < n; i++)
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
    inline void LUPSolve(T** A_, T* b, const uint_fast16_t n)
    {
        uint_fast16_t i, j, k, i_max;
        T max, absA;

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

    template<typename T>
    inline void triDiagonalSolve(T a[], T b[], T c[], T x[], const uint_fast16_t nDiag) 
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
    inline void getLU(const T* A, T* L, T* U, const uint_fast16_t& n) 
    {
        const auto n_sq = n * n;
        // returns L and U
        memset(U, A, n_sq*sizeof(T))

        uint_fast16_t i, j, k;
        for (k = 0; k < n; k++)
        {
            L[k][k] = 1;
        }
        for (k = 0; k < n - 1; k++) 
        {
            for (i = k + 1; i < n; i++) 
            {
                L[i][k] = U[i][k] / U[k][k];
            }
            for (i = k + 1; i < n; i++) 
            {
                for (j = k + 1; j < n; j++) 
                {
                    U[i][j] -= L[i][k] * U[k][j];
                }
            }
        }
    }

    template<typename T>
    inline T* cholesky(T* A, const uint_fast16_t n) {
        T* L = (T*)calloc(n * n, sizeof(T));

        uint_fast32_t i_row = 0;
        for (auto i = 0; i < n; i++) 
        {
            uint_fast32_t j_row = 0;
            for (auto j = 0; j < (i + 1); j++) 
            {
                T s = 0;
                for (auto k = 0; k < j; k++) {
                    s += L[i_row + k] * L[j_row + k];
                }
                L[i_row + j] = (i == j) ? sqrt(A[i_row + i] - s) : (1.0 / L[j_row + j] * (A[i_row + j] - s));

                j_row += n;
            }
            i_row += n;
        }
        return L;
    }

    constexpr unsigned long long factorial[] = { 1,1,2,6,24,120,720,5040,40320,362880,3628800,39916800,479001600,6227020800,87178291200,1307674368000,20922789888000,355687428096000,6402373705728000,121645100408832000,2432902008176640000 };

    constexpr int arithmetics_primes[] = { 2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,0 };

    constexpr int C[11][11] = { {1},{1,1},{1,2,1},{1,3,3,1},{1,4,6,4,1},{1,5,10,10,5,1},{1,6,15,20,15,6,1},{1,7,21,35,35,21,7,1},{1,8,28,56,70,56,28,8,1},{1,9,36,84,126,126,36,9,1},{1,10,45,120,210,252,210,120,45,10,1} };

    inline double generalBinomial(double alpha, int k) {
        // this can be further optimized for half values required by legendre
        double res = 1;
        for (int i = 1; i <= k; ++i)
            res = res * (alpha - (k + i)) / i;
        return res;
    }

    inline unsigned long binomial(unsigned long& n, unsigned long& k) {
        unsigned long c = 1, i;

        if (k > n - k) // take advantage of symmetry
            k = n - k;

        for (i = 1; i <= k; i++, n--) {
            if (c / i > 4294967295 / n) // return 0 on overflow
                return 0;

            c = c / i * n + c % i * n / i;  // split c * n / i into (c / i * i + c % i) * n / i
        }

        return c;
    }

    inline int combination(const int& n, const int& k) {
        if (n <= 10) {
            return C[n][k];
        }
        if (k > n / 2) {
            return combination(n, n - k);
        }
        int num = n;
        int den = k;
        //vectorizable
        for (int i = 1; i < k; i++) {
            den *= i;
            num *= (n - i);
        }

        return num / den;
    }

    inline double legendrePoly(const int n, const double x) {
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
        return (1 << n) * sums;
    }

    inline double assocLegendrePoly(int l, int m, double x) {
        int sums = 0;
        for (int k = m; k <= l; k++) {
            int prod = k;
            for (int j = m; m < k; m++)
                prod *= j;
            sums += prod * pow(x, k - m) * combination(l, k) * generalBinomial((l + k - 1) * 0.5, l);
        }
        if (m % 2 == 0)
            return (1 << l) * pow((1 - x * x), m / 2) * sums;
        else
            return -1 * (1 << l) * pow((1 - x * x), m * 0.5) * sums;
    }

    inline double evalPoly(const double coef[], const int& n, const double& x) {
        int idx = n - 1;
        double y = coef[idx];
        while (idx-- > 0) {
            y = y * x + coef[idx];
        }
        return y;
    }

    inline void polyFit(double* x, double* y, const uint_fast16_t& n, const uint_fast16_t& deg, double* coef) {
        const uint_fast16_t m = deg + 1;
        const uint_fast16_t mm = m + m;
        uint_fast16_t i, j;

        double** A = zeros(m);
        double* sums = new double[mm];

        std::fill_n(coef, m, 0.0);
        std::fill_n(sums, mm, 0.0);
        sums[0] = n;
        for (i = 0; i < n; i++) {
            coef[0] += y[i];
            double xn = x[i];
            for (j = 1; j < m; j++) {
                sums[j] += xn;
                coef[j] += xn * y[i];
                xn *= x[i];
            }

            for (; j < mm; j++) {
                sums[j] += xn;
                xn *= x[i];
            }
        }

        for (i = 0; i < m; i++) {
            for (j = 0; j < m; j++) {
                A[i][j] = sums[i + j];
            }
        }

        LUPSolve(A, coef, m);

        del(A, m);
        delete[] sums;
    }

    inline void weightedPolyFit(const double x[], const double y[], const double w[], const uint_fast16_t& n, const uint_fast16_t& deg, double* coef) {
        const uint_fast16_t m = deg + 1;
        const uint_fast16_t mm = m + m;
        uint_fast16_t i, j;

        double** A = zeros(m);
        double* sums = new double[mm];

        std::fill_n(coef, m, 0.0);
        std::fill_n(sums, mm, 0.0);

        for (i = 0; i < n; i++) {
            sums[0] += w[i];
            coef[0] += w[i] * y[i];
            double xn = x[i];
            for (j = 1; j < m; j++) {
                double tmp = xn * w[i];
                sums[j] += tmp;
                coef[j] += tmp * y[i];
                xn *= x[i];
            }

            for (; j < mm; j++) {
                sums[j] += xn * w[j];
                xn *= x[i];
            }
        }

        for (i = 0; i < m; i++) {
            for (j = 0; j < m; j++) {
                A[i][j] = sums[i + j];
            }
        }

        LUPSolve(A, coef, m);

        del(A, m);
        delete[] sums;
    }

}
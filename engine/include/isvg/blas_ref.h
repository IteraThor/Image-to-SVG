// The handful of reference-BLAS (netlib, column-major) routines the engine calls through the BLAS library,
// with the reference loop order so floating-point results match bit-for-bit. Only the argument
// combinations the engine uses are implemented (incx = incy = 1).
#pragma once

namespace isvg::blas {

// y := alpha*op(A)*x + beta*y, A is m×n column-major with leading dimension lda.
inline void dgemv(char trans, int m, int n, double alpha, const double* A, int lda, const double* x,
                  double beta, double* y) {
    if (m == 0 || n == 0 || (alpha == 0.0 && beta == 1.0)) return;
    const int leny = trans == 'N' ? m : n;
    if (beta != 1.0) {
        if (beta == 0.0)
            for (int i = 0; i < leny; ++i) y[i] = 0.0;
        else
            for (int i = 0; i < leny; ++i) y[i] = beta * y[i];
    }
    if (alpha == 0.0) return;
    if (trans == 'N') {
        for (int j = 0; j < n; ++j) {
            const double temp = alpha * x[j];
            for (int i = 0; i < m; ++i) y[i] = y[i] + temp * A[i + j * lda];
        }
    } else {
        for (int j = 0; j < n; ++j) {
            double temp = 0.0;
            for (int i = 0; i < m; ++i) temp = temp + A[i + j * lda] * x[i];
            y[j] = y[j] + alpha * temp;
        }
    }
}

// C := alpha*A^T*B + beta*C (transa = 'T', transb = 'N'), C is m×n, inner dimension k.
inline void dgemmTN(int m, int n, int k, double alpha, const double* A, int lda, const double* B, int ldb,
                    double beta, double* C, int ldc) {
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < m; ++i) {
            double temp = 0.0;
            for (int l = 0; l < k; ++l) temp = temp + A[l + i * lda] * B[l + j * ldb];
            if (beta == 0.0)
                C[i + j * ldc] = alpha * temp;
            else
                C[i + j * ldc] = alpha * temp + beta * C[i + j * ldc];
        }
    }
}

// Solve A·X = B for a small n×n system (the BLAS library's LAPACK dgesv: dgetrf → unblocked dgetf2 since
// n < NB = 64, then dgetrs). A is overwritten by its LU factors, B (n×nrhs, ldb) by X.
// Returns LAPACK's INFO (> 0: singular, B untouched).
inline int dgesv(int n, int nrhs, double* A, int lda, int* ipiv, double* B, int ldb) {
    auto a = [&](int i, int j) -> double& { return A[i + j * lda]; };
    int info = 0;
    // dgetf2
    const double sfmin = 2.2250738585072014e-308;  // dlamch('S')
    for (int j = 0; j < n; ++j) {
        int jp = j;  // idamax over a(j..n-1, j)
        double dmax = a(j, j) < 0 ? -a(j, j) : a(j, j);
        for (int i = j + 1; i < n; ++i) {
            const double v = a(i, j) < 0 ? -a(i, j) : a(i, j);
            if (v > dmax) {
                jp = i;
                dmax = v;
            }
        }
        ipiv[j] = jp;
        if (a(jp, j) != 0.0) {
            if (jp != j)
                for (int c = 0; c < n; ++c) {
                    const double t = a(j, c);
                    a(j, c) = a(jp, c);
                    a(jp, c) = t;
                }
            if (j < n - 1) {
                const double ajj = a(j, j);
                if ((ajj < 0 ? -ajj : ajj) >= sfmin) {
                    const double r = 1.0 / ajj;
                    for (int i = j + 1; i < n; ++i) a(i, j) = r * a(i, j);
                } else {
                    for (int i = j + 1; i < n; ++i) a(i, j) = a(i, j) / ajj;
                }
            }
        } else if (info == 0) {
            info = j + 1;
        }
        if (j < n - 1) {  // dger(-1): trailing rank-1 update, column by column
            for (int c = j + 1; c < n; ++c) {
                if (a(j, c) != 0.0) {
                    const double temp = -1.0 * a(j, c);
                    for (int i = j + 1; i < n; ++i) a(i, c) = a(i, c) + a(i, j) * temp;
                }
            }
        }
    }
    if (info != 0) return info;
    // dgetrs: dlaswp, then dtrsm (lower, unit) and dtrsm (upper, non-unit)
    for (int c = 0; c < nrhs; ++c) {
        double* b = B + c * ldb;
        for (int k = 0; k < n; ++k)
            if (ipiv[k] != k) {
                const double t = b[k];
                b[k] = b[ipiv[k]];
                b[ipiv[k]] = t;
            }
        for (int k = 0; k < n; ++k)
            if (b[k] != 0.0)
                for (int i = k + 1; i < n; ++i) b[i] = b[i] - b[k] * a(i, k);
        for (int k = n - 1; k >= 0; --k)
            if (b[k] != 0.0) {
                b[k] = b[k] / a(k, k);
                for (int i = 0; i < k; ++i) b[i] = b[i] - b[k] * a(i, k);
            }
    }
    return 0;
}

}  // namespace isvg::blas

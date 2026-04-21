/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 * Copyright (C) 2026, Si-Qiang Luo <luosq15@lzu.edu.cn>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/math/eigen.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

void eigen_tridiagonal(double **a, int n, double *d, double *e, double *et, int lt)
{
    /*
     * ALGORITHM: Implicit QR Method with Householder Reduction
     * 
     * This function finds eigenvalues and eigenvectors of a symmetric matrix using:
     * 1. Householder transformation to reduce matrix to tridiagonal form
     * 2. Implicit QR iteration to find eigenvalues
     * 3. Back-substitution to find selected eigenvectors
     * 
     * Input:  a[n][n] = symmetric matrix (input matrix is modified during computation)
     * Output: d[n]    = eigenvalues in ascending order
     *         e[n]    = workspace (subdiagonal elements during computation)
     *         et[lt]  = subdiagonal elements for selected eigenvectors
     *         vt      = selected eigenvectors (stored in rows of a after computation)
     * Parameters: lt = number of eigenvectors to compute (if lt < n-1, only lt eigenvectors computed)
     */
    int i, j, k;
    double sigma, beta, uu, eps, vmax;

    double r, c, s, t, ks;
    double d0, d1, e0, e1;
    double vc, vs;
    double *dd, *ee, *ds, *es, *as;
    int *od;
    double temp, g;

    eps = 1E-17;
    e[0] = 0;
    
    /* ========== PHASE 1: HOUSEHOLDER TRIDIAGONALIZATION ==========
     * Reduce matrix a to tridiagonal form using Householder reflections.
     * Process rows from bottom to top (i = n-1 down to 1).
     * 
     * For each row i, we eliminate all elements below the subdiagonal:
     * - Compute Householder vector v from a[i][0..i-1]
     * - Apply reflection: a = (I - 2*v*v^T)*a*(I - 2*v*v^T)^T
     * - Store the Householder coefficients for later eigenvector computation
     */
    for (i = n - 1; i >= 1; i--) {
        /* Compute norm of row i below the diagonal: sqrt(a[i][0]^2 + ... + a[i][i-1]^2) */
        sigma = 0;
        for (j = 0; j <= i - 1; j++) {
            sigma += a[i][j] * a[i][j];
        }
        sigma = sqrt(sigma);
        
        /* Choose sign to avoid cancellation */
        if (a[i][i - 1] < 0) {
            sigma *= -1;
        }
        
        /* Store diagonal and subdiagonal elements */
        d[i] = a[i][i];
        e[i] = -sigma;
        
        /* Householder vector scalar (beta = sigma * (sigma + a[i][i-1])) */
        beta = sigma * (sigma + a[i][i - 1]);
        a[i][i - 1] += sigma;

        /* Apply Householder transformation if beta is significant */
        if (fabs(beta) > eps * fabs(a[i][i])) {
            /* Compute a[0..i-1, i] = A * v, where v = [a[i][0], ..., a[i][i-1]]^T */
            for (j = 0; j < i; j++) {
                a[j][i] = 0;
                for (k = 0; k < i; k++)
                {
                    a[j][i] += a[j][k] * a[i][k];
                }
            }

            /* Compute u = (a * v) / (2 * beta) */
            uu = 0;
            for (j = 0; j < i; j++) {
                uu += a[i][j] * a[j][i];
            }
            uu = uu / beta / beta;

            /* Apply transformation: a = a - v*(a*v)^T/beta - (a*v)*v^T/beta + (a*v)*(u*v^T + v*u^T) */
            for (j = 0; j < i; j++) {
                for (k = 0; k < i; k++) {
                    a[j][k] = a[j][k] - a[i][j] * a[k][i] / beta - a[i][k] * a[j][i] / beta + a[i][j] * a[i][k] * uu;
                }
            }
            a[i][i] = beta;
        }
        else {
            a[i][i] = 0;
        }
    }
    d[0] = a[0][0];

    /* ========== PHASE 2: EIGENVECTOR COMPUTATION (if lt > 0) ==========
     * Back-calculate eigenvectors from Householder transformations.
     * 
     * If lt > 0, transform the stored Householder vectors into eigenvectors
     * of the original matrix by applying accumulated Householder reflections.
     * Initialize identity matrix and apply reflections in reverse order.
     */
    if (lt > 0) {
        a[0][0] = 1;  /* Initialize first row to identity */
        for (i = 1; i < n; i++) {
            if (0 != a[i][i]) {
                /* Apply Householder reflection for row i */
                for (j = 0; j < i; j++) {
                    a[j][i] = 0;
                    for (k = 0; k < i; k++) {
                        a[j][i] += a[j][k] * a[i][k];
                    }
                    a[j][i] /= a[i][i];
                }
                for (j = 0; j < i; j++) {
                    for (k = 0; k < i; k++) {
                        a[j][k] -= a[j][i] * a[i][k];
                    }
                }
                for (j = 0; j < i; j++) {
                    a[i][j] = 0;
                    a[j][i] = 0;
                }
            }
            a[i][i] = 1;  /* Identity matrix row */
        }
    }

    /* Allocate working arrays for QR iteration */
    dd = (double *)malloc(sizeof(double) * n);
    ee = (double *)malloc(sizeof(double) * n);
    ds = (double *)malloc(sizeof(double) * n);
    es = (double *)malloc(sizeof(double) * n);
    as = (double *)malloc(sizeof(double) * n);
    od = (int *)malloc(sizeof(int) * n);
    for (i = 0; i < n; i++) {
        ds[i] = d[i];  /* Save original d */
        es[i] = e[i];  /* Save original e */
        od[i] = i;     /* Initialize ordering array */
    }

    t = 0;
    e[0] = 0;
    
    /* ========== PHASE 3: IMPLICIT QR ITERATION ==========
     * Find eigenvalues using implicit QR algorithm.
     * 
     * For each tridiagonal block, perform QR iterations:
     * - Compute shift ks using Rayleigh quotient acceleration
     * - Perform Givens rotations to zero subdiagonal elements
     * - Update eigenvalues d[] and subdiagonal e[]
     * - Accumulate transformations for eigenvector computation
     */
    for (i = 0; i < n - 1; i++) {
        while (fabs(e[i + 1]) > eps * (fabs(d[i]) + fabs(d[i + 1]))) {
            /* Compute shift using Rayleigh quotient acceleration */
            g = (d[i + 1] - d[i]) / (2 * e[i + 1]);
            if (g >= 0) {
                ks = d[i] - e[i + 1] / (g + sqrt(1 + g * g));
            }
            else {
                ks = d[i] - e[i + 1] / (g - sqrt(1 + g * g));
            }
            
            /* Implicit QR step: apply Givens rotations to eliminate subdiagonal */
            for (j = n - 1; j > i; j--) {
                /* Compute Givens rotation (c, s) for 2x2 block */
                r = sqrt((d[j] - ks) * (d[j] - ks) + e[j] * e[j]);
                if (0 == r) {
                    continue;
                }
                c = (d[j] - ks) / r;      /* cos(theta) */
                s = -e[j] / r;             /* sin(theta) */

                /* Save old diagonal and subdiagonal */
                d0 = d[j - 1];
                d1 = d[j];
                e0 = e[j - 1];
                e1 = e[j];

                /* Apply Givens rotation: T' = G^T * T * G */
                d[j - 1] = c * c * d0 + s * s * d1 + 2 * c * s * e1;
                d[j] = s * s * d0 + c * c * d1 - 2 * c * s * e1;
                e[j - 1] = c * e0;
                e[j] = c * s * (-d0 + d1) + (c * c - s * s) * e1;

                if (j < n - 1) {
                    e[j + 1] = -s * t + c * e[j + 1];
                }
                t = -s * e0;

                /* If computing all eigenvectors, accumulate transformation */
                if (lt >= n - 1) {
                    for (k = 0; k < n; k++) {
                        vc = a[j - 1][k];
                        vs = a[j][k];
                        a[j - 1][k] = c * vc + s * vs;
                        a[j][k] = -s * vc + c * vs;
                    }
                }
            }
        }
    }

    /* ========== PHASE 4: SORT EIGENVALUES ==========
     * Sort eigenvalues in ascending order and track corresponding eigenvectors.
     * Use bubble sort on eigenvalues.
     */
    for (i = 0; i < n - 1; i++) {
        for (j = 0; j < n - i - 1; j++) {
            if (d[j] > d[j + 1]) {
                /* Swap eigenvalues */
                temp = d[j];
                d[j] = d[j + 1];
                d[j + 1] = temp;

                /* Update ordering array */
                k = od[j];
                od[j] = od[j + 1];
                od[j + 1] = k;
            }
        }
    }

    /* Save sorted eigenvalues and subdiagonal */
    for (i = 0; i < n; i++) {
        dd[i] = d[i];
        ee[i] = e[i];
    }

    /* ========== PHASE 5: PARTIAL EIGENVECTOR COMPUTATION ==========
     * If lt < n-1, compute only the first lt eigenvectors using iterative refinement.
     * This is more efficient than computing all eigenvectors.
     */
    if (lt >= n - 1) {
        /* Reorder eigenvectors according to sorted eigenvalues */
        for (j = 0; j < n; j++) {
            for (i = 0; i < n; i++) {
                as[i] = a[i][j];
            }
            for (i = 0; i < n; i++) {
                a[i][j] = as[od[i]];
            }
        }
    }
    else {
        /* For partial eigenvectors: reset and use inverse iteration */
        for (i = 0; i < n; i++) {
            od[i] = -1;
            d[i] = ds[i];   /* Restore original unsorted eigenvalues */
            e[i] = es[i];   /* Restore original subdiagonal */
        }

        t = 0;
        e[0] = 0;
        
        /* Inverse iteration to find specific eigenvectors */
        for (i = 0; i < lt; i++) {
            /* Shift the tridiagonal matrix by -dd[i] to find eigenvector for eigenvalue dd[i] */
            for (j = 0; j < n; j++) {
                d[j] -= dd[i];
            }
            
            /* Iterate until we find the eigenvalue closest to zero */
            while (1) {
                int bk = 0;
                
                /* Check if any diagonal element is nearly zero */
                for (j = 0; j < n - 1; j++) {
                    if (fabs(e[j + 1]) <= eps * (fabs(d[j]) + fabs(d[j + 1])) && 
                        fabs(d[j]) <= eps * pow(10, 8) * n * fabs(dd[i]) && od[j] < 0) {
                        od[j] = i;
                        bk = 1;
                        break;
                    }
                }
                if (fabs(d[n - 1]) <= eps * pow(10, 8) * n * fabs(dd[i]) && od[n - 1] < 0 && 0 == bk) {
                    od[n - 1] = i;
                    bk = 1;
                }
                if (1 == bk) {
                    break;
                }

                /* Apply more QR steps to converge */
                ks = 0;
                for (j = n - 1; j > 0; j--) {
                    if (od[j - 1] < 0) {
                        r = sqrt((d[j] - ks) * (d[j] - ks) + e[j] * e[j]);
                        if (0 == r) {
                            continue;
                        }
                        c = (d[j] - ks) / r;
                        s = -e[j] / r;

                        d0 = d[j - 1];
                        d1 = d[j];
                        e0 = e[j - 1];
                        e1 = e[j];

                        d[j - 1] = c * c * d0 + s * s * d1 + 2 * c * s * e1;
                        d[j] = s * s * d0 + c * c * d1 - 2 * c * s * e1;
                        e[j - 1] = c * e0;
                        e[j] = c * s * (-d0 + d1) + (c * c - s * s) * e1;

                        if (j < n - 1) {
                            e[j + 1] = -s * t + c * e[j + 1];
                        }
                        t = -s * e0;

                        for (k = 0; k < n; k++) {
                            vc = a[j - 1][k];
                            vs = a[j][k];
                            a[j - 1][k] = c * vc + s * vs;
                            a[j][k] = -s * vc + c * vs;
                        }
                    }
                }
            }
            /* Unshift: restore eigenvalues by adding back dd[i] */
            for (j = 0; j < n; j++) {
                d[j] += dd[i];
            }
        }

        /* Collect eigenvalues found */
        for (i = 0; i < n; i++) {
            as[i] = d[i];
        }

        /* Place eigenvalues in correct sorted positions */
        for (i = 0; i < n; i++) {
            if (od[i] >= 0) {
                d[od[i]] = as[i];
            }
        }
        
        /* Reorder eigenvector matrix columns */
        for (j = 0; j < n; j++) {
            for (i = 0; i < n; i++) {
                as[i] = a[i][j];
            }
            for (i = 0; i < n; i++) {
                if (od[i] >= 0) {
                    a[od[i]][j] = as[i];
                }
            }
        }
    }

    /* ========== PHASE 6: CLEANUP ==========
     * Zero out rows beyond lt (not needed eigenvectors)
     * Copy subdiagonal elements for first lt eigenvectors
     * Restore original eigenvalues
     * Normalize eigenvectors (ensure largest element is positive)
     */
    for (i = lt; i < n; i++){
        for (j = 0; j < n; j++) {
            a[i][j] = 0;
        }
    }
    for (i = 0; i < lt; i++) {
        et[i] = e[i];
    }
    for (i = 0; i < n; i++) {
        d[i] = dd[i];
        e[i] = ee[i];
    }
    
    /* Normalize eigenvectors: flip sign if largest element is negative */
    for (i = 0; i < lt; i++) {
        vmax = 0;
        for (j = 0; j < n; j++) {
            if (fabs(a[i][j]) > fabs(vmax)) {
                vmax = a[i][j];
            }
        }
        if (vmax < 0) {
            for (j = 0; j < n; j++) {
                a[i][j] *= -1;
            }
        }
    }

    free(dd);
    free(ee);
    free(ds);
    free(es);
    free(as);
    free(od);
}

void eigen_standard(double **a, int n, double *d, double **vt, int lt)
{
    /*
     * STANDARD EIGENVALUE PROBLEM: A * x = λ * x
     * 
     * Wrapper function that solves the standard symmetric eigenvalue problem.
     * 
     * Input:  a[n][n]  = symmetric matrix
     *         n        = matrix dimension
     *         lt       = number of eigenvectors to compute
     * Output: d[n]     = eigenvalues in ascending order
     *         vt[lt][n] = first lt eigenvectors (one per row)
     * 
     * Algorithm:
     * 1. Create a working copy of matrix a (to avoid modifying the input)
     * 2. Call eigen_tridiagonal() which uses Householder reduction + implicit QR
     * 3. Extract eigenvectors from the first lt rows of the modified matrix
     * 4. Free temporary arrays
     */
    if (vt == NULL) {
        lt = 0;
    }
    if (lt <= 0) {
        vt = NULL;
    }

    double **aa, *e, *et;
    int i, j;
    
    /* Allocate working copy of matrix and temporary arrays */
    aa = (double **)malloc(sizeof(double *) * n);
    e = (double *)malloc(sizeof(double) * n);
    et = (double *)malloc(sizeof(double) * lt);
    
    for (i = 0; i < n; i++) {
        aa[i] = (double *)malloc(sizeof(double) * n);
        for (j = 0; j < n; j++) {
            aa[i][j] = a[i][j];
        }
    }

    /* Call the core eigenvalue solver */
    eigen_tridiagonal(aa, n, d, e, et, lt);

    /* Extract computed eigenvectors from working matrix */
    if (lt > 0 && vt != NULL) {
        for (i = 0; i < lt; i++) {
            for (j = 0; j < n; j++) {
                vt[i][j] = aa[i][j];
            }
        }
    }

    /* Free temporary memory */
    for (i = 0; i < n; i++) {
        free(aa[i]);
    }
    free(aa);
    free(e);
    free(et);
}

void eigen_general(double **a, double **b, int n, double *d, double **vt, int lt)
{
    /*
     * GENERALIZED EIGENVALUE PROBLEM: A * x = λ * B * x
     * 
     * Solves the generalized symmetric definite eigenvalue problem using:
     * - Cholesky decomposition of B
     * - Congruence transformation to reduce to standard problem
     * - Implicit QR iteration to find eigenvalues
     * 
     * Input:  a[n][n]  = symmetric matrix A
     *         b[n][n]  = symmetric positive-definite matrix B
     *         n        = matrix dimension
     *         lt       = number of eigenvectors to compute
     * Output: d[n]     = eigenvalues in ascending order
     *         vt[lt][n] = first lt eigenvectors (one per row)
     * 
     * Mathematical Basis:
     * We convert the generalized problem A*x = λ*B*x into a standard problem by:
     * 1. Decompose B = G^T * G  (Cholesky decomposition, G is lower triangular)
     * 2. Define y = G*x, then x = G^(-1)*y
     * 3. Substitute into A*x = λ*B*x:
     *    A*G^(-1)*y = λ*G^T*G*G^(-1)*y = λ*G^T*y
     * 4. Multiply by (G^(-T)) from left:
     *    (G^(-T))*A*G^(-1)*y = λ*y  (standard problem!)
     * 5. Define S = (G^(-T))*A*G^(-1)
     * 6. Solve standard problem: S*z = λ*z
     * 7. Transform back: x = G^(-1)*z
     */
    if (vt == NULL) {
        lt = 0;
    }
    if (lt <= 0) {
        vt = NULL;
    }

    double **G, **IG, **IGA, **S, *e, *et;
    int i, j, k, ii;
    double s, ds;

    /* Allocate working matrices */
    G = (double **)malloc(sizeof(double *) * n);
    IG = (double **)malloc(sizeof(double *) * n);
    IGA = (double **)malloc(sizeof(double *) * n);
    S = (double **)malloc(sizeof(double *) * n);
    e = (double *)malloc(sizeof(double) * n);
    et = (double *)malloc(sizeof(double) * lt);

    for (i = 0; i < n; i++) {
        G[i] = (double *)malloc(sizeof(double) * n);
        IG[i] = (double *)malloc(sizeof(double) * n);
        IGA[i] = (double *)malloc(sizeof(double) * n);
        S[i] = (double *)malloc(sizeof(double) * n);
    }

    /* ========== STEP 1: CHOLESKY DECOMPOSITION OF B ==========
     * Decompose B = G^T * G where G is lower triangular matrix.
     * 
     * Algorithm: For j = 0 to n-1:
     *   G[j][j] = sqrt(B[j][j] - sum_{k=0}^{j-1} G[j][k]^2)
     *   For i > j:
     *     G[i][j] = (B[i][j] - sum_{k=0}^{j-1} G[i][k]*G[j][k]) / G[j][j]
     */
    for (j = 0; j < n; j++) {
        s = 0;
        for (k = 0; k <= j - 1; k++) {
            s = s + G[j][k] * G[j][k];
        }
        ds = b[j][j] - s;
        if (ds <= 0) {
            printf("error_cholesky\n");
            exit(1);
            return;
        }
        ds = fabs(ds);
        G[j][j] = sqrt(ds);

        for (i = j + 1; i < n; i++) {
            s = 0;
            for (k = 0; k <= j - 1; k++) {
                s = s + G[i][k] * G[j][k];
            }
            G[i][j] = (b[i][j] - s) / G[j][j];
        }
    }

    /* ========== STEP 2: COMPUTE INVERSE OF G ==========
     * Solve G * IG = I column by column using back-substitution.
     * 
     * For each column ii of the identity matrix:
     * - Set IG[ii][ii] = 1, all other IG[i][ii] = 0
     * - For i = ii+1 to n-1:
     *   IG[i][ii] -= sum_{j=ii}^{i-1} G[i][j] * IG[j][ii]
     *   IG[i][ii] /= G[i][i]
     * 
     * This solves the triangular system: G[i][ii] = delta_{ii,i}
     */
    for (ii = 0; ii < n; ii++) {
        for (i = 0; i < n; i++) {
            IG[i][ii] = 0;
        }
        IG[ii][ii] = 1;
        for (i = 0; i < n; i++) {
            for (j = i + 1; j < n; j++) {
                IG[j][ii] -= G[j][i] / G[i][i] * IG[i][ii];
            }
        }
        for (i = ii; i < n; i++) {
            IG[i][ii] /= G[i][i];
        }
    }

    /* ========== STEP 3: COMPUTE IGA = IG^T * A ==========
     * Compute the product IG^T * A.
     * Since we only use lower part due to symmetry, optimize:
     *   IGA[i][j] = sum_{k=0}^{i} IG[i][k] * A[k][j]
     */
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            IGA[i][j] = 0;
            for (k = 0; k <= i; k++) {
                IGA[i][j] += IG[i][k] * a[k][j];
            }
        }
    }

    /* ========== STEP 4: COMPUTE S = IGA * IG ==========
     * Form the standard eigenvalue problem matrix:
     *   S = (IG^T * A) * IG = IG^T * A * IG
     * 
     * Since S should be symmetric, we compute only lower triangle:
     *   S[i][j] = sum_{k=0}^{j} IGA[i][k] * IG[j][k] for i >= j
     *   S[j][i] = S[i][j]  (exploit symmetry)
     */
    for (i = 0; i < n; i++) {
        for (j = 0; j <= i; j++) {
            S[i][j] = 0;
            for (k = 0; k <= j; k++) {
                S[i][j] += IGA[i][k] * IG[j][k];
            }
            S[j][i] = S[i][j];
        }
    }

    /* ========== STEP 5: SOLVE STANDARD EIGENVALUE PROBLEM ==========
     * Call eigen_tridiagonal to solve: S * z = λ * z
     * This gives eigenvalues λ and eigenvectors z of the standard problem
     */
    eigen_tridiagonal(S, n, d, e, et, lt);

    /* ========== STEP 6: TRANSFORM EIGENVECTORS BACK ==========
     * Transform eigenvectors from standard problem back to original problem:
     *   x = G^(-1) * z = IG * z
     * 
     * Since eigenvectors are stored as rows in S after eigen_tridiagonal:
     *   vt[i][j] = sum_{k=j}^{n-1} S[i][k] * IG[k][j]
     */
    if (lt > 0 && vt != NULL) {
        for (j = 0; j < n; j++) {
            for (i = 0; i < lt; i++) {
                vt[i][j] = 0;
                for (k = j; k < n; k++) {
                    vt[i][j] += S[i][k] * IG[k][j];
                }
            }
        }
    }

    /* ========== CLEANUP: FREE MEMORY ==========  */
    for (i = 0; i < n; i++) {
        free(G[i]);
        free(IG[i]);
        free(IGA[i]);
        free(S[i]);
    }
    free(G);
    free(IG);
    free(IGA);
    free(S);
    free(e);
    free(et);
}

#ifdef LAPACKE

#include <stdio.h>
#include <lapacke.h>

void lapack_general(double **a, double **b, int n, double *e, double **vt, int lt)
{
    lapack_int N = n;          /* dimension of the matrices */
    double *A = (double *)malloc(sizeof(double) * N * N);
    double *B = (double *)malloc(sizeof(double) * N * N);
    lapack_int info;

    /* copy 2D arrays to 1D arrays in column-major order */
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
            A[j * n + i] = a[i][j];
            B[j * n + i] = b[i][j];
        }
    }

    /* calculate eigenvalues and right eigenvectors */
    info = LAPACKE_dsygv(LAPACK_COL_MAJOR,     /* storage： col major */
                         1,                    /* ITYPE=1： A x = λ B x */
                         'V',                  /* calculate eigenvector 'V' / or just eigenvalue 'N' */
                         'U',                  /* use up triangular 'U' or 'L' */
                         N,
                         A, N,
                         B, N,
                         e);

    if (info == 0) {
        /* copy the eigenvectors from the 1D array to the 2D array */
        if (vt != NULL) {
             for (int i = 0; i < lt; i++) {
                for (int j = 0; j < N; j++) {
                    vt[i][j] = A[i * N + j];
                }
            }
        }
    } else if (info > 0) {
        printf("DSYGV fault，info = %d\n", info);
    } else {
        printf("Illegal argument，info = %d\n", info);
    }
}

#endif
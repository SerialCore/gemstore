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

void eigen_tridiagonal_complex(double complex **a, int n, double *d, double *e, double *et, int lt)
{
    /* 
     * Find eigenvalues and eigenvectors using Jacobi method on Hermitian matrix.
     * For complex Hermitian matrices, we use complex Givens rotations
     * derived from the eigenvectors of the 2x2 blocks.
     */
    int i, j, p, q, iter, k;
    double complex **mat;
    double complex **U;  /* Cumulative transformation matrix: eigenvectors columns */
    double complex c, s;  /* Complex Givens coefficients */
    double eps = 1E-12;
    double complex apq, app, aqq;
    double max_elem;
    int pmax, qmax;
    
    /* Allocate working copy for matrix */
    mat = (double complex **)malloc(n * sizeof(double complex *));
    for (i = 0; i < n; i++) {
        mat[i] = (double complex *)malloc(n * sizeof(double complex));
        for (j = 0; j < n; j++) {
            mat[i][j] = a[i][j];
        }
    }
    
    /* Allocate and initialize eigenvector matrix U to identity */
    U = (double complex **)malloc(n * sizeof(double complex *));
    for (i = 0; i < n; i++) {
        U[i] = (double complex *)malloc(n * sizeof(double complex));
        for (j = 0; j < n; j++) {
            U[i][j] = (i == j) ? 1.0 + 0.0*I : 0.0 + 0.0*I;
        }
    }
    
    /* Initialize */
    e[0] = 0.0;
    for (i = 1; i < n; i++) {
        e[i] = 0.0;
    }
    
    /* Jacobi eigenvalue iteration */
    for (iter = 0; iter < 1000; iter++) {
        /* Find maximum off-diagonal element */
        max_elem = 0.0;
        pmax = 0;
        qmax = 1;
        
        for (p = 0; p < n; p++) {
            for (q = p + 1; q < n; q++) {
                double elem = cabs(mat[p][q]);
                if (elem > max_elem) {
                    max_elem = elem;
                    pmax = p;
                    qmax = q;
                }
            }
        }
        
        if (max_elem < 1E-14) break;  /* Converged to high precision */
        
        p = pmax;
        q = qmax;
        
        /* Get the 2x2 block */
        apq = mat[p][q];
        app = mat[p][p];
        aqq = mat[q][q];
        
        if (cabs(apq) < eps) {
            /* Already diagonal block, skip */
            continue;
        }
        
        /*
         * For complex Hermitian Jacobi with eigendecomposition:
         * Compute the eigenvalue of the 2x2 block using quadratic formula.
         * Then extract the eigenvector which becomes the first column of G.
         * 
         * For 2x2 Hermitian with diagonal app, aqq and off-diagonal apq:
         * trace = app + aqq
         * det = app*aqq - |apq|^2
         * disc = sqrt((trace/2)^2 - det)
         * lambda_1 = trace/2 + disc
         * 
         * Eigenvector for lambda_1 is proportional to [apq, aqq - lambda_1]
         * After normalization, this becomes [c, -s*] where G = [c s; -s* c*]
         */
        
        double complex trace = app + aqq;
        double complex det = app * aqq - apq * conj(apq);
        double complex disc_term = (trace * trace) / 4.0 - det;
        double complex disc = csqrt(disc_term);
        double complex lambda1 = trace / 2.0 + disc;
        
        /* Eigenvector for lambda_1: [apq, aqq - lambda_1] */
        double complex v1_p = apq;
        double complex v1_q = aqq - lambda1;
        
        /* Normalize */
        double norm_v1 = sqrt(creal(v1_p)*creal(v1_p) + cimag(v1_p)*cimag(v1_p)
                            + creal(v1_q)*creal(v1_q) + cimag(v1_q)*cimag(v1_q));
        
        if (norm_v1 > eps) {
            v1_p /= norm_v1;
            v1_q /= norm_v1;
        } else {
            /* Degenerate case, skip */
            continue;
        }
        
        /* Givens matrix has first column [v1_p, v1_q]
         * So: c = v1_p, -s* = v1_q
         * Therefore: s = -conj(v1_q)
         */
        c = v1_p;
        s = -conj(v1_q);
        
        /* Apply full rotation: A' = G^H A G */
        /* Build full G matrix */
        double complex G[n][n];
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                if (i == p && j == p) G[i][j] = c;
                else if (i == p && j == q) G[i][j] = s;
                else if (i == q && j == p) G[i][j] = -conj(s);
                else if (i == q && j == q) G[i][j] = conj(c);
                else if (i == j) G[i][j] = 1.0 + 0.0*I;
                else G[i][j] = 0.0 + 0.0*I;
            }
        }
        
        /* Compute G^H A */
        double complex GH_A[n][n];
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                GH_A[i][j] = 0.0 + 0.0*I;
                for (k = 0; k < n; k++) {
                    GH_A[i][j] += conj(G[k][i]) * mat[k][j];
                }
            }
        }
        
        /* Compute (G^H A) G */
        double complex new_mat[n][n];
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                new_mat[i][j] = 0.0 + 0.0*I;
                for (k = 0; k < n; k++) {
                    new_mat[i][j] += GH_A[i][k] * G[k][j];
                }
            }
        }
        
        /* Copy result back */
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                mat[i][j] = new_mat[i][j];
            }
        }
        
        /* Accumulate transformation: U_new = U_old * G */
        /* Where G is the Givens rotation matrix for this step */
        double complex U_new[n][n];
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                U_new[i][j] = 0.0 + 0.0*I;
                for (k = 0; k < n; k++) {
                    double complex G_kj;
                    if (k == p && j == p) G_kj = c;
                    else if (k == p && j == q) G_kj = s;
                    else if (k == q && j == p) G_kj = -conj(s);
                    else if (k == q && j == q) G_kj = conj(c);
                    else if (k == j) G_kj = 1.0 + 0.0*I;
                    else G_kj = 0.0 + 0.0*I;
                    
                    U_new[i][j] += U[i][k] * G_kj;
                }
            }
        }
        
        /* Copy U_new back to U */
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                U[i][j] = U_new[i][j];
            }
        }
    }
    
    /* Extract eigenvalues from diagonal */
    for (i = 0; i < n; i++) {
        d[i] = creal(mat[i][i]);
    }
    
    /* Eigenvectors are the columns of U (eigenvector i is U[:,i]) */
    /* Store them in the et output buffer if lt > 0 */
    /* Note: et is passed as (double *) to maintain C89 compatibility, 
     * but we'll reinterpret it as (double complex **) */
    if (lt > 0 && et != NULL) {
        double complex **vecs = (double complex **)et;
        for (i = 0; i < n && i < lt; i++) {
            for (j = 0; j < n; j++) {
                vecs[i][j] = U[j][i];  /* Column i of U */
            }
        }
    }
    
    /* Free working copies */
    for (i = 0; i < n; i++) {
        free(mat[i]);
        free(U[i]);
    }
    free(mat);
    free(U);
}

void eigen_standard_complex(double complex **a, int n, double *d, double complex **vt, int lt)
{
    if (vt == NULL) {
        lt = 0;
    }
    if (lt <= 0) {
        vt = NULL;
    }
    
    double complex **aa;
    double *e;
    double complex **et_vecs;  /* Temporary storage for eigenvectors */
    int i, j;
    
    aa = (double complex **)malloc(sizeof(double complex *) * n);
    e = (double *)malloc(sizeof(double) * n);
    
    /* Allocate eigenvector matrix if needed */
    if (lt > 0 && vt != NULL) {
        et_vecs = (double complex **)malloc(sizeof(double complex *) * lt);
        for (i = 0; i < lt; i++) {
            et_vecs[i] = (double complex *)malloc(sizeof(double complex) * n);
        }
    } else {
        et_vecs = NULL;
    }
    
    for (i = 0; i < n; i++) {
        aa[i] = (double complex *)malloc(sizeof(double complex) * n);
        for (j = 0; j < n; j++) {
            aa[i][j] = a[i][j];
        }
    }
    
    /* Pass et_vecs reinterpreted as (double *) for the C signature */
    eigen_tridiagonal_complex(aa, n, d, e, (double *)et_vecs, lt);
    
    /* Copy eigenvectors to output if provided */
    if (lt > 0 && vt != NULL && et_vecs != NULL) {
        for (i = 0; i < lt; i++) {
            for (j = 0; j < n; j++) {
                vt[i][j] = et_vecs[i][j];
            }
        }
    }
    
    for (i = 0; i < n; i++) {
        free(aa[i]);
    }
    free(aa);
    free(e);
    
    if (et_vecs != NULL) {
        for (i = 0; i < lt; i++) {
            free(et_vecs[i]);
        }
        free(et_vecs);
    }
}

void eigen_general_complex(double complex **a, double complex **b, int n, double *d, double complex **vt, int lt)
{
    if (vt == NULL) {
        lt = 0;
    }
    if (lt <= 0) {
        vt = NULL;
    }
    
    /* Fallback: Check if B is identity, then use standard solver */
    int is_identity = 1;
    for (int i = 0; i < n && is_identity; i++) {
        for (int j = 0; j < n && is_identity; j++) {
            double complex expected = (i == j) ? 1.0 : 0.0;
            if (cabs(b[i][j] - expected) > 1E-14) {
                is_identity = 0;
            }
        }
    }
    
    if (is_identity) {
        /* B is identity, use standard eigenvalue solver */
        eigen_standard_complex(a, n, d, vt, lt);
    } else {
        /* Need generalized solver - not implemented without LAPACKE */
        printf("ERROR: Generalized complex eigenvalue solver requires LAPACKE. Recompile with -DLAPACKE\n");
        for (int i = 0; i < n; i++) {
            d[i] = 0.0;
        }
    }
}


#ifdef LAPACKE

#include <stdio.h>
#include <lapacke.h>

void lapack_general_complex(double complex **a, double complex **b, int n, double *e, double complex **vt, int lt)
{
    lapack_int N = n;
    double complex *A = (double complex *)malloc(sizeof(double complex) * N * N);
    double complex *B = (double complex *)malloc(sizeof(double complex) * N * N);
    lapack_int info;
    
    /* copy 2D arrays to 1D arrays in column-major order */
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
            A[j * n + i] = a[i][j];
            B[j * n + i] = b[i][j];
        }
    }
    
    /* calculate eigenvalues and right eigenvectors using ZHEGV */
    info = LAPACKE_zhegv(LAPACK_COL_MAJOR,     /* storage: col major */
                         1,                    /* ITYPE=1: A x = λ B x */
                         'V',                  /* calculate eigenvector 'V' / or just eigenvalue 'N' */
                         'U',                  /* use upper triangular 'U' or 'L' */
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
        printf("ZHEGV fault, info = %d\n", info);
    } else {
        printf("Illegal argument, info = %d\n", info);
    }
    
    free(A);
    free(B);
}

#endif
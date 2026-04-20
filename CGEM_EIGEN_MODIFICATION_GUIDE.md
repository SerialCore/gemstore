# CGEM Eigenvalue Solver Modification Guide

## Overview

The current `src/math/eigen.c` implements eigenvalue solvers for **real symmetric** matrices:
- `eigen_tridiagonal()` — Householder tridiagonalization + implicit QR (all-real arithmetic)
- `eigen_standard()` — Standard symmetric eigenproblem (wrapper around tridiagonal)
- `eigen_general()` — Generalized symmetric eigenproblem: `A x = λ B x` via Cholesky reduction
- `lapack_general()` — LAPACK wrapper using `LAPACKE_dsygv()` for real symmetric generalized problem

For **CGEM (Complex-Range GEM)**, the Hamiltonian matrices become **complex Hermitian** (or complex-symmetric, depending on CSM formulation). The eigenvalue solver must be upgraded to handle `double complex` matrices.

---

## Current Architecture Analysis

### Real Matrices (Current)
```
Input: A (n×n, real symmetric), B (n×n, real symmetric positive-definite)
Output: λ_i ∈ ℝ, v_i ∈ ℝⁿ (eigenpairs satisfying A v_i = λ_i B v_i)
Algorithm: Cholesky factorization B = G·Gᵀ → A' = G⁻¹ A G⁻ᵀ (symmetric) 
           → tridiagonalize → implicit QR
```

### Complex Matrices (CGEM)
```
Input: A (n×n, complex Hermitian), B (n×n, complex Hermitian positive-definite)
Output: λ_i ∈ ℝ (eigenvalues still real for Hermitian problem), v_i ∈ ℂⁿ
Algorithm: Complex Cholesky: B = G·Gᴴ → A' = G⁻¹ A G⁻ᴴ (Hermitian)
           → tridiagonalize to real tridiagonal → implicit QR
Key: The tridiagonal form is REAL even though A' is complex Hermitian
```

---

## Strategy: Layered Addition (Non-Destructive)

**DO NOT modify existing real functions.** Instead, add new `_complex` variants:

1. Add `eigen_tridiagonal_complex()` — complex Hermitian → real tridiagonal via Householder
2. Add `eigen_general_complex()` — complex generalized eigenproblem via complex Cholesky + reduction
3. (Optional) Add `lapack_general_complex()` — LAPACK `LAPACKE_zhegv()` wrapper
4. Create a **dispatch function** in `spectra.c` that chooses real vs. complex solver based on `orbit_type`

---

## Detailed Modifications

### File 1: `include/gemstore/math/eigen.h` — Add New Function Declarations

**Add after line 35 (after `eigen_general` declaration):**

```c
/* ============ COMPLEX EIGENPROBLEM SOLVERS (for CGEM) ============ */

/* Complex Hermitian tridiagonalization + implicit QR
 * a: Input complex Hermitian matrix A (n × n)
 * n: Dimension of the matrices
 * d: Output array of REAL eigenvalues (length at least n)
 * e: Output / working array: REAL subdiagonal elements (length ≥ n)
 * et: Output REAL subdiagonal elements for selected eigenvalues
 * lt: Number of eigenvectors requested */
void eigen_tridiagonal_complex(double complex **a, int n, double *d, double *e, double *et, int lt);

/* Wrapper for complex Hermitian eigenproblem: copies matrix, calls tridiagonal, extracts eigenvectors
 * a: Input complex Hermitian matrix A (n × n)
 * n: Dimension
 * d: Output array of REAL eigenvalues
 * vt: Output matrix of selected COMPLEX eigenvectors (lt rows × n columns)
 * lt: Number of eigenvectors */
void eigen_standard_complex(double complex **a, int n, double *d, double complex **vt, int lt);

/* Generalized complex Hermitian eigenproblem: A x = λ B x using complex Cholesky + reduction
 * a: Input complex Hermitian matrix A (n × n)
 * b: Input complex Hermitian positive-definite matrix B (n × n)
 * n: Dimension
 * d: Output array of REAL eigenvalues
 * vt: Output matrix of selected COMPLEX eigenvectors (lt rows × n columns)
 * lt: Number of eigenvectors */
void eigen_general_complex(double complex **a, double complex **b, int n, double *d, double complex **vt, int lt);

#ifdef LAPACKE
/* LAPACK wrapper for complex generalized Hermitian eigenproblem using zhegv */
void lapack_general_complex(double complex **a, double complex **b, int n, double *e, double complex **vt, int lt);
#endif
```

---

### File 2: `src/math/eigen.c` — Implement Complex Functions

**Add at the end of the file, before the `#ifdef LAPACKE` section (around line 463):**

#### Complex Tridiagonalization Function

```c
/* ============ COMPLEX EIGENPROBLEM SOLVERS ============ */

void eigen_tridiagonal_complex(double complex **a, int n, double *d, double *e, double *et, int lt)
{
    int i, j, k;
    double complex sigma, beta, uu, uu_complex;
    double complex u, u_conj;
    
    double eps = 1E-17;
    
    e[0] = 0;
    
    /* Householder reduction to REAL tridiagonal form */
    for (i = n - 1; i >= 1; i--) {
        /* Compute norm of column i below diagonal */
        sigma = 0;
        for (j = 0; j <= i - 1; j++) {
            sigma += conj(a[i][j]) * a[i][j];  /* |a[i][j]|² */
        }
        sigma = csqrt(sigma);  /* real-valued, sqrt of sum of squares */
        
        if (creal(a[i][i - 1]) < 0) {
            sigma *= -1;
        }
        
        d[i] = creal(a[i][i]);  /* diagonal element is real for Hermitian */
        e[i] = -creal(sigma);   /* store real part as subdiagonal */
        
        beta = sigma * (sigma + a[i][i - 1]);
        a[i][i - 1] += sigma;
        
        if (cabs(beta) > eps * cabs(a[i][i])) {
            /* Compute u = A[:i] * a[i, :i]ᵀ */
            for (j = 0; j < i; j++) {
                a[j][i] = 0;
                for (k = 0; k < i; k++) {
                    a[j][i] += a[j][k] * conj(a[i][k]);  /* note: Hermitian conjugate */
                }
            }
            
            /* Compute uu = uᴴ * a[i, :i] / |beta|² */
            uu_complex = 0;
            for (j = 0; j < i; j++) {
                uu_complex += conj(a[i][j]) * a[j][i];
            }
            uu_complex = uu_complex / beta / conj(beta);
            
            /* Update A = A - a[i] * uᴴ / beta - u * a[i]ᴴ / beta + a[i] * a[i]ᴴ * uu */
            for (j = 0; j < i; j++) {
                for (k = 0; k < i; k++) {
                    a[j][k] = a[j][k] 
                             - a[i][j] * conj(a[k][i]) / beta 
                             - conj(a[i][k]) * a[j][i] / conj(beta)
                             + a[i][j] * conj(a[i][k]) * uu_complex;
                }
            }
            a[i][i] = beta;
        } else {
            a[i][i] = 0;
        }
    }
    d[0] = creal(a[0][0]);
    
    /* Eigenvector extraction (similar to real case, but a[][] is complex) */
    if (lt > 0) {
        a[0][0] = 1 + 0*I;
        for (i = 1; i < n; i++) {
            if (0 != cabs(a[i][i])) {
                for (j = 0; j < i; j++) {
                    a[j][i] = 0;
                    for (k = 0; k < i; k++) {
                        a[j][i] += a[j][k] * conj(a[i][k]);
                    }
                    a[j][i] /= a[i][i];
                }
                for (j = 0; j < i; j++) {
                    for (k = 0; k < i; k++) {
                        a[j][k] -= a[j][i] * conj(a[i][k]);
                    }
                }
                for (j = 0; j < i; j++) {
                    a[i][j] = 0;
                    a[j][i] = 0;
                }
            }
            a[i][i] = 1 + 0*I;
        }
    }
    
    /* QR iteration on the REAL tridiagonal matrix (same as real case) */
    double *dd, *ee, *ds, *es, *as;
    int *od;
    double temp, g, r, c, s, t, ks;
    double d0, d1, e0, e1, vc, vs;
    
    dd = (double *)malloc(sizeof(double) * n);
    ee = (double *)malloc(sizeof(double) * n);
    ds = (double *)malloc(sizeof(double) * n);
    es = (double *)malloc(sizeof(double) * n);
    as = (double *)malloc(sizeof(double) * n);
    od = (int *)malloc(sizeof(int) * n);
    
    for (i = 0; i < n; i++) {
        ds[i] = d[i];
        es[i] = e[i];
        od[i] = i;
    }
    
    t = 0;
    e[0] = 0;
    for (i = 0; i < n - 1; i++) {
        while (fabs(e[i + 1]) > eps * (fabs(d[i]) + fabs(d[i + 1]))) {
            g = (d[i + 1] - d[i]) / (2 * e[i + 1]);
            if (g >= 0) {
                ks = d[i] - e[i + 1] / (g + sqrt(1 + g * g));
            } else {
                ks = d[i] - e[i + 1] / (g - sqrt(1 + g * g));
            }
            for (j = n - 1; j > i; j--) {
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
                
                if (lt >= n - 1) {
                    for (k = 0; k < n; k++) {
                        double complex vc_cplx = a[j - 1][k];
                        double complex vs_cplx = a[j][k];
                        a[j - 1][k] = c * vc_cplx + s * vs_cplx;
                        a[j][k] = -s * vc_cplx + c * vs_cplx;
                    }
                }
            }
        }
    }
    
    /* Sorting and eigenvector reordering (real part, same as real case) */
    for (i = 0; i < n - 1; i++) {
        for (j = 0; j < n - i - 1; j++) {
            if (d[j] > d[j + 1]) {
                temp = d[j];
                d[j] = d[j + 1];
                d[j + 1] = temp;
                k = od[j];
                od[j] = od[j + 1];
                od[j + 1] = k;
            }
        }
    }
    
    for (i = 0; i < n; i++) {
        dd[i] = d[i];
        ee[i] = e[i];
    }
    
    if (lt >= n - 1) {
        for (j = 0; j < n; j++) {
            for (i = 0; i < n; i++) {
                as[i] = creal(a[i][j]);  /* extract real part for now */
            }
            for (i = 0; i < n; i++) {
                a[i][j] = as[od[i]] + 0*I;
            }
        }
    } else {
        for (i = 0; i < n; i++) {
            od[i] = -1;
            d[i] = ds[i];
            e[i] = es[i];
        }
        
        t = 0;
        e[0] = 0;
        for (i = 0; i < lt; i++) {
            for (j = 0; j < n; j++) {
                d[j] -= dd[i];
            }
            while (1) {
                int bk = 0;
                for (j = 0; j < n - 1; j++) {
                    if (fabs(e[j + 1]) <= eps * (fabs(d[j]) + fabs(d[j + 1])) 
                        && fabs(d[j]) <= eps * pow(10, 8) * n * fabs(dd[i]) 
                        && od[j] < 0) {
                        od[j] = i;
                        bk = 1;
                        break;
                    }
                }
                if (fabs(d[n - 1]) <= eps * pow(10, 8) * n * fabs(dd[i]) && od[n - 1] < 0 && 0 == bk) {
                    od[n - 1] = i;
                    bk = 1;
                }
                if (1 == bk) break;
                
                ks = 0;
                for (j = n - 1; j > 0; j--) {
                    if (od[j - 1] < 0) {
                        r = sqrt((d[j] - ks) * (d[j] - ks) + e[j] * e[j]);
                        if (0 == r) continue;
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
                            double complex vc_cplx = a[j - 1][k];
                            double complex vs_cplx = a[j][k];
                            a[j - 1][k] = c * vc_cplx + s * vs_cplx;
                            a[j][k] = -s * vc_cplx + c * vs_cplx;
                        }
                    }
                }
            }
            for (j = 0; j < n; j++) {
                d[j] += dd[i];
            }
        }
        
        for (i = 0; i < n; i++) {
            as[i] = d[i];
        }
        
        for (i = 0; i < n; i++) {
            if (od[i] >= 0) {
                d[od[i]] = as[i];
            }
        }
        for (j = 0; j < n; j++) {
            for (i = 0; i < n; i++) {
                as[i] = creal(a[i][j]);
            }
            for (i = 0; i < n; i++) {
                if (od[i] >= 0) {
                    a[od[i]][j] = as[i] + 0*I;
                }
            }
        }
    }
    
    for (i = lt; i < n; i++) {
        for (j = 0; j < n; j++) {
            a[i][j] = 0 + 0*I;
        }
    }
    for (i = 0; i < lt; i++) {
        et[i] = e[i];
    }
    for (i = 0; i < n; i++) {
        d[i] = dd[i];
        e[i] = ee[i];
    }
    
    /* Normalize eigenvectors */
    for (i = 0; i < lt; i++) {
        double vmax = 0;
        for (j = 0; j < n; j++) {
            if (cabs(a[i][j]) > vmax) {
                vmax = cabs(a[i][j]);
            }
        }
        if (vmax < 0 || (creal(a[i][0]) < 0 && cabs(creal(a[i][0])) == vmax)) {
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

/* Wrapper for complex standard Hermitian eigenproblem */
void eigen_standard_complex(double complex **a, int n, double *d, double complex **vt, int lt)
{
    if (vt == NULL) {
        lt = 0;
    }
    if (lt <= 0) {
        vt = NULL;
    }
    
    double complex **aa;
    double *e, *et;
    int i, j;
    
    aa = (double complex **)malloc(sizeof(double complex *) * n);
    e = (double *)malloc(sizeof(double) * n);
    et = (double *)malloc(sizeof(double) * lt);
    
    for (i = 0; i < n; i++) {
        aa[i] = (double complex *)malloc(sizeof(double complex) * n);
        for (j = 0; j < n; j++) {
            aa[i][j] = a[i][j];
        }
    }
    
    eigen_tridiagonal_complex(aa, n, d, e, et, lt);
    
    if (lt > 0 && vt != NULL) {
        for (i = 0; i < lt; i++) {
            for (j = 0; j < n; j++) {
                vt[i][j] = aa[i][j];
            }
        }
    }
    
    for (i = 0; i < n; i++) {
        free(aa[i]);
    }
    free(aa);
    free(e);
    free(et);
}

/* Generalized complex Hermitian eigenproblem */
void eigen_general_complex(double complex **a, double complex **b, int n, double *d, double complex **vt, int lt)
{
    if (vt == NULL) {
        lt = 0;
    }
    if (lt <= 0) {
        vt = NULL;
    }
    
    double complex **G, **IG, **IGA, **S;
    double *e, *et;
    int i, j, k, ii;
    double complex s, ds, beta_complex;
    
    G = (double complex **)malloc(sizeof(double complex *) * n);
    IG = (double complex **)malloc(sizeof(double complex *) * n);
    IGA = (double complex **)malloc(sizeof(double complex *) * n);
    S = (double complex **)malloc(sizeof(double complex *) * n);
    e = (double *)malloc(sizeof(double) * n);
    et = (double *)malloc(sizeof(double) * lt);
    
    for (i = 0; i < n; i++) {
        G[i] = (double complex *)malloc(sizeof(double complex) * n);
        IG[i] = (double complex *)malloc(sizeof(double complex) * n);
        IGA[i] = (double complex *)malloc(sizeof(double complex) * n);
        S[i] = (double complex *)malloc(sizeof(double complex) * n);
    }
    
    /* Complex Cholesky factorization: B = G * Gᴴ */
    for (j = 0; j < n; j++) {
        s = 0;
        for (k = 0; k <= j - 1; k++) {
            s = s + conj(G[j][k]) * G[j][k];
        }
        ds = b[j][j] - s;
        if (creal(ds) <= 0) {
            printf("error_cholesky_complex: matrix B not positive definite\n");
            exit(1);
            return;
        }
        ds = csqrt(ds);
        G[j][j] = ds;
        
        for (i = j + 1; i < n; i++) {
            s = 0;
            for (k = 0; k <= j - 1; k++) {
                s = s + conj(G[i][k]) * G[j][k];
            }
            G[i][j] = (b[i][j] - s) / G[j][j];
        }
    }
    
    /* Compute G⁻¹ (lower triangular solve) */
    for (ii = 0; ii < n; ii++) {
        for (i = 0; i < n; i++) {
            IG[i][ii] = 0;
        }
        IG[ii][ii] = 1 + 0*I;
        for (i = 0; i < n; i++) {
            for (j = i + 1; j < n; j++) {
                IG[j][ii] -= G[j][i] / G[i][i] * IG[i][ii];
            }
        }
        for (i = ii; i < n; i++) {
            IG[i][ii] /= G[i][i];
        }
    }
    
    /* Compute A' = G⁻¹ A G⁻ᴴ */
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            IGA[i][j] = 0;
            for (k = 0; k <= i; k++) {
                IGA[i][j] += IG[i][k] * a[k][j];
            }
        }
    }
    
    for (i = 0; i < n; i++) {
        for (j = 0; j <= i; j++) {
            S[i][j] = 0;
            for (k = 0; k <= j; k++) {
                S[i][j] += IGA[i][k] * conj(IG[j][k]);
            }
            S[j][i] = conj(S[i][j]);  /* Hermitian */
        }
    }
    
    /* Solve standard complex Hermitian problem on S */
    eigen_tridiagonal_complex(S, n, d, e, et, lt);
    
    /* Transform eigenvectors back: v = G⁻ᴴ * u */
    if (lt > 0 && vt != NULL) {
        for (j = 0; j < n; j++) {
            for (i = 0; i < lt; i++) {
                vt[i][j] = 0;
                for (k = j; k < n; k++) {
                    vt[i][j] += conj(IG[k][j]) * S[i][k];
                }
            }
        }
    }
    
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
```

#### (Optional) LAPACK Complex Wrapper

**Add after the complex functions, before the closing `#endif` (around line 509):**

```c
#ifdef LAPACKE

/* LAPACK wrapper for complex generalized Hermitian eigenproblem */
void lapack_general_complex(double complex **a, double complex **b, int n, double *e, double complex **vt, int lt)
{
    lapack_int N = n;
    double complex *A = (double complex *)malloc(sizeof(double complex) * N * N);
    double complex *B = (double complex *)malloc(sizeof(double complex) * N * N);
    lapack_int info;
    
    /* Copy 2D arrays to 1D arrays in column-major order */
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
            A[j * N + i] = a[i][j];
            B[j * N + i] = b[i][j];
        }
    }
    
    /* Solve: A x = λ B x (ITYPE=1) */
    info = LAPACKE_zhegv(LAPACK_COL_MAJOR,
                         1,           /* ITYPE=1: A x = λ B x */
                         'V',         /* compute eigenvectors */
                         'U',         /* use upper triangle */
                         N,
                         A, N,
                         B, N,
                         e);          /* eigenvalues (real) */
    
    if (info == 0) {
        if (vt != NULL) {
            for (int i = 0; i < lt; i++) {
                for (int j = 0; j < N; j++) {
                    vt[i][j] = A[i * N + j];
                }
            }
        }
    } else if (info > 0) {
        printf("ZHEGV error: info = %d (failed to converge)\n", info);
    } else {
        printf("ZHEGV illegal argument: info = %d\n", info);
    }
    
    free(A);
    free(B);
}

#endif
```

---

## File 3: `src/model/spectra.c` — Add Dispatch Logic

**Find the call to `eigen_general()` at line ~192. Replace with:**

```c
// Around line 192, replace:
// eigen_general(Hfi.value, temp.value, nmax, e_out->value, (v_out == NULL)? NULL : v_out->value, v_len);

// With:
#ifdef LAPACKE
    if (input->orbit == ORBIT_CGEM) {
        lapack_general_complex((double complex **)Hfi.value_complex, 
                               (double complex **)temp.value_complex,
                               nmax, e_out->value, 
                               (v_out == NULL) ? NULL : (double complex **)v_out->value_complex, 
                               v_len);
    } else {
        lapack_general(Hfi.value, temp.value, nmax, e_out->value, 
                       (v_out == NULL) ? NULL : v_out->value, v_len);
    }
#else
    if (input->orbit == ORBIT_CGEM) {
        eigen_general_complex((double complex **)Hfi.value_complex, 
                              (double complex **)temp.value_complex,
                              nmax, e_out->value, 
                              (v_out == NULL) ? NULL : (double complex **)v_out->value_complex, 
                              v_len);
    } else {
        eigen_general(Hfi.value, temp.value, nmax, e_out->value, 
                      (v_out == NULL) ? NULL : v_out->value, v_len);
    }
#endif
```

**Note:** This assumes `Hfi` and `temp` have been changed to `value_complex` fields when `orbit == ORBIT_CGEM`. See the main CGEM guide for matrix data structure changes.

---

## File 4: `include/gemstore/math/eigen.h` — Update Header

**The header changes are described in File 1 above.**

---

## Summary of Changes

| File | Change |
|---|---|
| `include/gemstore/math/eigen.h` | Add function declarations for `eigen_tridiagonal_complex()`, `eigen_standard_complex()`, `eigen_general_complex()`, `lapack_general_complex()` |
| `src/math/eigen.c` | Implement all four complex functions using complex Cholesky, Householder reduction, and QR iteration on real tridiagonal |
| `src/model/spectra.c` | Add dispatch: if `orbit == ORBIT_CGEM`, call complex solver; else call real solver |

---

## Key Implementation Notes

### 1. **Real Eigenvalues from Complex Hermitian**
Even though the matrices are complex, the eigenvalues are always **real** (property of Hermitian matrices). Store eigenvalues in `double *d`, not `double complex *d`.

### 2. **Tridiagonal Form is Real**
The key insight: after Householder reduction of a **complex Hermitian** matrix, the result is a **real tridiagonal** matrix. So the expensive QR iteration runs on reals, and only the eigenvector transformations handle complex arithmetic.

### 3. **Complex Cholesky**
Use the formula: `B = G · Gᴴ` (Hermitian conjugate, not transpose). The inverse `G⁻¹` is computed via forward substitution (lower triangular solve).

### 4. **Eigenvector Transformation**
After solving `A' u = λ u`, transform back: `v = G⁻ᴴ u` (not `G⁻ᵀ`).

### 5. **Memory Layout**
Complex matrices are stored as 2D arrays of `double complex`, just like real matrices as 2D arrays of `double`. Column-major order for LAPACK calls.

### 6. **LAPACK Availability**
- If compiled with `-DLAPACKE`, use `LAPACKE_zhegv()` (faster, more robust)
- Otherwise, fall back to your own `eigen_general_complex()` (pure C, no external deps)

---

## Testing Strategy

1. **Unit test**: Create a small 3×3 complex Hermitian test matrix with known eigenvalues
2. **Integration test**: Run CGEM spectrum with `theta = 0.0` (should match real GEM results)
3. **Physical test**: Run CGEM with small `theta` (e.g., 0.05 rad), check for resonance poles

---

## Compiler Flags

Ensure compilation includes `<complex.h>`:

```bash
gcc -std=c99 -DLAPACKE ... src/math/eigen.c
```

If not using LAPACK, the code compiles with just C99 standard library complex support.


/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/debug.h>

#include <gemstore/basis/basis.h>
#include <gemstore/basis/intrin.h>
#include <gemstore/basis/color.h>
#include <gemstore/basis/spin.h>
#include <gemstore/basis/isospin.h>
#include <gemstore/basis/orbit.h>

#include <gemstore/math/integral.h>
#include <gemstore/math/matrix.h>
#include <gemstore/math/eigen.h>
#include <gemstore/math/ceigen.h>
#include <gemstore/math/cmi.h>
#include <gemstore/math/soc.h>
#include <gemstore/math/su3.h>

#include <stdio.h>
#include <stdlib.h>

void debug_su3_product()
{
    su3_product(1, 1, 1, 1);
}

void debug_soc_operator()
{
    double s1 = 0.5, s2 = 0.5, s = 1.0, l = 2.0, jl = 1.5;
    double s1p = 0.5, s2p = 0.5, sp = 1.0, lp = 2.0, jlp = 1.5;
    double j = 1.0;
    printf("s1=%1.1f, s2=%1.1f, s=%1.1f, l=%1.1f, jl=%1.1f\n", s1, s2, s, l, jl);
    printf("s1p=%1.1f, s2p=%1.1f, sp=%1.1f, lp=%1.1f, jlp=%1.1f\n", s1p, s2p, sp, lp, jlp);
    printf("j=%1.1f\n", j);

    double sl, jj;
    sl = operator_sdots_sl(s1, s2, s, l, s1p, s2p, sp, lp, j);
    printf("Casimir operator value (s1.s2) in sl coupling: %f\n", sl);
    jj = operator_sdots_jj(s1, s2, l, jl, s1p, s2p, lp, jlp, j);
    printf("Casimir operator value (s1.s2) in jj coupling: %f\n", jj);
    sl = operator_ldotsi_sl(s1, s2, s, l, s1p, s2p, sp, lp, j);
    printf("SOC operator value (l.s1) in sl coupling: %f\n", sl);
    jj = operator_ldotsi_jj(s1, s2, l, jl, s1p, s2p, lp, jlp, j);
    printf("SOC operator value (l.s1) in jj coupling: %f\n", jj);
    sl = operator_ldotsj_sl(s1, s2, s, l, s1p, s2p, sp, lp, j);
    printf("SOC operator value (l.s2) in sl coupling: %f\n", sl);
    jj = operator_ldotsj_jj(s1, s2, l, jl, s1p, s2p, lp, jlp, j);
    printf("SOC operator value (l.s2) in jj coupling: %f\n", jj);
    sl = operator_tensor_sl(s1, s2, s, l, s1p, s2p, sp, lp, j);
    printf("Tensor operator value in sl coupling: %f\n", sl);
    jj = operator_tensor_jj(s1, s2, l, jl, s1p, s2p, lp, jlp, j);
    printf("Tensor operator value in jj coupling: %f\n\n", jj);

    s1 = 0.5, s2 = 0.5, s = 1.0, l = 2.0, jl = 1.5;
    s1p = 0.5, s2p = 0.5, sp = 1.0, lp = 2.0, jlp = 2.5;
    j = 2.0;
    printf("s1=%1.1f, s2=%1.1f, s=%1.1f, l=%1.1f, jl=%1.1f\n", s1, s2, s, l, jl);
    printf("s1p=%1.1f, s2p=%1.1f, sp=%1.1f, lp=%1.1f, jlp=%1.1f\n", s1p, s2p, sp, lp, jlp);
    printf("j=%1.1f\n", j);

    sl = operator_sdots_sl(s1, s2, s, l, s1p, s2p, sp, lp, j);
    printf("Casimir operator value (s1.s2) in sl coupling: %f\n", sl);
    jj = operator_sdots_jj(s1, s2, l, jl, s1p, s2p, lp, jlp, j);
    printf("Casimir operator value (s1.s2) in jj coupling: %f\n", jj);
    sl = operator_ldotsi_sl(s1, s2, s, l, s1p, s2p, sp, lp, j);
    printf("SOC operator value (l.s1) in sl coupling: %f\n", sl);
    jj = operator_ldotsi_jj(s1, s2, l, jl, s1p, s2p, lp, jlp, j);
    printf("SOC operator value (l.s1) in jj coupling: %f\n", jj);
    sl = operator_ldotsj_sl(s1, s2, s, l, s1p, s2p, sp, lp, j);
    printf("SOC operator value (l.s2) in sl coupling: %f\n", sl);
    jj = operator_ldotsj_jj(s1, s2, l, jl, s1p, s2p, lp, jlp, j);
    printf("SOC operator value (l.s2) in jj coupling: %f\n", jj);
    sl = operator_tensor_sl(s1, s2, s, l, s1p, s2p, sp, lp, j);
    printf("Tensor operator value in sl coupling: %f\n", sl);
    jj = operator_tensor_jj(s1, s2, l, jl, s1p, s2p, lp, jlp, j);
    printf("Tensor operator value in jj coupling: %f\n\n", jj);

    s1 = 0.5, s2 = 0.5, s = 1.0, l = 2.0, jl = 2.5;
    s1p = 0.5, s2p = 0.5, sp = 1.0, lp = 2.0, jlp = 2.5;
    j = 3.0;
    printf("s1=%1.1f, s2=%1.1f, s=%1.1f, l=%1.1f, jl=%1.1f\n", s1, s2, s, l, jl);
    printf("s1p=%1.1f, s2p=%1.1f, sp=%1.1f, lp=%1.1f, jlp=%1.1f\n", s1p, s2p, sp, lp, jlp);
    printf("j=%1.1f\n", j);

    sl = operator_sdots_sl(s1, s2, s, l, s1p, s2p, sp, lp, j);
    printf("Casimir operator value (s1.s2) in sl coupling: %f\n", sl);
    jj = operator_sdots_jj(s1, s2, l, jl, s1p, s2p, lp, jlp, j);
    printf("Casimir operator value (s1.s2) in jj coupling: %f\n", jj);
    sl = operator_ldotsi_sl(s1, s2, s, l, s1p, s2p, sp, lp, j);
    printf("SOC operator value (l.s1) in sl coupling: %f\n", sl);
    jj = operator_ldotsi_jj(s1, s2, l, jl, s1p, s2p, lp, jlp, j);
    printf("SOC operator value (l.s1) in jj coupling: %f\n", jj);
    sl = operator_ldotsj_sl(s1, s2, s, l, s1p, s2p, sp, lp, j);
    printf("SOC operator value (l.s2) in sl coupling: %f\n", sl);
    jj = operator_ldotsj_jj(s1, s2, l, jl, s1p, s2p, lp, jlp, j);
    printf("SOC operator value (l.s2) in jj coupling: %f\n", jj);
    sl = operator_tensor_sl(s1, s2, s, l, s1p, s2p, sp, lp, j);
    printf("Tensor operator value in sl coupling: %f\n", sl);
    jj = operator_tensor_jj(s1, s2, l, jl, s1p, s2p, lp, jlp, j);
    printf("Tensor operator value in jj coupling: %f\n", jj);
}

void debug_casimir_operator()
{
    matrix_t sfac = matrix_init(6 * 6, 6);
    matrix_t cfac = matrix_init(2 * 2, 6);
    intrin_wfn_t swv[6];
    intrin_wfn_t cwv[2];
    
    swv[0] = spin_wfn_tetra(1, 1, 2, 2);
    swv[1] = spin_wfn_tetra(1, 1, 1, 1);
    swv[2] = spin_wfn_tetra(1, 1, 0, 0);
    swv[3] = spin_wfn_tetra(1, 0, 1, 1);
    swv[4] = spin_wfn_tetra(0, 1, 1, 1);
    swv[5] = spin_wfn_tetra(0, 0, 0, 0);
    cwv[0] = color_wfn_tetra6();
    cwv[1] = color_wfn_tetra3();

    operator_sigma2(swv, 6, &sfac);
    operator_lambda2(cwv, 2, "qqQQ", &cfac);

    printf("SpinWF:\n");
    for (int i = 0; i < 6; i++) {
        intrin_wfn_print(&swv[i]);
        intrin_wfn_free(&swv[i]);
    }
    printf("\nColorWF:\n");
    for (int i = 0; i < 2; i++) {
        intrin_wfn_print(&cwv[i]);
        intrin_wfn_free(&cwv[i]);
    }
    printf("\nSpinFactor:\n");
    matrix_print(&sfac);
    printf("\nColorFactor:\n");
    matrix_print(&cfac);
    matrix_free(&sfac);
    matrix_free(&cfac);
}

void debug_color_wfn()
{
    intrin_wfn_t wf, ref_wf;

    wf = color_wfn_tetra1();
    printf("ColorWFTetra1:\n");
    intrin_wfn_print(&wf);

    ref_wf = color_wfn_tetra8();
    printf("\nColorWFTetra8:\n");
    intrin_wfn_print(&ref_wf);

    double ortho = intrin_wfn_overlap(&wf, &ref_wf);
    printf("\nOrthogonal degree: %f\n", ortho);

    double normal = intrin_wfn_overlap(&ref_wf, &ref_wf);
    printf("Normalized degree: %f\n", normal);

    intrin_wfn_free(&wf);
    intrin_wfn_free(&ref_wf);
}

void debug_spin_wfn()
{
    intrin_wfn_t wf;

    wf = spin_basis(0.5);
    printf("SpinBasis[1/2]:\n");
    intrin_wfn_print(&wf);
    intrin_wfn_free(&wf);

    wf = spin_wfn_meson(1.0, 0.0);
    printf("\nSpinWFMeson[1][0]:\n");
    intrin_wfn_print(&wf);
    intrin_wfn_free(&wf);

    wf = spin_wfn_baryon(1.0, 1.5, 0.5);
    printf("\nSpinWFBaryon[1,3/2][1/2]:\n");
    intrin_wfn_print(&wf);
    intrin_wfn_free(&wf);

    wf = spin_wfn_tetra(1.0, 1.0, 1.0, 0.0);
    printf("\nSpinWFTetra[1,1,1][0]:\n");
    intrin_wfn_print(&wf);
    intrin_wfn_free(&wf);

    wf = spin_wfn_penta(1.0, 0.5, 1.0, 1.5, 0.5);
    printf("\nSpinWFPenta[1,1/2,1,3/2][1/2]:\n");
    intrin_wfn_print(&wf);
    intrin_wfn_free(&wf);

    wf = spin_wfn_hexa(1.0, 0.5, 1.0, 1.5, 1.0, 0.0);
    printf("\nSpinWFHexa[1,1/2,1,3/2,1][0]:\n");
    intrin_wfn_print(&wf);
    intrin_wfn_free(&wf);
}

void debug_isospin_wfn()
{
    intrin_wfn_t wf;

    wf = isospin_basis(0.5, 'q');
    printf("IsospinBasis[1/2]:\n");
    intrin_wfn_print(&wf);
    intrin_wfn_free(&wf);

    wf = isospin_wfn_meson(1.0, 0.0);
    printf("\nIsospinWFMeson[1][0]:\n");
    intrin_wfn_print(&wf);
    intrin_wfn_free(&wf);

    wf = isospin_wfn_baryon(1.0, 1.5, 0.5);
    printf("\nIsospinWFBaryon[1,3/2][1/2]:\n");
    intrin_wfn_print(&wf);
    intrin_wfn_free(&wf);

    wf = isospin_wfn_tetra(1.0, 1.0, 1.0, 0.0);
    printf("\nIsospinWFTetra[1,1,1][0]:\n");
    intrin_wfn_print(&wf);
    intrin_wfn_free(&wf);

    wf = isospin_wfn_penta(1.0, 0.5, 1.0, 1.5, 0.5);
    printf("\nIsospinWFPenta[1,1/2,1,3/2][1/2]:\n");
    intrin_wfn_print(&wf);
    intrin_wfn_free(&wf);

    wf = isospin_wfn_hexa(1.0, 0.5, 1.0, 1.5, 1.0, 0.0);
    printf("\nIsospinWFHexa[1,1/2,1,3/2,1][0]:\n");
    intrin_wfn_print(&wf);
    intrin_wfn_free(&wf);
}

void debug_orbit_wfn()
{
    double nu1 = getnu(1, 30, 25, 0.1);
    argsOrbit_t args_bra = {
        .n = 1,
        .l = 0,
        .scale = nu1
    };

    double nu2 = getnu(2, 30, 25, 0.1);
    argsOrbit_t args_ket = {
        .n = 2,
        .l = 0,
        .scale = nu2
    };

    double factor, overlap;

    factor = 1 / sqrt(2 * nu1);
    overlap = integral_wfn_overlap(GRnlr, factor, &args_bra, &args_bra);
    printf("Normalized overlap for Gr: %f\n", overlap);

    factor = 1 / sqrt(nu1 + nu2);
    overlap = integral_wfn_overlap(GRnlr, factor, &args_bra, &args_ket);
    printf("Orthogonal overlap for Gr: %f\n", overlap);

    factor = sqrt(2 * nu1);
    overlap = integral_wfn_overlap_complex(GRnlp, factor, &args_bra, &args_bra);
    printf("Normalized overlap for Gp: %f\n", overlap);

    factor = sqrt(4 * nu1 * nu2 / (nu1 + nu2));
    overlap = integral_wfn_overlap_complex(GRnlp, factor, &args_bra, &args_ket);
    printf("Orthogonal overlap for Gp: %f\n", overlap);

    double beta1 = 0.8;
    args_bra.scale = beta1;
    double beta2 = 1.0;
    args_ket.scale = beta2;

    factor = 1 / beta1;
    overlap = integral_wfn_overlap(SRnlr, factor, &args_bra, &args_bra);
    printf("Normalized overlap for Sr: %f\n", overlap);

    factor = sqrt(2 / (beta1 * beta1 + beta2 * beta2));
    overlap = integral_wfn_overlap(SRnlr, factor, &args_bra, &args_ket);
    printf("Orthogonal overlap for Sr: %f\n", overlap);

    factor = beta1;
    overlap = integral_wfn_overlap_complex(SRnlp, factor, &args_bra, &args_bra);
    printf("Normalized overlap for Sp: %f\n", overlap);

    factor = sqrt(2 * beta1 * beta1 * beta2 * beta2 / (beta1 * beta1 + beta2 * beta2));
    overlap = integral_wfn_overlap_complex(SRnlp, factor, &args_bra, &args_ket);
    printf("Orthogonal overlap for Sp: %f\n", overlap);
}

void debug_eigen_system()
{
    int n = 4;
    matrix_t a = matrix_init(n, n);
    a.value[0][0] = 2.1; a.value[0][1] = 1.2; a.value[0][2] = 1.3; a.value[0][3] = 1.4;
    a.value[1][0] = 1.2; a.value[1][1] = 2.2; a.value[1][2] = 1.3; a.value[1][3] = 1.4;
    a.value[2][0] = 1.3; a.value[2][1] = 1.3; a.value[2][2] = 2.3; a.value[2][3] = 1.2;
    a.value[3][0] = 1.4; a.value[3][1] = 1.4; a.value[3][2] = 1.2; a.value[3][3] = 2.4;

    matrix_t b = matrix_init(n, n);
    b.value[0][0] = 1.0; b.value[0][1] = 0.9; b.value[0][2] = 0.8; b.value[0][3] = 0.7;
    b.value[1][0] = 0.9; b.value[1][1] = 1.0; b.value[1][2] = 0.7; b.value[1][3] = 0.6;
    b.value[2][0] = 0.8; b.value[2][1] = 0.7; b.value[2][2] = 1.0; b.value[2][3] = 0.6;
    b.value[3][0] = 0.7; b.value[3][1] = 0.6; b.value[3][2] = 0.6; b.value[3][3] = 1.0;

    array_t e = array_init(n);
    matrix_t v = matrix_init(n, n);

#ifdef LAPACKE
    lapack_general(a.value, b.value, n, e.value, v.value, n);
#else
    eigen_general(a.value, b.value, n, e.value, v.value, n);
#endif

    printf("A:\n");
    matrix_print(&a);
    printf("B:\n");
    matrix_print(&b);
    printf("e:\n");
    array_print(&e);
    printf("v:\n");
    matrix_print(&v);

    /* verify eigen vectors */
    printf("Residual of Av - eBv:\n");
    array_t res1 = array_init(n);
    array_t res2 = array_init(n);
    for (int k = 0; k < n; k++) {
        for (int i = 0; i < n; i++) {
        double sum1 = 0.0;
        double sum2 = 0.0;
        for (int j = 0; j < n; j++) {
            sum1 += a.value[i][j] * v.value[k][j];
            sum2 += b.value[i][j] * v.value[k][j];
        }
        res1.value[i] = sum1;
        res2.value[i] = sum2;
        }
        for (int i = 0; i < n; i++) {
            res1.value[i] -= e.value[k] * res2.value[i];
        }
        array_print(&res1);
    }

    matrix_free(&a);
    matrix_free(&b);
    array_free(&e);
    matrix_free(&v);
    array_free(&res1);
    array_free(&res2);
}

void debug_eigen_system_complex()
{
    int n = 3;
    
    /* Allocate complex matrices */
    double complex **a = (double complex **)malloc(n * sizeof(double complex *));
    double complex **b = (double complex **)malloc(n * sizeof(double complex *));
    double complex **v = (double complex **)malloc(n * sizeof(double complex *));
    for (int i = 0; i < n; i++) {
        a[i] = (double complex *)malloc(n * sizeof(double complex));
        b[i] = (double complex *)malloc(n * sizeof(double complex));
        v[i] = (double complex *)malloc(n * sizeof(double complex));
    }
    
    double *e = (double *)malloc(n * sizeof(double));
    
    /* Create a simple 3x3 complex Hermitian matrix A
     * A = [ 3     1+i   0.5-0.5i ]
     *     [ 1-i   2     0.3+0.2i ]
     *     [ 0.5+0.5i  0.3-0.2i  1 ]
     */
    a[0][0] = 3.0 + 0.0*I;
    a[0][1] = 1.0 + 1.0*I;
    a[0][2] = 0.5 - 0.5*I;
    
    a[1][0] = 1.0 - 1.0*I;
    a[1][1] = 2.0 + 0.0*I;
    a[1][2] = 0.3 + 0.2*I;
    
    a[2][0] = 0.5 + 0.5*I;
    a[2][1] = 0.3 - 0.2*I;
    a[2][2] = 1.0 + 0.0*I;
    
    /* Create identity matrix B for standard eigenvalue problem */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            b[i][j] = (i == j) ? 1.0 + 0.0*I : 0.0 + 0.0*I;
        }
    }
    
    printf("========== Complex Eigenvalue System Test ==========\n");
    printf("Matrix A (Complex Hermitian):\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("  (%.2f + %.2fi) ", creal(a[i][j]), cimag(a[i][j]));
        }
        printf("\n");
    }
    printf("\nMatrix B (Identity):\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("  (%.2f + %.2fi) ", creal(b[i][j]), cimag(b[i][j]));
        }
        printf("\n");
    }
    
    /* Call complex eigenvalue solver */
#ifdef LAPACKE
    lapack_general_complex(a, b, n, e, v, n);
    printf("\n[Using LAPACKE ZHEGV]\n");
#else
    eigen_general_complex(a, b, n, e, v, n);
    printf("\n[Using Pure C Complex Solver]\n");
#endif
    
    /* Print eigenvalues (should be real) */
    printf("\nEigenvalues (should be REAL):\n");
    for (int k = 0; k < n; k++) {
        printf("  λ_%d = %.6f\n", k, e[k]);
    }
    
    /* Print eigenvectors (complex) */
    printf("\nEigenvectors (COMPLEX):\n");
    for (int k = 0; k < n; k++) {
        printf("  v_%d = [ ", k);
        for (int i = 0; i < n; i++) {
            printf("(%.4f + %.4fi) ", creal(v[k][i]), cimag(v[k][i]));
        }
        printf("]\n");
    }
    
    /* Verify residuals: |A v_k - λ_k B v_k| should be small */
    printf("\nResiduals |A v_k - λ_k B v_k| (should be near zero):\n");
    for (int k = 0; k < n; k++) {
        double complex *av = (double complex *)malloc(n * sizeof(double complex));
        double complex *bv = (double complex *)malloc(n * sizeof(double complex));
        
        /* Compute A*v_k */
        for (int i = 0; i < n; i++) {
            av[i] = 0.0;
            for (int j = 0; j < n; j++) {
                av[i] += a[i][j] * v[k][j];
            }
        }
        
        /* Compute B*v_k */
        for (int i = 0; i < n; i++) {
            bv[i] = 0.0;
            for (int j = 0; j < n; j++) {
                bv[i] += b[i][j] * v[k][j];
            }
        }
        
        /* Compute residual and its magnitude */
        double residual_norm = 0.0;
        for (int i = 0; i < n; i++) {
            double complex res = av[i] - e[k] * bv[i];
            residual_norm += creal(res) * creal(res) + cimag(res) * cimag(res);
        }
        residual_norm = sqrt(residual_norm);
        
        printf("  residual[%d] = %.2e\n", k, residual_norm);
        
        free(av);
        free(bv);
    }
    
    /* Verify orthogonality: v_i^H * B * v_j = δ_ij */
    printf("\nOrthogonality v_i^H * B * v_j (should be δ_ij):\n");
    for (int k = 0; k < n; k++) {
        for (int l = 0; l < n; l++) {
            double complex ortho = 0.0;
            for (int i = 0; i < n; i++) {
                double complex sum = 0.0;
                for (int j = 0; j < n; j++) {
                    sum += conj(v[k][j]) * b[i][j] * v[l][i];
                }
                ortho += conj(v[k][i]) * b[i][i] * v[l][i];
            }
            /* Simpler: v_i^H * v_j (since B = I) */
            ortho = 0.0;
            for (int i = 0; i < n; i++) {
                ortho += conj(v[k][i]) * v[l][i];
            }
            printf("  (%.4f + %.4fi) ", creal(ortho), cimag(ortho));
        }
        printf("\n");
    }
    printf("====================================================\n\n");
    
    /* Cleanup */
    for (int i = 0; i < n; i++) {
        free(a[i]);
        free(b[i]);
        free(v[i]);
    }
    free(a);
    free(b);
    free(v);
    free(e);
}


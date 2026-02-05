/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/debug.h>

#include <gemstore/model/model.h>
#include <gemstore/model/NRScreen.h>

#include <gemstore/basis/intrin.h>
#include <gemstore/basis/color.h>
#include <gemstore/basis/spin.h>
#include <gemstore/basis/isospin.h>
#include <gemstore/basis/orbit.h>
#include <gemstore/basis/soc.h>
#include <gemstore/basis/su3.h>

#include <gemstore/numerical/integral.h>
#include <gemstore/numerical/matrix.h>
#include <gemstore/numerical/eigen.h>

#include <gemstore/thread.h>

#include <stdio.h>
#include <stdlib.h>

void debug_su3_product()
{
    su3_product(1, 1, 1, 1);
}

void debug_soc_operator()
{
    double s1 = 0.5, s2 = 0.5, s = 1.0, l = 1.0;
    double s1p = 0.5, s2p = 0.5, sp = 1.0, lp = 1.0;
    double j = 2.0;
    printf("s1=%1.1f, s2=%1.1f, s=%1.1f, l=%1.1f\n", s1, s2, s, l);
    printf("s1p=%1.1f, s2p=%1.1f, sp=%1.1f, lp=%1.1f\n", s1p, s2p, sp, lp);
    printf("j=%1.1f\n", j);

    double sds = operator_sdots(0.5, 0.5, 1.0, 1.0, 0.5, 0.5, 1.0, 1.0);
    printf("SOC operator value (s1.s2): %f\n", sds);

    double lds1 = operator_sdots(1.0, 0.5, 1.0, 1.0, 1.0, 0.5, 1.0, 1.0);
    printf("SOC operator value (l.s1): %f\n", lds1);

    double lds2 = operator_sdots(1.0, 0.5, 1.0, 1.0, 0.5, 0.5, 1.0, 1.0);
    printf("SOC operator value (l.s2): %f\n", lds2);

    double tens = operator_tensor(0.5, 0.5, 1.0, 1.0, 0.5, 0.5, 1.0, 1.0, 2.0);
    printf("Tensor operator value: %f\n", tens);
}

void debug_color_wfn()
{
    intrin_wfn_t wf, ref_wf;

    wf = color_wfn_penta38();
    printf("ColorWFPenta38:\n");
    intrin_wfn_print(&wf);

    ref_wf = color_wfn_penta68();
    printf("\nColorWFPenta68:\n");
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
    overlap = integral_wfn_overlap(GRnlr_nonexp, factor, &args_bra, &args_bra);
    printf("Normalized overlap for Gr: %f\n", overlap);

    factor = 1 / sqrt(nu1 + nu2);
    overlap = integral_wfn_overlap(GRnlr_nonexp, factor, &args_bra, &args_ket);
    printf("Orthogonal overlap for Gr: %f\n", overlap);

    factor = sqrt(2 * nu1);
    overlap = integral_wfn_overlap_complex(GRnlp_nonexp, factor, &args_bra, &args_bra);
    printf("Normalized overlap for Gp: %f\n", overlap);

    factor = sqrt(4 * nu1 * nu2 / (nu1 + nu2));
    overlap = integral_wfn_overlap_complex(GRnlp_nonexp, factor, &args_bra, &args_ket);
    printf("Orthogonal overlap for Gp: %f\n", overlap);

    double beta1 = 0.8;
    args_bra.scale = beta1;
    double beta2 = 1.0;
    args_ket.scale = beta2;

    factor = 1 / beta1;
    overlap = integral_wfn_overlap(SRnlr_nonexp, factor, &args_bra, &args_bra);
    printf("Normalized overlap for Sr: %f\n", overlap);

    factor = sqrt(2 / (beta1 * beta1 + beta2 * beta2));
    overlap = integral_wfn_overlap(SRnlr_nonexp, factor, &args_bra, &args_ket);
    printf("Orthogonal overlap for Sr: %f\n", overlap);

    factor = beta1;
    overlap = integral_wfn_overlap_complex(SRnlp_nonexp, factor, &args_bra, &args_bra);
    printf("Normalized overlap for Sp: %f\n", overlap);

    factor = sqrt(2 * beta1 * beta1 * beta2 * beta2 / (beta1 * beta1 + beta2 * beta2));
    overlap = integral_wfn_overlap_complex(SRnlp_nonexp, factor, &args_bra, &args_ket);
    printf("Orthogonal overlap for Sp: %f\n", overlap);
}

void debug_potential_model()
{
    double s1 = 0.5, s2 = 0.5;
    double S = 0, L = 0, J = 0;

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

    argsModel_t args_model = argsNRScreen_meson;
    double m1 = args_model.mc;
    double m2 = args_model.mc;
    argsFlavor_t args_flavor = {
        .m1 = m1,
        .m2 = m2
    };

    double cen = operator_center(s1, s2, S, L, s1, s2, S, L);
    double sds = operator_sdots(s1, s2, S, L, s1, s2, S, L);
    double ls1 = operator_ldots1(s1, s2, S, L, s1, s2, S, L, J);
    double ls2 = operator_ldots2(s1, s2, S, L, s1, s2, S, L, J);
    double ten = operator_tensor(s1, s2, S, L, s1, s2, S, L, J);
    argsSOC_t args_soc = {
        .OCent = 1,
        .OSdS = sds,
        .OLS1 = ls1,
        .OLS2 = ls2,
        .OTens = ten
    };
    printf("OCen = %f, OSdS = %f, OLS1 = %f, OLS2 = %f, OTens = %f\n", cen, sds, ls1, ls2, ten);

    double factor = 1 / sqrt(nu1 + nu2);
    double factor_complex =  sqrt(4 * nu1 * nu2 / (nu1 + nu2));
    double element;

    element = integral_matrix_element_complex(GRnlp_nonexp, NRScreen_T, factor_complex, &args_bra, &args_ket, &args_flavor, &args_soc, &args_model);
    printf("Matrix element for <1|T|2>: %f\n", element);

    element = integral_matrix_element(GRnlr_nonexp, NRScreen_Vconf, factor, &args_bra, &args_ket, &args_flavor, &args_soc, &args_model);
    printf("Matrix element for <1|Vconf|2>: %f\n", element);

    element = integral_matrix_element(GRnlr_nonexp, NRScreen_Vcont, factor, &args_bra, &args_ket, &args_flavor, &args_soc, &args_model);
    printf("Matrix element for <1|Vcont|2>: %f\n", element);

    element = integral_matrix_element(GRnlr_nonexp, NRScreen_Vsocm, factor, &args_bra, &args_ket, &args_flavor, &args_soc, &args_model);
    printf("Matrix element for <1|Vsocm|2>: %f\n", element);

    element = integral_matrix_element(GRnlr_nonexp, NRScreen_Vsotp, factor, &args_bra, &args_ket, &args_flavor, &args_soc, &args_model);
    printf("Matrix element for <1|Vsotp|2>: %f\n", element);

    element = integral_matrix_element(GRnlr_nonexp, NRScreen_Vtens, factor, &args_bra, &args_ket, &args_flavor, &args_soc, &args_model);
    printf("Matrix element for <1|Vtens|2>: %f\n", element);
}

void debug_eigen_system()
{
    int n = 10;
    matrix_t mat = matrix_random(n, n);
    matrix_print(&mat);

    double *e1 = (double *)malloc(n * sizeof(double));
    double *e2 = (double *)malloc(n * sizeof(double));
    double **v = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        v[i] = (double *)malloc(n * sizeof(double));
    }

    eigen_standard_thread(mat.value, n, e1, e2, v, n);
    //eigen_general_thread(Hfi, Nfi, n, e1, e2, v, n, &info);

    array_t val = {
        .con = n,
        .value = e1
    };
    array_print(&val);

    matrix_t vec = {
        .row = n,
        .col = n,
        .value = v
    };
    matrix_print(&vec);

    for (int i = 0; i < n; i++) {
        free(v[i]);
    }
    free(v);
    free(e1);
    free(e2);
    matrix_free(&mat);
}
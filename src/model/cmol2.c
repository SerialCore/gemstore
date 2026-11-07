/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Two-body molecule GEM. Radial path as meson_nrmodel_GEM: getnu,
 * GRnlr / GRnlp, integral_nlr_hamilton / integral_nlp_hamilton, H c = E N c.
 * Partial waves are L and L+2 with every total spin that couples to J.
 * Tensor S12 mixes those blocks (Liu ²S–⁴D at J=1/2, ⁴S–²D–⁴D at J=3/2).
 */

#include <gemstore/model/cmol2.h>
#include <gemstore/model/obemodel.h>
#include <gemstore/model/ykwmodel.h>
#include <gemstore/model/eftmodel.h>

#include <gemstore/basis/orbit.h>
#include <gemstore/math/soc.h>
#include <gemstore/math/eigen.h>
#include <gemstore/math/matrix.h>
#include <gemstore/math/integral.h>
#include <gemstore/input.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int mol2_triangle(double a, double b, double c)
{
    return (a + b + 1e-10 >= c) && (a + c + 1e-10 >= b) && (b + c + 1e-10 >= a);
}

static void mol2_term_name(int L, double S, char *out, int n)
{
    const char *lett = "SPDFGHI";
    int mult = (int)(2.0 * S + 1.0 + 1e-8);
    char ell = (L >= 0 && L <= 6) ? lett[L] : '?';
    snprintf(out, (size_t)n, "%d%c", mult, ell);
}

/* Two t=1/2: ⟨τ1·τ2⟩ = 2I(I+1) − 3 = iso_pi. */
static double mol2_pair_I(double iso_pi)
{
    double x = 0.5 * (iso_pi + 3.0);
    if (x < -1e-12) {
        return -1.0;
    }
    return 0.5 * (-1.0 + sqrt(1.0 + 4.0 * x));
}

/* (−1)^{ℓ+(2s−S)+(2t−T)} = −1 fermion, +1 boson. */
static int mol2_pair_sym_ok(int ell, double s, double t, double S, double T, int fermion)
{
    int two_s = (int)(2.0 * s + 0.5);
    int two_t = (int)(2.0 * t + 0.5);
    int iS = (int)(S + 0.5);
    int iT = (int)(T + 0.5);
    int expo = ell + (two_s - iS) + (two_t - iT);
    int odd = (expo % 2 != 0);
    return fermion ? odd : !odd;
}

int mol2_build_channels(const argsInput_t *input, int *Lout, double *Sout, int maxch)
{
    double s1, s2, J, S, totI;
    int L0, L, n = 0;
    int ident, fermion;

    if (input == NULL || maxch < 1) {
        return 0;
    }

    s1 = obe_leg_spin(input->type1);
    s2 = obe_leg_spin(input->type2);
    J = input->J;
    L0 = (int)(input->L + 1e-8);
    if (L0 < 0) {
        L0 = 0;
    }

    ident = (input->type1 == input->type2 && fabs(input->m1 - input->m2) < 1e-9);
    fermion = (input->type1 == OBE_LEG_DIRAC);
    totI = mol2_pair_I(input->iso_pi);

    for (L = L0; L <= L0 + 2; L += 2) {
        for (S = fabs(s1 - s2); S <= s1 + s2 + 1e-10; S += 1.0) {
            int t_half, t_zero;
            if (!mol2_triangle(S, (double)L, J)) {
                continue;
            }
            /* Two I=1/2 legs: T comes from iso_pi. Two I=0 legs (Λc): T=0
             * and the isospin factor in the exchange sign is absent. */
            t_half = fabs(input->iso_t1 - 0.5) < 1e-9 && fabs(input->iso_t2 - 0.5) < 1e-9;
            t_zero = fabs(input->iso_t1) < 1e-9 && fabs(input->iso_t2) < 1e-9;
            if (ident && t_half && totI >= -1e-12
                && !mol2_pair_sym_ok(L, s1, 0.5, S, totI, fermion)) {
                continue;
            }
            if (ident && t_zero
                && !mol2_pair_sym_ok(L, s1, 0.0, S, 0.0, fermion)) {
                continue;
            }
            if (ident && !t_half && !t_zero) {
                fprintf(stderr, "Error: identical-leg exchange supports t=1/2 or t=0 only.\n");
                exit(1);
            }
            if (Lout) Lout[n] = L;
            if (Sout) Sout[n] = S;
            n++;
            if (n >= maxch) {
                return n;
            }
        }
    }
    return n;
}

int mol2_channel_count(const argsInput_t *input)
{
    return mol2_build_channels(input, NULL, NULL, MOL2_CHMAX);
}

static void mol2_radius_RMS(const argsOrbit_t *basis, int nch, int nmax, const matrix_t *Nfi,
    matrix_t *v_out, array_t *rms_out)
{
    int n, a, i, j, ntot;
    matrix_t mR2;

    if (rms_out == NULL || v_out == NULL || Nfi == NULL) {
        return;
    }

    ntot = nch * nmax;
    mR2 = matrix_init(ntot, ntot);
    for (a = 0; a < nch; a++) {
        for (i = 0; i < nmax; i++) {
            for (j = 0; j < nmax; j++) {
                int ia = a * nmax + i;
                int ja = a * nmax + j;
                double factor = 1.0 / sqrt(basis[ia].scale + basis[ja].scale);
                mR2.value[ia][ja] = integral_nlr_radius(GRnlr, factor, &basis[ia], &basis[ja]);
            }
        }
    }

    *rms_out = array_init(v_out->row);
    for (n = 0; n < v_out->row; n++) {
        double nrm = matrix_expect(v_out, n, Nfi);
        double r2 = matrix_expect(v_out, n, &mR2);
        if (nrm > 1e-16 && r2 > 0.0) {
            rms_out->value[n] = sqrt(r2 / nrm) / GEMSTORE_FM;
        }
        else {
            rms_out->value[n] = 0.0;
        }
    }
    matrix_free(&mR2);
}

void mol2_obestat_GEM(const argsInput_t *args_input, const argsOBEModel_t *args_model, argsOBEModelDy_t *args_dynmc, 
    array_t *e_out, matrix_t *v_out, matrix_t *n_out, array_t *rms_out)
{
    int nmax = args_input->nmax;
    double rmax = args_input->rmax;
    double rmin = args_input->rmin;
    double J = args_input->J;
    double s1 = obe_leg_spin(args_model->type1);
    double s2 = obe_leg_spin(args_model->type2);
    double kss = (args_model->type1 == OBE_LEG_DIRAC ? 2.0 : 1.0)
        * (args_model->type2 == OBE_LEG_DIRAC ? 2.0 : 1.0);
    int Lch[MOL2_CHMAX];
    double Sch[MOL2_CHMAX];
    int nch, ntot, a, b, i, j;
    argsOrbit_t *basis;
    matrix_t Hfi, Nfi;
    obe_pot_ctx_t ctx;
    char name[8];

    nch = mol2_build_channels(args_input, Lch, Sch, MOL2_CHMAX);
    if (nch < 1) {
        fprintf(stderr, "Error: no molecular partial waves couple to J=%.1f at L=%.0f",
            J, args_input->L);
        if (args_input->type1 == args_input->type2
            && fabs(args_input->m1 - args_input->m2) < 1e-9) {
            fprintf(stderr, " (identical-leg Fermi/Bose emptied the list)");
        }
        fprintf(stderr, ".\n");
        exit(1);
    }
    ntot = nch * nmax;
    if (e_out == NULL || e_out->value == NULL || e_out->len < ntot) {
        fprintf(stderr, "Error: mol2_obestat_GEM needs e_out of length %d.\n", ntot);
        exit(1);
    }
    if (v_out != NULL && (v_out->row < ntot || v_out->col < ntot)) {
        fprintf(stderr, "Error: mol2_obestat_GEM needs v_out of size %d.\n", ntot);
        exit(1);
    }

    printf("  Molecular channels (J=%.1f):", J);
    for (a = 0; a < nch; a++) {
        mol2_term_name(Lch[a], Sch[a], name, (int)sizeof(name));
        printf("  %s", name);
    }
    printf("   (%d x %d GEM)\n", nch, nmax);
    if (args_input->type1 == args_input->type2
        && fabs(args_input->m1 - args_input->m2) < 1e-9) {
        printf("  identical %s: (−1)^{L+(2s−S)+(2t−I)}=%s\n",
            args_input->type1 == OBE_LEG_DIRAC ? "fermions" : "bosons",
            args_input->type1 == OBE_LEG_DIRAC ? "−1" : "+1");
    }

    basis = (argsOrbit_t *)malloc((size_t)ntot * sizeof(argsOrbit_t));
    for (a = 0; a < nch; a++) {
        for (i = 0; i < nmax; i++) {
            int ia = a * nmax + i;
            basis[ia].n = i + 1;
            basis[ia].l = Lch[a];
            basis[ia].scale = getnu(i + 1, nmax, rmax, rmin);
        }
    }

    Hfi = matrix_init(ntot, ntot);
    Nfi = matrix_init(ntot, ntot);

    args_dynmc->mi = args_input->m1;
    args_dynmc->mj = args_input->m2;
    args_dynmc->OLSi = 0.0;
    args_dynmc->OLSj = 0.0;
    ctx.model = args_model;
    ctx.dyn = args_dynmc;

    for (a = 0; a < nch; a++) {
        for (b = 0; b < nch; b++) {
            int same = (Lch[a] == Lch[b] && fabs(Sch[a] - Sch[b]) < 1e-10);
            args_dynmc->OCent = operator_center_sl(s1, s2, Sch[a], Lch[a], s1, s2, Sch[b], Lch[b], J);
            args_dynmc->OSdS = kss * operator_sdots_sl(s1, s2, Sch[a], Lch[a], s1, s2, Sch[b], Lch[b], J);
            args_dynmc->OTens = kss * operator_tensor_sl(s1, s2, Sch[a], Lch[a], s1, s2, Sch[b], Lch[b], J);

            for (i = 0; i < nmax; i++) {
                for (j = 0; j < nmax; j++) {
                    int ia = a * nmax + i;
                    int jb = b * nmax + j;
                    double factor_r = 1.0 / sqrt(basis[ia].scale + basis[jb].scale);
                    double factor_p = sqrt(4.0 * basis[ia].scale * basis[jb].scale
                        / (basis[ia].scale + basis[jb].scale));
                    double h = 0.0;
                    double ov = 0.0;

                    if (same) {
                        h += integral_nlp_hamilton(GRnlp, OBEVt, factor_p, &basis[ia], &basis[jb], &ctx);
                        h += integral_nlr_hamilton(GRnlr, OBEVc, factor_r, &basis[ia], &basis[jb], &ctx);
                        h += integral_nlr_hamilton(GRnlr, OBEVss, factor_r, &basis[ia], &basis[jb], &ctx);
                        ov = integral_nlr_overlap(GRnlr, factor_r, &basis[ia], &basis[jb]);
                    }
                    if (fabs(args_dynmc->OTens) > 1e-14) {
                        h += integral_nlr_hamilton(GRnlr, OBEVtens, factor_r, &basis[ia], &basis[jb], &ctx);
                    }
                    Hfi.value[ia][jb] = h;
                    Nfi.value[ia][jb] = ov;
                }
            }
        }
    }

    matrix_symmetrize(&Hfi);
    matrix_symmetrize(&Nfi);

#ifdef LAPACKE
    lapack_general(Hfi.value, Nfi.value, ntot, e_out->value, (v_out == NULL) ? NULL : v_out->value, ntot);
#else
    eigen_general(Hfi.value, Nfi.value, ntot, e_out->value, (v_out == NULL) ? NULL : v_out->value, ntot);
#endif
    e_out->len = ntot;

    if (n_out != NULL) {
        *n_out = matrix_init(Nfi.row, Nfi.col);
        matrix_copy(&Nfi, n_out);
    }
    mol2_radius_RMS(basis, nch, nmax, &Nfi, v_out, rms_out);

    free(basis);
    matrix_free(&Hfi);
    matrix_free(&Nfi);
}

void mol2_yukawa_GEM(const argsInput_t *args_input, const argsYKWModel_t *args_model, argsYKWModelDy_t *args_dynmc,
    array_t *e_out, matrix_t *v_out, matrix_t *n_out, array_t *rms_out)
{
    int nmax = args_input->nmax;
    double rmax = args_input->rmax;
    double rmin = args_input->rmin;
    double J = args_input->J;
    double s1 = obe_leg_spin(args_input->type1);
    double s2 = obe_leg_spin(args_input->type2);
    double kss = (args_input->type1 == OBE_LEG_DIRAC ? 2.0 : 1.0)
        * (args_input->type2 == OBE_LEG_DIRAC ? 2.0 : 1.0);
    int Lch[MOL2_CHMAX];
    double Sch[MOL2_CHMAX];
    int nch, ntot, a, b, i, j;
    argsOrbit_t *basis;
    matrix_t Hfi, Nfi;
    ykw_pot_ctx_t ctx;
    char name[8];

    nch = mol2_build_channels(args_input, Lch, Sch, MOL2_CHMAX);
    if (nch < 1) {
        fprintf(stderr, "Error: no molecular partial waves couple to J=%.1f at L=%.0f",
            J, args_input->L);
        if (args_input->type1 == args_input->type2
            && fabs(args_input->m1 - args_input->m2) < 1e-9) {
            fprintf(stderr, " (identical-leg Fermi/Bose emptied the list)");
        }
        fprintf(stderr, ".\n");
        exit(1);
    }
    ntot = nch * nmax;
    if (e_out == NULL || e_out->value == NULL || e_out->len < ntot) {
        fprintf(stderr, "Error: mol2_yukawa_GEM needs e_out of length %d.\n", ntot);
        exit(1);
    }
    if (v_out != NULL && (v_out->row < ntot || v_out->col < ntot)) {
        fprintf(stderr, "Error: mol2_yukawa_GEM needs v_out of size %d.\n", ntot);
        exit(1);
    }

    printf("  Molecular channels (J=%.1f):", J);
    for (a = 0; a < nch; a++) {
        mol2_term_name(Lch[a], Sch[a], name, (int)sizeof(name));
        printf("  %s", name);
    }
    printf("   (%d x %d GEM, YUKAWA)\n", nch, nmax);

    basis = (argsOrbit_t *)malloc((size_t)ntot * sizeof(argsOrbit_t));
    for (a = 0; a < nch; a++) {
        for (i = 0; i < nmax; i++) {
            int ia = a * nmax + i;
            basis[ia].n = i + 1;
            basis[ia].l = Lch[a];
            basis[ia].scale = getnu(i + 1, nmax, rmax, rmin);
        }
    }

    Hfi = matrix_init(ntot, ntot);
    Nfi = matrix_init(ntot, ntot);

    args_dynmc->mi = args_input->m1;
    args_dynmc->mj = args_input->m2;
    ctx.model = args_model;
    ctx.dyn = args_dynmc;

    for (a = 0; a < nch; a++) {
        for (b = 0; b < nch; b++) {
            int same = (Lch[a] == Lch[b] && fabs(Sch[a] - Sch[b]) < 1e-10);
            args_dynmc->OCent = operator_center_sl(s1, s2, Sch[a], Lch[a], s1, s2, Sch[b], Lch[b], J);
            args_dynmc->OSdS = kss * operator_sdots_sl(s1, s2, Sch[a], Lch[a], s1, s2, Sch[b], Lch[b], J);

            for (i = 0; i < nmax; i++) {
                for (j = 0; j < nmax; j++) {
                    int ia = a * nmax + i;
                    int jb = b * nmax + j;
                    double factor_r = 1.0 / sqrt(basis[ia].scale + basis[jb].scale);
                    double factor_p = sqrt(4.0 * basis[ia].scale * basis[jb].scale
                        / (basis[ia].scale + basis[jb].scale));
                    double h = 0.0;
                    double ov = 0.0;

                    if (same) {
                        h += integral_nlp_hamilton(GRnlp, YKWVt, factor_p, &basis[ia], &basis[jb], &ctx);
                        h += integral_nlr_hamilton(GRnlr, YKWVc, factor_r, &basis[ia], &basis[jb], &ctx);
                        h += integral_nlr_hamilton(GRnlr, YKWVss, factor_r, &basis[ia], &basis[jb], &ctx);
                        ov = integral_nlr_overlap(GRnlr, factor_r, &basis[ia], &basis[jb]);
                    }
                    Hfi.value[ia][jb] = h;
                    Nfi.value[ia][jb] = ov;
                }
            }
        }
    }

    matrix_symmetrize(&Hfi);
    matrix_symmetrize(&Nfi);

#ifdef LAPACKE
    lapack_general(Hfi.value, Nfi.value, ntot, e_out->value, (v_out == NULL) ? NULL : v_out->value, ntot);
#else
    eigen_general(Hfi.value, Nfi.value, ntot, e_out->value, (v_out == NULL) ? NULL : v_out->value, ntot);
#endif
    e_out->len = ntot;

    if (n_out != NULL) {
        *n_out = matrix_init(Nfi.row, Nfi.col);
        matrix_copy(&Nfi, n_out);
    }
    mol2_radius_RMS(basis, nch, nmax, &Nfi, v_out, rms_out);

    free(basis);
    matrix_free(&Hfi);
    matrix_free(&Nfi);
}

void mol2_eftcont_GEM(const argsInput_t *args_input, const argsEFTModel_t *args_model, argsEFTModelDy_t *args_dynmc,
    array_t *e_out, matrix_t *v_out, matrix_t *n_out, array_t *rms_out)
{
    int nmax = args_input->nmax;
    double rmax = args_input->rmax;
    double rmin = args_input->rmin;
    double J = args_input->J;
    double s1 = obe_leg_spin(args_input->type1);
    double s2 = obe_leg_spin(args_input->type2);
    int Lch[MOL2_CHMAX];
    double Sch[MOL2_CHMAX];
    int nch, ntot, a, b, i, j;
    argsOrbit_t *basis;
    matrix_t Hfi, Nfi;
    eft_pot_ctx_t ctx;
    char name[8];

    nch = mol2_build_channels(args_input, Lch, Sch, MOL2_CHMAX);
    if (nch < 1) {
        fprintf(stderr, "Error: no molecular partial waves couple to J=%.1f at L=%.0f",
            J, args_input->L);
        if (args_input->type1 == args_input->type2
            && fabs(args_input->m1 - args_input->m2) < 1e-9) {
            fprintf(stderr, " (identical-leg Fermi/Bose emptied the list)");
        }
        fprintf(stderr, ".\n");
        exit(1);
    }
    ntot = nch * nmax;
    if (e_out == NULL || e_out->value == NULL || e_out->len < ntot) {
        fprintf(stderr, "Error: mol2_eftcont_GEM needs e_out of length %d.\n", ntot);
        exit(1);
    }
    if (v_out != NULL && (v_out->row < ntot || v_out->col < ntot)) {
        fprintf(stderr, "Error: mol2_eftcont_GEM needs v_out of size %d.\n", ntot);
        exit(1);
    }

    printf("  Molecular channels (J=%.1f):", J);
    for (a = 0; a < nch; a++) {
        mol2_term_name(Lch[a], Sch[a], name, (int)sizeof(name));
        printf("  %s", name);
    }
    printf("   (%d x %d GEM, EFTCONT)\n", nch, nmax);

    basis = (argsOrbit_t *)malloc((size_t)ntot * sizeof(argsOrbit_t));
    for (a = 0; a < nch; a++) {
        for (i = 0; i < nmax; i++) {
            int ia = a * nmax + i;
            basis[ia].n = i + 1;
            basis[ia].l = Lch[a];
            basis[ia].scale = getnu(i + 1, nmax, rmax, rmin);
        }
    }

    Hfi = matrix_init(ntot, ntot);
    Nfi = matrix_init(ntot, ntot);

    args_dynmc->mi = args_input->m1;
    args_dynmc->mj = args_input->m2;
    ctx.model = args_model;
    ctx.dyn = args_dynmc;

    for (a = 0; a < nch; a++) {
        for (b = 0; b < nch; b++) {
            int same = (Lch[a] == Lch[b] && fabs(Sch[a] - Sch[b]) < 1e-10);
            args_dynmc->OCent = operator_center_sl(s1, s2, Sch[a], Lch[a], s1, s2, Sch[b], Lch[b], J);

            for (i = 0; i < nmax; i++) {
                for (j = 0; j < nmax; j++) {
                    int ia = a * nmax + i;
                    int jb = b * nmax + j;
                    double factor_r = 1.0 / sqrt(basis[ia].scale + basis[jb].scale);
                    double factor_p = sqrt(4.0 * basis[ia].scale * basis[jb].scale
                        / (basis[ia].scale + basis[jb].scale));
                    double h = 0.0;
                    double ov = 0.0;

                    if (same) {
                        h += integral_nlp_hamilton(GRnlp, EFTVt, factor_p, &basis[ia], &basis[jb], &ctx);
                        h += integral_nlr_hamilton(GRnlr, EFTVc, factor_r, &basis[ia], &basis[jb], &ctx);
                        ov = integral_nlr_overlap(GRnlr, factor_r, &basis[ia], &basis[jb]);
                    }
                    Hfi.value[ia][jb] = h;
                    Nfi.value[ia][jb] = ov;
                }
            }
        }
    }

    matrix_symmetrize(&Hfi);
    matrix_symmetrize(&Nfi);

#ifdef LAPACKE
    lapack_general(Hfi.value, Nfi.value, ntot, e_out->value, (v_out == NULL) ? NULL : v_out->value, ntot);
#else
    eigen_general(Hfi.value, Nfi.value, ntot, e_out->value, (v_out == NULL) ? NULL : v_out->value, ntot);
#endif
    e_out->len = ntot;

    if (n_out != NULL) {
        *n_out = matrix_init(Nfi.row, Nfi.col);
        matrix_copy(&Nfi, n_out);
    }
    mol2_radius_RMS(basis, nch, nmax, &Nfi, v_out, rms_out);

    free(basis);
    matrix_free(&Hfi);
    matrix_free(&Nfi);
}

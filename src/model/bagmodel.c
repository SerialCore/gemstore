/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/model/bagmodel.h>
#include <gemstore/basis/cmisym.h>
#include <gemstore/basis/intrin.h>
#include <gemstore/math/cmi.h>
#include <gemstore/math/eigen.h>

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static const double BAG_R_MIN = 1.0;
static const double BAG_R_MAX = 10.0;
static const int BAG_NQUARK_MAX = 6;
static const int BAG_VARIANT_ITER = 80;
static const double BAG_ST_EPS = 1e-12;
static const double BAG_GOLDEN = 0.6180339887498948482;

/* Si(x) = ∫_0^x sin(t)/t dt. Odd series Si(x) = Σ (-1)^k x^{2k+1} / ((2k+1) (2k+1)!). */
static double bag_si(double x)
{
    long double z, term, sum;

    if (x < 0.0) return -bag_si(-x);
    if (x == 0.0) return 0.0;

    z = (long double)x;
    term = z;
    sum = term;
    for (int k = 0; k < 100; k++) {
        long double n = (long double)k;
        long double den = (2.0L * n + 2.0L) * (2.0L * n + 3.0L) * (2.0L * n + 3.0L);

        term *= -(z * z) * (2.0L * n + 1.0L) / den;
        sum += term;
        if (fabsl(term) <= 1.0e-20L * fabsl(sum)) break;
    }
    return (double)sum;
}

/* Color-magnetic overlap A(x_i, x_j) of two ground-state bag modes. */
static double bag_overlap(double xi, double xj)
{
    double sxi = sin(xi), sxj = sin(xj);
    double cxi = cos(xi), cxj = cos(xj);
    double yi = xi - sxi * cxi;
    double yj = xj - sxj * cxj;
    double di = xi * sxi * sxi - 1.5 * yi;
    double dj = xj * sxj * sxj - 1.5 * yj;
    double si_sum, num;

    if (di == 0.0 || dj == 0.0) return NAN;

    si_sum = 2.0 * xi * bag_si(2.0 * xi)
           + 2.0 * xj * bag_si(2.0 * xj)
           - (xi + xj) * bag_si(2.0 * (xi + xj))
           - (xi - xj) * bag_si(2.0 * (xi - xj));
    num = -1.5 * yi * yj - 2.0 * xi * xj * sxi * sxi * sxj * sxj + 0.5 * xi * xj * si_sum;
    return 1.0 + num / (di * dj);
}

/* S_ij = 3 α_s(R) μ_i μ_j A(x_i, x_j) / R^3. */
static double bag_sij(double mi, double xi, double mj, double xj, double radius)
{
    double alpha = bag_alphas(radius);
    double mui = bag_momentum(mi, xi, radius);
    double muj = bag_momentum(mj, xj, radius);

    return 3.0 * alpha * mui * muj * bag_overlap(xi, xj) / (radius * radius * radius);
}

static double bag_pair_binding(int left, int right, const argsBAGModel_t *model)
{
    int a = left < right ? left : right;
    int b = left < right ? right : left;

    if (a == 2 && b == 3) return model->Bcs;
    if (a == 2 && b == 4) return model->Bbs;
    if (a == 3 && b == 3) return model->Bcc;
    if (a == 3 && b == 4) return model->Bbc;
    if (a == 4 && b == 4) return model->Bbb;
    return 0.0;
}

double bag_frequency(double m, double r)
{
    double beta, low, high;

    /* g(x) = x cos x - sin x (1 - β - sqrt(β^2 + x^2)), β = mR.
       g(2.03) > 0, g(π) < 0 for β >= 0. Massless root x ≈ 2.0428. */
    beta = m * r;
    low = 2.03;
    high = M_PI;
    for (int iter = 0; iter < 60; iter++) {
        double x = 0.5 * (low + high);
        double energy = sqrt(beta * beta + x * x);
        double g = x * cos(x) - sin(x) * (1.0 - beta - energy);

        if (g == 0.0) return x;
        if (g > 0.0) low = x;
        else high = x;
    }
    return 0.5 * (low + high);
}

void bag_mass(const argsInput_t *input, const argsBAGModel_t *model, const matrix_t *color, const matrix_t *spin,
              int nprod, const int *icolor, const int *ispin, double radius, int state, bag_mit_t *mit)
{
    double quark_mass[BAG_NQUARK_MAX];
    double pair_s[BAG_NQUARK_MAX * (BAG_NQUARK_MAX - 1) / 2];
    double pair_b[BAG_NQUARK_MAX * (BAG_NQUARK_MAX - 1) / 2];
    int flavor[BAG_NQUARK_MAX];
    int nq[4];
    int ncolor = 1, nspin = 1, nbasis, npair, nf, p;
    matrix_t hamilton;
    double energy, result;

    mit->nvec = 0;
    for (int a = 0; a < CMI_COLOR_MAX * CMI_SPIN_MAX; a++) mit->coeff[a] = 0.0;

    nq[0] = input->nn;
    nq[1] = input->ns;
    nq[2] = input->nc;
    nq[3] = input->nb;
    nf = nq[0] + nq[1] + nq[2] + nq[3];

    /* bag energy without CMI/CEI interaction */
    energy = 0.0;
    for (int f = 0; f < 4; f++) {
        double m = getmq(f + 1, model);

        mit->x[f] = bag_frequency(m, radius);
        mit->w[f] = bag_kinematic(m, mit->x[f], radius);
        mit->mu[f] = bag_momentum(m, mit->x[f], radius);
        energy += nq[f] * mit->w[f];
    }
    mit->vac = bag_vacuum(radius, model->B, model->Z);
    mit->R = radius;
    energy += mit->vac;
    npair = nf * (nf - 1) / 2;

    if (color != NULL) {
        flavor[0] = input->f1;
        flavor[1] = input->f2;
        flavor[2] = input->f3;
        flavor[3] = input->f4;
        flavor[4] = input->f5;
        flavor[5] = input->f6;
        for (int q = 0; q < nf; q++) quark_mass[q] = getmq(flavor[q], model);
        ncolor = (int)(sqrt((double)color->row) + 0.5);
        nspin = (int)(sqrt((double)spin->row) + 0.5);
        p = 0;
        for (int i = 0; i < nf; i++) {
            for (int j = i + 1; j < nf; j++) {
                pair_s[p] = bag_sij(quark_mass[i], mit->x[flavor[i] - 1], quark_mass[j], mit->x[flavor[j] - 1], radius);
                pair_b[p] = bag_pair_binding(flavor[i], flavor[j], model);
                p++;
            }
        }
    }

    nbasis = (color != NULL) ? nprod : 1;

    /* hamiltonian for CMI/CEI interaction matrix */
    hamilton = matrix_init(nbasis, nbasis);
    for (int a = 0; a < nbasis; a++) {
        for (int b = 0; b < nbasis; b++) {
            double hij = (a == b) ? energy : 0.0;

            if (color != NULL) {
                int ca = icolor[a], cb = icolor[b];
                int sa = ispin[a], sb = ispin[b];

                for (p = 0; p < npair; p++) {
                    double cfac = color->value[ca * ncolor + cb][p];
                    double sfac = spin->value[sa * nspin + sb][p];

                    hij -= pair_s[p] * cfac * sfac;
                    if (sa == sb) hij += (-3.0 / 8.0) * pair_b[p] * cfac;
                }
            }
            hamilton.value[a][b] = hij;
        }
    }
    matrix_symmetrize(&hamilton);

    /* final eigen system */
    mit->nvec = nbasis;
    if (nbasis == 1) {
        result = hamilton.value[0][0];
        mit->coeff[0] = 1.0;
    }
    else {
        double *ev = (double *)malloc((size_t)nbasis * sizeof(double));
        matrix_t vt = matrix_init(nbasis, nbasis);

        eigen_standard(hamilton.value, nbasis, ev, vt.value, nbasis);
        result = ev[state];
        for (int a = 0; a < nbasis; a++) mit->coeff[a] = vt.value[state][a];
        matrix_free(&vt);
        free(ev);
    }
    matrix_free(&hamilton);
    mit->mass = result;
}

void bag_variant(const argsInput_t *input, const argsBAGModel_t *model, const matrix_t *color, const matrix_t *spin,
                 int nprod, const int *icolor, const int *ispin, int state, bag_mit_t *mit)
{
    bag_mit_t left_mit, right_mit;
    double low = BAG_R_MIN, high = BAG_R_MAX;
    double left, right, fl, fr;

    left = high - (high - low) * BAG_GOLDEN;
    right = low + (high - low) * BAG_GOLDEN;
    bag_mass(input, model, color, spin, nprod, icolor, ispin, left, state, &left_mit);
    bag_mass(input, model, color, spin, nprod, icolor, ispin, right, state, &right_mit);
    fl = left_mit.mass;
    fr = right_mit.mass;

    for (int iter = 0; iter < BAG_VARIANT_ITER; iter++) {
        if (fl < fr) {
            high = right;
            right = left;
            fr = fl;
            left = high - (high - low) * BAG_GOLDEN;
            bag_mass(input, model, color, spin, nprod, icolor, ispin, left, state, &left_mit);
            fl = left_mit.mass;
        }
        else {
            low = left;
            left = right;
            fl = fr;
            right = low + (high - low) * BAG_GOLDEN;
            bag_mass(input, model, color, spin, nprod, icolor, ispin, right, state, &right_mit);
            fr = right_mit.mass;
        }
    }
    bag_mass(input, model, color, spin, nprod, icolor, ispin, 0.5 * (low + high), state, mit);
}

void bag_compute(const argsInput_t *input, const argsBAGModel_t *model, bag_mit_t *state, int *nstate, char *order)
{
    intrin_wfn_t color_wfn[CMI_COLOR_MAX];
    intrin_wfn_t spin_wfn[CMI_SPIN_MAX];
    matrix_t color_fac = {0};
    matrix_t spin_fac = {0};
    char config[7] = {0};
    int nquark, npair, used[CMI_SPIN_MAX];

    *nstate = 0;
    order[0] = '\0';

    /* calculation without CMI/CEI interaction */
    if (input->basis == BASIS_NONE) {
        state[0].n[0] = input->nn;
        state[0].n[1] = input->ns;
        state[0].n[2] = input->nc;
        state[0].n[3] = input->nb;
        state[0].st = NAN;
        bag_variant(input, model, NULL, NULL, 0, NULL, NULL, 0, &state[0]);
        *nstate = 1;
        return;
    }

    nquark = input->nn + input->ns + input->nc + input->nb;
    npair = nquark * (nquark - 1) / 2;
    if (cmi_load_basis(input, color_wfn, spin_wfn) != 0) {
        fprintf(stderr, "Error: Failed to load the CMI basis for bag spectra.\n");
        exit(1);
    }

    int ok = (color_wfn[0].num_terms > 0 && color_wfn[0].num_configs == nquark);

    /* configuration of quark constituent for CMI */
    for (int k = 0; ok && k < nquark; k++) {
        unsigned char slot = (unsigned char)color_wfn[0].configs[0][k];

        if (islower(slot)) config[k] = 'q';
        else if (isupper(slot)) config[k] = 'Q';
        else ok = 0;
    }
    config[nquark] = '\0';
    for (int s = 0; ok && s < input->ncolor; s++) {
        if (color_wfn[s].num_configs != nquark) ok = 0;
        for (int t = 0; ok && t < color_wfn[s].num_terms; t++) {
            for (int k = 0; k < nquark; k++) {
                unsigned char slot = (unsigned char)color_wfn[s].configs[t][k];

                if ((config[k] == 'q' && !islower(slot)) || (config[k] == 'Q' && !isupper(slot))) ok = 0;
            }
        }
    }
    if (!ok) {
        for (int i = 0; i < input->ncolor; i++) intrin_wfn_free(&color_wfn[i]);
        for (int i = 0; i < input->nspin; i++) intrin_wfn_free(&spin_wfn[i]);
        fprintf(stderr, "Error: Color channels do not share one quark/antiquark order.\n");
        exit(1);
    }

    /* calculation of color and spin factor matrix in basis space */
    color_fac = matrix_init(input->ncolor * input->ncolor, npair);
    spin_fac = matrix_init(input->nspin * input->nspin, npair);
    operator_lambda2(color_wfn, input->ncolor, config, &color_fac);
    operator_sigma2(spin_wfn, input->nspin, &spin_fac);
    for (int i = 0; i < input->ncolor; i++) intrin_wfn_free(&color_wfn[i]);
    for (int i = 0; i < input->nspin; i++) intrin_wfn_free(&spin_wfn[i]);
    strcpy(order, config);

    /* construct color-spin subspace and perform bag-model calculation */
    for (int i = 0; i < input->nspin; i++) used[i] = 0;
    for (int i = 0; i < input->nspin; i++) {
        double st;
        int nprod = 0;
        int icolor[CMI_COLOR_MAX * CMI_SPIN_MAX];
        int ispin[CMI_COLOR_MAX * CMI_SPIN_MAX];
        int flavor[6] = {input->f1, input->f2, input->f3, input->f4, input->f5, input->f6};

        if (used[i]) continue;
        st = input->spin_ch[i].st;
        for (int j = i; j < input->nspin; j++) {
            if (fabs(input->spin_ch[j].st - st) > BAG_ST_EPS) continue;
            used[j] = 1;
            for (int c = 0; c < input->ncolor; c++) {
                if (nquark >= 5 && !cmi_color_spin_allowed(input->color_id[c], &input->spin_ch[j], nquark, flavor)) continue;
                icolor[nprod] = c;
                ispin[nprod] = j;
                nprod++;
            }
        }
        if (nprod < 1) {
            matrix_free(&spin_fac);
            matrix_free(&color_fac);
            fprintf(stderr, "Error: No Pauli-allowed color-spin product for st=%.1f.\n", st);
            exit(1);
        }
        for (int k = 0; k < nprod; k++) {
            bag_mit_t *mit = &state[*nstate];

            mit->n[0] = input->nn;
            mit->n[1] = input->ns;
            mit->n[2] = input->nc;
            mit->n[3] = input->nb;
            bag_variant(input, model, &color_fac, &spin_fac, nprod, icolor, ispin, k, mit);
            mit->st = st;
            for (int a = 0; a < mit->nvec; a++) {
                mit->color_id[a] = input->color_id[icolor[a]];
                mit->spin_ch[a] = input->spin_ch[ispin[a]];
            }
            (*nstate)++;
        }
    }
    matrix_free(&spin_fac);
    matrix_free(&color_fac);
}

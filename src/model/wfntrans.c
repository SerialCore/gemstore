/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/model/wfntrans.h>
#include <gemstore/basis/orbit.h>
#include <gemstore/math/matrix.h>
#include <gemstore/math/integral.h>
#include <gemstore/param/argset.h>

#include <math.h>
#include <stdlib.h>

double get_normalized_factor(const argsInput_t *input, const double *vector)
{
    if (input == NULL || vector == NULL) {
        return 1.0;
    }

    int L = (int)input->L;
    int nmax = input->nmax;
    double normalized = 1.0;
    double overlap_sum = 0.0;
    argsOrbit_t *basis = (argsOrbit_t *)malloc(nmax * sizeof(argsOrbit_t));

    if (basis == NULL) {
        return normalized;
    }

    for (int i = 0; i < nmax; i++) {
        basis[i].n = i + 1;
        basis[i].l = L;
        basis[i].scale = getnu(i + 1, nmax, input->rmax, input->rmin);
        basis[i].param = input->omega;
    }

    for (int i = 0; i < nmax; i++) {
        for (int j = 0; j < nmax; j++) {
            double overlap_ij = 0.0;
            double node_factor = 1.0 / sqrt(basis[i].scale + basis[j].scale);

            if (input->orbit == ORBIT_GEM) {
                overlap_ij = integral_nlr_overlap(GRnlr, node_factor, &basis[i], &basis[j]);
            }
            else if (input->orbit == ORBIT_CRG) {
                overlap_ij = integral_crg_overlap(CGRnlr, node_factor, &basis[i], &basis[j]);
            }
            else if (input->orbit == ORBIT_SHO) {
                overlap_ij = (i == j) ? 1.0 : 0.0;
            }

            overlap_sum += vector[i] * overlap_ij * vector[j];
        }
    }

    if (fabs(overlap_sum) > 1e-12) {
        normalized = sqrt(1.0 / overlap_sum);
    }

    free(basis);
    return normalized;
}

double get_state_wfn_value(const argsInput_t *input, const double *vector, double normalized, double r)
{
    if (input == NULL || vector == NULL) {
        return 0.0;
    }

    int L = (int)input->L;
    int nmax = input->nmax;
    double fm = 5.06773093854369882649;
    double rGeV = r * fm;
    double psi_r = 0.0;

    for (int n = 0; n < nmax; n++) {
        double basis_func = 0.0;
        double N = n + 1;

        if (input->orbit == ORBIT_GEM) {
            double nu = getnu(N, nmax, input->rmax, input->rmin);
            basis_func = GRnlr(rGeV, N, L, nu) * exp(-nu * rGeV * rGeV);
        }
        else if (input->orbit == ORBIT_CRG) {
            double nu = getnu(N, nmax, input->rmax, input->rmin);
            complex basis_func_complex = CGRnlr(rGeV, N, L, nu, input->omega) * exp(-nu * rGeV * rGeV);
            basis_func = creal(basis_func_complex);
        }
        else if (input->orbit == ORBIT_SHO) {
            basis_func = SRnlr(rGeV, n, L, input->beta) * exp(-0.5 * input->beta * input->beta * rGeV * rGeV);
        }

        psi_r += vector[n] * basis_func;
    }

    return psi_r * normalized;
}

void get_meson_rmsradii(const argsInput_t *input, const matrix_t *vector, array_t *radius, int len)
{
    if (input == NULL || vector == NULL || radius == NULL) {
        return;
    }

    int nmax = input->nmax;
    int L = (int)input->L;
    argsOrbit_t *basis = (argsOrbit_t *)malloc(nmax * sizeof(argsOrbit_t));
    for (int i = 0; i < nmax; i++) {
        basis[i].n = i + 1;
        basis[i].l = L;
        basis[i].scale = getnu(i + 1, nmax, input->rmax, input->rmin);
        basis[i].param = input->omega;
    }

    matrix_t mR2 = matrix_init(nmax, nmax);
    for (int i = 0; i < nmax; i++) {
        for (int j = 0; j < nmax; j++) {
            double factor = 1.0 / sqrt(basis[i].scale + basis[j].scale);

            if (input->orbit == ORBIT_GEM) {
                mR2.value[i][j] = integral_nlr_radius(GRnlr, factor, &basis[i], &basis[j]);
            }
            else if (input->orbit == ORBIT_CRG) {
                mR2.value[i][j] = integral_crg_radius(CGRnlr, factor, &basis[i], &basis[j]);
            }
            else {
                mR2.value[i][j] = (i == j) ? 1.0 : 0.0;
            }
        }
    }

    double fm = 5.06773093854369882649;
    for (int n = 0; n < len; n++) {
        double r2sum = 0.0;
        double normalized = get_normalized_factor(input, vector->value[n]);

        for (int i = 0; i < nmax; i++) {
            for (int j = 0; j < nmax; j++) {
                double coef = vector->value[n][i] * vector->value[n][j];
                r2sum += coef * mR2.value[i][j];
            }
        }

        radius->value[n] = sqrt(r2sum) * normalized / fm;
    }

    free(basis);
    matrix_free(&mR2);
}

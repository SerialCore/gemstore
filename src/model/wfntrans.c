/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Unit and normalization conventions used in this file:
 *
 * The coordinate-space meson wavefunction is treated as the radial function
 * R(r), not the reduced wavefunction u(r) = r R(r). In this convention,
 * normalization and radius moments are
 *
 *     integral dr r^2 |R(r)|^2 = 1,
 *     <r^2> = integral dr r^4 |R(r)|^2.
 *
 * This means get_normalized_factor() and radius_meson_rms() work with R(r)
 * together with the radial measure. The overlap integrals in
 * src/math/integral.c therefore include the extra r^2 factor, and the radius
 * integrals include r^4. This is mathematically equivalent to the reduced
 * wavefunction convention
 *
 *     u(r) = r R(r),
 *     integral dr |u(r)|^2 = 1,
 *     <r^2> = integral dr r^2 |u(r)|^2,
 *
 * but the code here does not build u(r) explicitly; it constructs and exports
 * R(r), and derives u(r) later only for output when needed.
 *
 * Internally, the orbital basis parameters are defined in natural units. The
 * basis scale nu returned by getnu() is in GeV^2, the radial argument used in
 * the Gaussian/CRG/SHO basis functions is in GeV^-1, and beta is combined with
 * the radius in the same natural-unit convention. For that reason,
 * get_state_wfn_value() converts the plotting radius from fm to GeV^-1 before
 * evaluating the basis.
 *
 * After the basis sum is formed, the overlap-based normalization factor
 * normalized = 1 / sqrt(c^T S c) is applied to the pointwise wavefunction.
 * The final return value is then converted back to an fm-based radial
 * wavefunction so that the normalization condition above holds with r measured
 * in fm. The exported quantities should therefore be interpreted as
 *
 *     R(r)            [fm^(-3/2)],
 *     u(r) = r R(r)   [fm^(-1/2)],
 *     r^2 |R(r)|^2    [fm^(-1)].
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

    return psi_r * normalized * pow(fm, 1.5);
}

void radius_meson_rms(const argsInput_t *input, const matrix_t *vector, array_t *radius, int len)
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

void effective_beta_sho(const argsInput_t *input, const array_t *radius, array_t *ebeta, int len)
{
    if (input == NULL || radius == NULL || ebeta == NULL) {
        return;
    }

    int L = (int)input->L;
    double fm = 5.06773093854369882649;

    for (int n = 0; n < len; n++) {
        double target = radius->value[n];

        if (target <= 0.0) {
            ebeta->value[n] = 0.0;
            continue;
        }

        ebeta->value[n] = sqrt(2.0 * n + L + 1.5) / (target * fm);
    }
}

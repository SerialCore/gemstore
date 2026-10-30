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
 * This means the meson RMS output and the overlap-normalized radial dump
 * work with R(r) together with the radial measure. The overlap integrals in
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
 * the Gaussian/SHO basis functions is in GeV^-1, and beta is combined with
 * the radius in the same natural-unit convention. For that reason,
 * get_gem_channel_wfn() and get_sho_channel_wfn() convert the plotting radius
 * from fm to GeV^-1 before evaluating the basis.
 *
 * After the basis sum is formed, the meson writer multiplies by
 * 1 / sqrt(c^T N c) from matrix_expect() on the solver overlap N. The channel
 * functions return an fm-based radial wavefunction so that the normalization
 * condition above holds with r measured in fm once that factor is applied.
 * The exported quantities should therefore be interpreted as
 *
 *     R(r)            [fm^(-3/2)],
 *     u(r) = r R(r)   [fm^(-1/2)],
 *     r^2 |R(r)|^2    [fm^(-1)].
 */

#include <gemstore/model/wfntrans.h>
#include <gemstore/basis/orbit.h>
#include <gemstore/input.h>

#include <math.h>
#include <stddef.h>

double get_gem_channel_wfn(int L, int nmax, double rmax, double rmin, const double *coeff, double r)
{
    int n;
    double fm = GEMSTORE_FM;
    double rGeV = r * fm;
    double psi_r = 0.0;

    if (coeff == NULL || nmax < 1) {
        return 0.0;
    }

    for (n = 0; n < nmax; n++) {
        double nu = getnu(n + 1, nmax, rmax, rmin);
        double basis_func = GRnlr(rGeV, n + 1, L, nu) * exp(-nu * rGeV * rGeV);
        psi_r += coeff[n] * basis_func;
    }
    return psi_r * pow(fm, 1.5);
}

double get_sho_channel_wfn(int L, int nmax, double beta, const double *coeff, double r)
{
    int n;
    double fm = GEMSTORE_FM;
    double rGeV = r * fm;
    double psi_r = 0.0;

    if (coeff == NULL || nmax < 1) {
        return 0.0;
    }

    for (n = 0; n < nmax; n++) {
        double basis_func = SRnlr(rGeV, n, L, beta) * exp(-0.5 * beta * beta * rGeV * rGeV);
        psi_r += coeff[n] * basis_func;
    }
    return psi_r * pow(fm, 1.5);
}

void get_effective_beta_sho(const argsInput_t *input, const array_t *radius, array_t *ebeta, int len)
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

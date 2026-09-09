/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MATH_INTEGRAL
#define GEMSTORE_MATH_INTEGRAL

#include <gemstore/basis/orbit.h>
#include <gemstore/param/argset.h>

/* Integrand for hamilton matrix elements. Model-specific parameters live in ctx. */
typedef double (*potential_t)(double x, void *ctx);

/* Integrate rms radius with NLR basis */
double integral_nlr_radius(
    orbit_nlr_t wfn,
    double node_factor,
    const argsOrbit_t *args_bra,
    const argsOrbit_t *args_ket);

/* Integrate rms radius with CRG basis */
double integral_crg_radius(
    orbit_crg_t wfn,
    double node_factor,
    const argsOrbit_t *args_bra,
    const argsOrbit_t *args_ket);

/* Integrate wavefunction overlaps with NLR basis */
double integral_nlr_overlap(
    orbit_nlr_t wfn,
    double node_factor,
    const argsOrbit_t *args_bra,
    const argsOrbit_t *args_ket);

/* Integrate wavefunction overlaps with NLP basis */
double integral_nlp_overlap(
    orbit_nlp_t wfn,
    double node_factor,
    const argsOrbit_t *args_bra,
    const argsOrbit_t *args_ket);

/* Integrate wavefunction overlaps with CRG basis */
double integral_crg_overlap(
    orbit_crg_t wfn,
    double node_factor,
    const argsOrbit_t *args_bra,
    const argsOrbit_t *args_ket);

/* Integrate hamiltonian with given potential and NLR basis */
double integral_nlr_hamilton(
    orbit_nlr_t wfn,
    potential_t pot,
    double node_factor,
    const argsOrbit_t *args_bra,
    const argsOrbit_t *args_ket,
    void *ctx);

/* Integrate hamiltonian with given potential and NLP basis */
double integral_nlp_hamilton(
    orbit_nlp_t wfn,
    potential_t pot,
    double node_factor,
    const argsOrbit_t *args_bra,
    const argsOrbit_t *args_ket,
    void *ctx);

/* Integrate hamiltonian with given potential and CRG basis */
double integral_crg_hamilton(
    orbit_crg_t wfn,
    potential_t pot,
    double node_factor,
    const argsOrbit_t *args_bra,
    const argsOrbit_t *args_ket,
    void *ctx);

/* ∫_0^∞ r^2 exp(-b11 r^2) pot(r) dr; Gaussian is in the quadrature weight */
double integral_exp_r2(potential_t pot, double b11, void *ctx);

/* ∫_0^∞ r^n exp(−b11 r²) pot(r) dr. Not a GRnlr integral: no ν^{l/2+3/4}. */
double integral_exp_rn(potential_t pot, double b11, int n, void *ctx);

#endif
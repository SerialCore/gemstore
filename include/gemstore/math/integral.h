/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MATH_INTEGRAL
#define GEMSTORE_MATH_INTEGRAL

#include <gemstore/basis/orbit.h>
#include <gemstore/model/gimodel.h>

/* Integrate wavefunction overlaps */
double integral_wfn_overlap(
    orbit_wfn_t wfn,
    double node_factor,
    const argsOrbit_t *args_bra,
    const argsOrbit_t *args_ket);

/* Integrate wavefunction overlaps (complex) */
double integral_wfn_overlap_complex(
    orbit_wfn_complex_t wfn,
    double node_factor,
    const argsOrbit_t *args_bra,
    const argsOrbit_t *args_ket);

/* Integrate matrix elements with given potential */
double integral_matrix_element(
    orbit_wfn_t wfn,
    potential_t pot,
    double node_factor,
    const argsOrbit_t *args_bra,
    const argsOrbit_t *args_ket,
    const argsGIModel_t *args_model,
    const argsGIModelDy_t *args_dynmc);

/* Integrate matrix elements with given potential (complex) */
double integral_matrix_element_complex(
    orbit_wfn_complex_t wfn,
    potential_t pot,
    double node_factor,
    const argsOrbit_t *args_bra,
    const argsOrbit_t *args_ket,
    const argsGIModel_t *args_model,
    const argsGIModelDy_t *args_dynmc);

#endif
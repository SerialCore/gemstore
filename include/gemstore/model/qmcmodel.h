/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_MODEL_QMCMODEL
#define GEMSTORE_MODEL_QMCMODEL

#include <gemstore/input.h>
#include <gemstore/param/argset.h>

/* Symmetric nuclear matter plus the bag that carries the spectrum.
 * The source is three light quarks in a color-free bag. Only those quarks
 * couple to σ and ω. QMCSTAD keeps B fixed. QMCMOD uses
 * B/B0 = (1 - 4 g_sigma_B σ / (δ M_N))^δ, with M_N the free source mass.
 * Temperature enters the nucleon Fermi–Dirac distributions only. */
typedef struct qmc_matter {
    argsQMCModelDy_t field;
    double m_free;          /* source mass at σ = 0, GeV */
    double m_nucleon;       /* in-medium source mass, GeV */
    double R_nucleon;       /* in-medium source radius, GeV^{-1} */
    double rho_s;           /* scalar density, GeV^3 */
    double mu;              /* baryon chemical potential, GeV */
    double eps;             /* energy density, GeV^4 */
    double press;           /* pressure, GeV^4 */
    double ebind;           /* ε/ρ_B - m_free, GeV; 0 at ρ_B = 0 */
} qmc_matter_t;

/* Solve σ and ω at input->rho_b and input->temp. On return, bag has the
 * in-medium light mass and bag constant used for the spectrum. */
void qmc_solve(const argsInput_t *input, const argsQMCModel_t *model, qmc_matter_t *matter, argsBAGModel_t *bag);

#endif

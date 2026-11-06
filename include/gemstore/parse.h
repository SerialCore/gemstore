/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_PARSE
#define GEMSTORE_PARSE

#include <gemstore/input.h>
#include <gemstore/param/argset.h>

/* parse input file */
void parse_input_file(const char *filename, argsInput_t *input);

/* parse GI parameters. mu is optional and defaults to 0 (GISTRING). */
void parse_param_GI(const char *filename, argsGIModel_t *args_model);

/* parse NR parameters. mu is optional and defaults to 0 (NRSTRING). */
void parse_param_NR(const char *filename, argsNRModel_t *args_model);

/* parse MITBAG parameters */
void parse_param_BAG(const char *filename, argsBAGModel_t *args_model);

/* parse QMC param file. QMCSTAD and QMCMOD share this schema. */
void parse_param_QMC(const char *filename, argsQMCModel_t *args_model);

/* parse OBESTAT param file: lambda, tensor_mass1/2, g_*, f_*.
 * Legs, cc, eta, iso_pi, iso_rho, m1/m2 stay in the system JSON. */
void parse_param_OBESTAT(const char *filename, argsOBEModel_t *args_model);

/* parse YUKAWA param file: param.terms[] of V (MeV·fm), mu (fm^{-1}), wigner, bartlett. */
void parse_param_YUKAWA(const char *filename, argsYKWModel_t *args_model);

#endif

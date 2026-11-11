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

/* parse QMC param file. g_sigma_B and delta are read only for QMCMOD. */
void parse_param_QMC(const char *filename, argsQMCModel_t *args_model, model_type_t model);

/* parse OBESTAT param file: lambda, tensor_mass1/2, g_*, f_*,
 * and required cc1, cc2, eta. Legs, masses, and I stay in the system JSON.
 * Pair ⟨τ·τ⟩ is set from I (two-body t=1/2, or three-body recoupling). */
void parse_param_OBESTAT(const char *filename, argsOBEModel_t *args_model);

/* parse YUKAWA param file: param.terms[] of V (MeV·fm), mu (fm^{-1}), wigner, bartlett. */
void parse_param_YUKAWA(const char *filename, argsYKWModel_t *args_model);

/* parse EFTCONT param file: param.Ca (MeV), param.Rc (fm). */
void parse_param_EFTCONT(const char *filename, argsEFTModel_t *args_model);

/* parse HALQCD param file. Flat form is param.V0, param.Vsig, param.VT.
 * Isospin form keys those objects by t_ij: "0", "1", "0.5", "1.5".
 * Each component is {r: fm, V: MeV}, strictly increasing in r. A missing array is zero. */
void parse_param_HALQCD(const char *filename, argsHALSet_t *set);

/* parse QPC3P0 param file: param.gamma, param.mq1, param.mq2, param.mq3 (GeV). */
void parse_param_QPC3P0(const char *filename, argsQPCModel_t *args_model);

#endif

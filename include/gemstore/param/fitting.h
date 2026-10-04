/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_PARAM_FITTING
#define GEMSTORE_PARAM_FITTING

#include <gemstore/param/argset.h>

/* Experimental masses and errors are MeV. The spectrum kernel returns GeV. */
typedef struct fit_state {
    const char *name;
    int f1;
    int f2;
    int n;                  /* 1-based radial index into eigenvalue[n - 1] */
    double S;
    double L;
    double J;
    double mass;
    double error;
} fit_state_t;

/* One row of the Minuit parameter vector, in Add() order. */
typedef struct fit_param {
    const char *name;
    double value;
    double step;
    double min;
    double max;
    int fixed;
} fit_param_t;

typedef struct fit_target {
    model_type_t model;
    const fit_state_t *states;
    int nstates;
    const fit_param_t *params;
    int nparams;
} fit_target_t;

typedef struct fit_result {
    double chi2;
    double edm;
    int valid;
    double *value;
    double *error;
    double *mass_mev;
    double *residual;
} fit_result_t;

void fitting_run(const argsInput_t *input);

#ifdef __cplusplus
extern "C" {
#endif

/* theta[i] matches target->params[i]. Returns the mass in GeV. */
double fit_predict_mass(const argsInput_t *input, const fit_target_t *target,
    const fit_state_t *state, const double *theta);

/* Weighted chi2 in MeV. mass_mev and residual may be NULL. */
double fit_chi2(const argsInput_t *input, const fit_target_t *target,
    const double *theta, double *mass_mev, double *residual);

void fit_minuit_run(const argsInput_t *input, const fit_target_t *target, fit_result_t *out);

#ifdef __cplusplus
}
#endif

extern const fit_target_t fit_target_giscreen_meson;
extern const fit_target_t fit_target_giscreen_ccbar;
extern const fit_target_t fit_target_giscreen_bbbar;
extern const fit_target_t fit_target_giscreen_bcbar;
extern const fit_target_t fit_target_giscreen_bsbar;
extern const fit_target_t fit_target_giscreen_csbar;

#endif

/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/param/fitting.h>
#include <gemstore/model/cmeson.h>
#include <gemstore/math/matrix.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const char *const k_gi_names[] = {
    "mn", "ms", "mc", "mb", "b", "mu", "c",
    "sigma_0", "s", "epsilon_cont", "epsilon_sov", "epsilon_sos", "epsilon_tens"
};

static const char *const k_nr_names[] = {
    "mn", "ms", "mc", "mb", "b", "mu", "c", "alpha_s", "sigma"
};

static const fit_target_t *const k_targets[] = {
    [FITTING_GISCREEN_MESON] = &fit_target_giscreen_meson,
    [FITTING_GISCREEN_BBBAR] = &fit_target_giscreen_bbbar,
    [FITTING_GISCREEN_BCBAR] = &fit_target_giscreen_bcbar,
    [FITTING_GISCREEN_BSBAR] = &fit_target_giscreen_bsbar,
    [FITTING_GISCREEN_CCBAR] = &fit_target_giscreen_ccbar,
    [FITTING_GISCREEN_CSBAR] = &fit_target_giscreen_csbar
};

static int set_gi_field(argsGIModel_t *model, const char *name, double value)
{
    if (strcmp(name, "mn") == 0) model->mn = value;
    else if (strcmp(name, "ms") == 0) model->ms = value;
    else if (strcmp(name, "mc") == 0) model->mc = value;
    else if (strcmp(name, "mb") == 0) model->mb = value;
    else if (strcmp(name, "b") == 0) model->b = value;
    else if (strcmp(name, "mu") == 0) model->mu = value;
    else if (strcmp(name, "c") == 0) model->c = value;
    else if (strcmp(name, "sigma_0") == 0) model->sigma_0 = value;
    else if (strcmp(name, "s") == 0) model->s = value;
    else if (strcmp(name, "epsilon_cont") == 0) model->epsilon_cont = value;
    else if (strcmp(name, "epsilon_sov") == 0) model->epsilon_sov = value;
    else if (strcmp(name, "epsilon_sos") == 0) model->epsilon_sos = value;
    else if (strcmp(name, "epsilon_tens") == 0) model->epsilon_tens = value;
    else return -1;
    return 0;
}

static int set_nr_field(argsNRModel_t *model, const char *name, double value)
{
    if (strcmp(name, "mn") == 0) model->mn = value;
    else if (strcmp(name, "ms") == 0) model->ms = value;
    else if (strcmp(name, "mc") == 0) model->mc = value;
    else if (strcmp(name, "mb") == 0) model->mb = value;
    else if (strcmp(name, "b") == 0) model->b = value;
    else if (strcmp(name, "mu") == 0) model->mu = value;
    else if (strcmp(name, "c") == 0) model->c = value;
    else if (strcmp(name, "alpha_s") == 0) model->alpha_s = value;
    else if (strcmp(name, "sigma") == 0) model->sigma = value;
    else return -1;
    return 0;
}

static void fit_params_check(model_type_t model, const fit_param_t *params, int nparams)
{
    const char *const *names;
    int nnames;
    int seen[16];
    int i;
    int j;

    if (model == MODEL_GISTRING || model == MODEL_GISCREEN) {
        names = k_gi_names;
        nnames = (int)(sizeof k_gi_names / sizeof k_gi_names[0]);
    }
    else if (model == MODEL_NRSTRING || model == MODEL_NRSCREEN) {
        names = k_nr_names;
        nnames = (int)(sizeof k_nr_names / sizeof k_nr_names[0]);
    }
    else {
        fprintf(stderr, "FITTING does not support this model\n");
        exit(1);
    }

    if (nparams != nnames) {
        fprintf(stderr, "FITTING parameter list for %s must contain %d parameters\n",
            model_type_str[model], nnames);
        exit(1);
    }

    memset(seen, 0, sizeof seen);
    for (i = 0; i < nparams; i++) {
        int found = -1;

        if (params[i].name == NULL) {
            fprintf(stderr, "Fit parameter %d has no name\n", i);
            exit(1);
        }
        for (j = 0; j < nnames; j++) {
            if (strcmp(params[i].name, names[j]) == 0) {
                found = j;
                break;
            }
        }
        if (found < 0) {
            fprintf(stderr, "Unknown parameter for %s: %s\n", model_type_str[model], params[i].name);
            exit(1);
        }
        if (seen[found]) {
            fprintf(stderr, "Duplicate parameter: %s\n", params[i].name);
            exit(1);
        }
        seen[found] = 1;
        if (!(params[i].min < params[i].max)) {
            fprintf(stderr, "Parameter %s has min >= max\n", params[i].name);
            exit(1);
        }
        if (params[i].value < params[i].min || params[i].value > params[i].max) {
            fprintf(stderr, "Parameter %s start is outside [min, max]\n", params[i].name);
            exit(1);
        }
        if (!params[i].fixed && !(params[i].step > 0.0)) {
            fprintf(stderr, "Parameter %s needs a positive step\n", params[i].name);
            exit(1);
        }
        if ((model == MODEL_GISTRING || model == MODEL_NRSTRING)
            && strcmp(params[i].name, "mu") == 0 && !params[i].fixed) {
            fprintf(stderr, "%s does not use mu. Keep mu fixed, or choose GISCREEN or NRSCREEN.\n",
                model_type_str[model]);
            exit(1);
        }
    }
}

static void fit_states_check(const argsInput_t *input, const fit_state_t *states, int nstates)
{
    int i;

    if (nstates < 1) {
        fprintf(stderr, "Fit data has no states\n");
        exit(1);
    }
    if (input->nmax < 1) {
        fprintf(stderr, "FITTING needs nmax >= 1\n");
        exit(1);
    }
    for (i = 0; i < nstates; i++) {
        if (states[i].n < 1 || states[i].n > input->nmax) {
            fprintf(stderr, "State %s radial index N=%d is outside 1..nmax (%d)\n",
                states[i].name != NULL ? states[i].name : "", states[i].n, input->nmax);
            exit(1);
        }
    }
}

double fit_predict_mass(const argsInput_t *input, const fit_target_t *target,
    const fit_state_t *state, const double *theta)
{
    argsInput_t one = *input;
    array_t eigenvalue;
    double mass;
    int i;

    one.f1 = state->f1;
    one.f2 = state->f2;
    one.S = state->S;
    one.L = state->L;
    one.J = state->J;
    eigenvalue = array_init(one.nmax);

    if (one.model == MODEL_GISTRING || one.model == MODEL_GISCREEN) {
        argsGIModel_t model;
        argsGIModelDy_t dyn = {0};

        memset(&model, 0, sizeof model);
        for (i = 0; i < target->nparams; i++) {
            if (set_gi_field(&model, target->params[i].name, theta[i]) != 0) {
                fprintf(stderr, "Unknown GI parameter: %s\n", target->params[i].name);
                exit(1);
            }
        }
        dyn.model = one.model;
        dyn.system = one.system;
        if (one.orbit == ORBIT_GEM) {
            meson_gimodel_GEM(&one, &model, &dyn, &eigenvalue, NULL, 0);
        }
        else if (one.orbit == ORBIT_SHO) {
            meson_gimodel_SHO(&one, &model, &dyn, &eigenvalue, NULL, 0);
        }
        else {
            fprintf(stderr, "FITTING meson supports GEM and SHO\n");
            exit(1);
        }
    }
    else if (one.model == MODEL_NRSTRING || one.model == MODEL_NRSCREEN) {
        argsNRModel_t model;
        argsNRModelDy_t dyn = {0};

        memset(&model, 0, sizeof model);
        for (i = 0; i < target->nparams; i++) {
            if (set_nr_field(&model, target->params[i].name, theta[i]) != 0) {
                fprintf(stderr, "Unknown NR parameter: %s\n", target->params[i].name);
                exit(1);
            }
        }
        dyn.model = one.model;
        dyn.system = one.system;
        if (one.orbit == ORBIT_GEM) {
            meson_nrmodel_GEM(&one, &model, &dyn, &eigenvalue, NULL, 0);
        }
        else if (one.orbit == ORBIT_SHO) {
            meson_nrmodel_SHO(&one, &model, &dyn, &eigenvalue, NULL, 0);
        }
        else {
            fprintf(stderr, "FITTING meson supports GEM and SHO\n");
            exit(1);
        }
    }
    else {
        fprintf(stderr, "FITTING does not support this model\n");
        exit(1);
    }

    mass = eigenvalue.value[state->n - 1];
    array_free(&eigenvalue);
    return mass;
}

double fit_chi2(const argsInput_t *input, const fit_target_t *target,
    const double *theta, double *mass_mev, double *residual)
{
    static int ncall = 0;
    double chi2 = 0.0;
    int i;

    for (i = 0; i < target->nstates; i++) {
        double calc = 1000.0 * fit_predict_mass(input, target, &target->states[i], theta);
        double diff = calc - target->states[i].mass;
        double weight = (target->states[i].error > 0.0)
            ? 1.0 / (target->states[i].error * target->states[i].error)
            : 1.0;

        chi2 += diff * diff * weight;
        if (mass_mev != NULL) {
            mass_mev[i] = calc;
        }
        if (residual != NULL) {
            residual[i] = diff;
        }
    }
    ncall += 1;
    printf("call %d  chi2 = %.8g\n", ncall, chi2);
    fflush(stdout);
    return chi2;
}

static void write_fit_result(const argsInput_t *input, const fit_target_t *target, const fit_result_t *result)
{
    char path[512];
    FILE *fp;
    int nw;
    int i;

    nw = snprintf(path, sizeof path, "%s.fit.json", input->project);
    if (nw < 0 || nw >= (int)sizeof path) {
        fprintf(stderr, "Fit output path is too long\n");
        exit(1);
    }
    fp = fopen(path, "w");
    if (fp == NULL) {
        fprintf(stderr, "Cannot open %s\n", path);
        exit(1);
    }

    fprintf(fp, "{\n");
    fprintf(fp, "  \"target\": \"%s\",\n", fitting_type_str[input->target]);
    fprintf(fp, "  \"model\": \"%s\",\n", model_type_str[input->model]);
    fprintf(fp, "  \"system\": \"%s\",\n", system_type_str[input->system]);
    fprintf(fp, "  \"basis\": \"%s\",\n", orbit_type_str[input->orbit]);
    fprintf(fp, "  \"mass_unit\": \"MeV\",\n");
    fprintf(fp, "  \"chi2\": %.16g,\n", result->chi2);
    fprintf(fp, "  \"edm\": %.16g,\n", result->edm);
    fprintf(fp, "  \"valid\": %s,\n", result->valid ? "true" : "false");
    fprintf(fp, "  \"parameters\": [\n");
    for (i = 0; i < target->nparams; i++) {
        fprintf(fp, "    {\"name\": \"%s\", \"value\": %.16g, \"error\": %.16g, \"fixed\": %s}%s\n",
            target->params[i].name, result->value[i], result->error[i],
            target->params[i].fixed ? "true" : "false",
            (i + 1 == target->nparams) ? "" : ",");
    }
    fprintf(fp, "  ],\n");
    fprintf(fp, "  \"states\": [\n");
    for (i = 0; i < target->nstates; i++) {
        const char *name = target->states[i].name != NULL ? target->states[i].name : "";
        fprintf(fp, "    {\"name\": \"%s\", \"mass_calc\": %.16g, \"mass_exp\": %.16g, \"residual\": %.16g}%s\n",
            name, result->mass_mev[i], target->states[i].mass, result->residual[i],
            (i + 1 == target->nstates) ? "" : ",");
    }
    fprintf(fp, "  ]\n}\n");
    fclose(fp);
    printf("Wrote %s\n", path);
}

void fitting_run(const argsInput_t *input)
{
    const fit_target_t *target;
    fit_result_t result = {0};

    if (input->system != SYSTEM_MESON) {
        fprintf(stderr, "FITTING currently supports only MESON\n");
        exit(1);
    }
    if (input->orbit != ORBIT_GEM && input->orbit != ORBIT_SHO) {
        fprintf(stderr, "FITTING meson supports GEM and SHO\n");
        exit(1);
    }

    target = k_targets[input->target];
    if (target->model != input->model) {
        fprintf(stderr, "Fit target %s is defined for %s, not %s\n",
            fitting_type_str[input->target], model_type_str[target->model], model_type_str[input->model]);
        exit(1);
    }

    fit_params_check(input->model, target->params, target->nparams);
    fit_states_check(input, target->states, target->nstates);
    fit_minuit_run(input, target, &result);
    printf("chi2 = %.8g  valid = %s  edm = %.8g\n",
        result.chi2, result.valid ? "true" : "false", result.edm);
    write_fit_result(input, target, &result);

    if (result.value == NULL) {
        return;
    }
    free(result.value);
    free(result.error);
    free(result.mass_mev);
    free(result.residual);
}

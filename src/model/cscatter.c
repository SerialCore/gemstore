/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Two-body S-wave scattering length. Projects one partial wave onto a
 * central potential and calls the effective-range solver.
 */

#include <gemstore/model/cscatter.h>
#include <gemstore/model/cmol2.h>
#include <gemstore/model/eremodel.h>

#include <gemstore/basis/orbit.h>
#include <gemstore/math/soc.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Projected central potential in MeV. r_fm is fm. */
static double scatter_job_v(double r_fm, void *ctx)
{
    scatter_job_t *job = (scatter_job_t *)ctx;
    double r = r_fm * GEMSTORE_FM;
    double gev;

    if (job->model == MODEL_EFTCONT) {
        gev = EFTVc(r, &job->eft);
    }
    else if (job->model == MODEL_HALQCD) {
        gev = HALVc(r, &job->hal) + HALVss(r, &job->hal) + HALVten(r, &job->hal);
    }
    else if (job->model == MODEL_OBESTAT) {
        gev = OBEVc(r, &job->obe) + OBEVss(r, &job->obe) + OBEVtens(r, &job->obe);
    }
    else {
        gev = YKWVc(r, &job->ykw) + YKWVss(r, &job->ykw);
    }
    return 1000.0 * gev;
}

void scatter_ere_mol2(const argsInput_t *input, double *mu_gev, double *a0_fm, double *r0_fm)
{
    int Lch[MOL2_CHMAX];
    double Sch[MOL2_CHMAX];
    int nch, a, found, extra;
    double s1, s2, kss, S, J;
    double mu;
    ere_t ere;
    scatter_job_t job;
    argsEFTModel_t eft_model;
    argsEFTModelDy_t eft_dyn = {0};
    argsHALModel_t hal_model;
    argsHALModelDy_t hal_dyn = {0};
    argsOBEModel_t obe_model;
    argsOBEModelDy_t obe_dyn = {0};
    argsYKWModel_t ykw_model;
    argsYKWModelDy_t ykw_dyn = {0};
    int rc;

    if (input == NULL || mu_gev == NULL || a0_fm == NULL || r0_fm == NULL) {
        fprintf(stderr, "Error: cscatter needs an input and result slots.\n");
        exit(1);
    }
    if (input->model != MODEL_OBESTAT && input->model != MODEL_YUKAWA
        && input->model != MODEL_EFTCONT && input->model != MODEL_HALQCD) {
        fprintf(stderr, "Error: SCATTER supports OBESTAT, YUKAWA, EFTCONT, and HALQCD.\n");
        exit(1);
    }
    if (input->m1 <= 0.0 || input->m2 <= 0.0) {
        fprintf(stderr, "Error: SCATTER needs positive physical masses m1 and m2.\n");
        exit(1);
    }
    if (fabs(input->L) > 1e-8) {
        fprintf(stderr, "Error: SCATTER solves the l=0 radial equation. Set L to 0.\n");
        exit(1);
    }

    nch = mol2_build_channels(input, Lch, Sch, MOL2_CHMAX);
    found = 0;
    extra = 0;
    for (a = 0; a < nch; a++) {
        if (Lch[a] == 0 && fabs(Sch[a] - input->S) < 1e-8) {
            found = 1;
        }
        else {
            extra = 1;
        }
    }
    if (!found) {
        fprintf(stderr, "Error: SCATTER has no l=0 channel for S=%.1f and J=%.1f",
            input->S, input->J);
        if (input->type1 == input->type2 && fabs(input->m1 - input->m2) < 1e-9) {
            fprintf(stderr, " (identical-leg Fermi/Bose emptied the list)");
        }
        fprintf(stderr, ".\n");
        exit(1);
    }
    if (extra) {
        printf("  SCATTER uses the diagonal l=0 potential. Coupling to L=2 is not solved.\n");
    }

    s1 = obe_leg_spin(input->type1);
    s2 = obe_leg_spin(input->type2);
    S = input->S;
    J = input->J;
    kss = (input->type1 == OBE_LEG_DIRAC ? 2.0 : 1.0)
        * (input->type2 == OBE_LEG_DIRAC ? 2.0 : 1.0);

    memset(&job, 0, sizeof job);
    job.model = input->model;
    if (input->model == MODEL_EFTCONT) {
        eft_model = argsEFTModel_from(input, 0);
        eft_dyn.OCent = operator_center_sl(s1, s2, S, 0.0, s1, s2, S, 0.0, J);
        job.eft.model = &eft_model;
        job.eft.dyn = &eft_dyn;
    }
    else if (input->model == MODEL_HALQCD) {
        hal_model = argsHALModel_from(input, 0);
        hal_dyn.OCent = operator_center_sl(s1, s2, S, 0.0, s1, s2, S, 0.0, J);
        hal_dyn.OSdS = kss * operator_sdots_sl(s1, s2, S, 0.0, s1, s2, S, 0.0, J);
        hal_dyn.OTens = kss * operator_tensor_sl(s1, s2, S, 0.0, s1, s2, S, 0.0, J);
        job.hal.model = &hal_model;
        job.hal.dyn = &hal_dyn;
    }
    else if (input->model == MODEL_OBESTAT) {
        obe_model = argsOBEModel_from(input, 0);
        obe_dyn.mi = input->m1;
        obe_dyn.mj = input->m2;
        obe_dyn.OCent = operator_center_sl(s1, s2, S, 0.0, s1, s2, S, 0.0, J);
        obe_dyn.OSdS = kss * operator_sdots_sl(s1, s2, S, 0.0, s1, s2, S, 0.0, J);
        obe_dyn.OTens = kss * operator_tensor_sl(s1, s2, S, 0.0, s1, s2, S, 0.0, J);
        job.obe.model = &obe_model;
        job.obe.dyn = &obe_dyn;
    }
    else {
        ykw_model = argsYKWModel_from(input, 0);
        ykw_dyn.mi = input->m1;
        ykw_dyn.mj = input->m2;
        ykw_dyn.OCent = operator_center_sl(s1, s2, S, 0.0, s1, s2, S, 0.0, J);
        ykw_dyn.OSdS = kss * operator_sdots_sl(s1, s2, S, 0.0, s1, s2, S, 0.0, J);
        job.ykw.model = &ykw_model;
        job.ykw.dyn = &ykw_dyn;
    }

    mu = input->m1 * input->m2 / (input->m1 + input->m2);
    rc = ere_swave(1000.0 * mu, scatter_job_v, &job, &ere);
    if (rc != 0) {
        fprintf(stderr, "Error: SCATTER failed to extract the phase shift (%d).\n", rc);
        exit(1);
    }
    *mu_gev = mu;
    *a0_fm = ere.a0_fm;
    *r0_fm = ere.r0_fm;
}

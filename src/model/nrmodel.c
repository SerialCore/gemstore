/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/model/nrmodel.h>

#include <math.h>

#define NR_UNPACK(ctx) \
    const argsNRModel_t *args_model = ((const nr_pot_ctx_t *)(ctx))->model; \
    const argsNRModelDy_t *args_dynmc = ((const nr_pot_ctx_t *)(ctx))->dyn; \
    (void)args_model; \
    (void)args_dynmc

double NRVt(double p, void *ctx)
{
    NR_UNPACK(ctx);
    double mi = args_dynmc->mi;
    double mj = args_dynmc->mj;
    double cent = args_dynmc->OCent;

    return cent * (mi + p * p / (2.0 * mi) + mj + p * p / (2.0 * mj));
}

double NRVt_quark(double p, void *ctx)
{
    /* Single quark m + p²/2m. Baryon spectator kinetic energy on λ. */
    NR_UNPACK(ctx);
    double mi = args_dynmc->mi;
    double cent = args_dynmc->OCent;

    return cent * (mi + p * p / (2.0 * mi));
}

double NRVcoul(double r, void *ctx)
{
    NR_UNPACK(ctx);
    if (r == 0.0) return 0.0;

    double Cij = args_dynmc->Cij;
    double cent = args_dynmc->OCent;
    double alpha_s = args_model->alpha_s;

    return Cij * cent * alpha_s / r;
}

double NRVconf(double r, void *ctx)
{
    NR_UNPACK(ctx);

    model_type_t model = args_dynmc->model;
    double b = args_model->b;
    double mu = args_model->mu;
    double c = args_model->c;
    double Cij = args_dynmc->Cij;
    double cent = args_dynmc->OCent;
    double confine;

    if (model == MODEL_NRSTRING) {
        confine = b * r + c;
    }
    else if (model == MODEL_NRSCREEN) {
        confine = b * (1.0 - exp(-mu * r)) / mu + c;
    }
    else {
        return 0.0;
    }

    return -0.75 * Cij * cent * confine;
}

double NRVcont(double r, void *ctx)
{
    NR_UNPACK(ctx);
    double mi = args_dynmc->mi;
    double mj = args_dynmc->mj;
    double Cij = args_dynmc->Cij;
    double sds = args_dynmc->OSdS;
    double alpha_s = args_model->alpha_s;
    double sigma = args_model->sigma;

    double gauss = sigma / sqrt(M_PI);
    double smeared = gauss * gauss * gauss * exp(-sigma * sigma * r * r);

    /* Meson Cij = -4/3 recovers 32 π αs / (9 m_i m_j). */
    return -(8.0 * M_PI / 3.0) * Cij * alpha_s * smeared * sds / (mi * mj);
}

double NRVsov(double r, void *ctx)
{
    NR_UNPACK(ctx);
    if (r == 0.0) return 0.0;

    system_type_t system = args_dynmc->system;
    double mi = args_dynmc->mi;
    double mj = args_dynmc->mj;
    double Cij = args_dynmc->Cij;
    double ldsi = args_dynmc->OLSi;
    double ldsj = args_dynmc->OLSj;
    double alpha_s = args_model->alpha_s;

    /* Meson exchange keeps quark j the same sign as quark i.
     * Baryon: −(r×p_j·s_j) and −(r×p_j·s_i − r×p_i·s_j). */
    double diag = ldsi / (mi * mi);
    double cross = 0.0;
    if (system == SYSTEM_MESON) {
        diag += ldsj / (mj * mj);
        cross = (ldsi + ldsj) / (mi * mj);
    }
    else if (system == SYSTEM_BARYON) {
        diag -= ldsj / (mj * mj);
        cross = -(ldsi - ldsj) / (mi * mj);
    }

    /* dG/dr, G = Cij αs / r. Meson Cij = -4/3 gives +4 αs / (3 r²). */
    double dG = -Cij * alpha_s / (r * r);

    return (diag + cross) * dG / r;
}

double NRVsos(double r, void *ctx)
{
    NR_UNPACK(ctx);
    if (r == 0.0) return 0.0;

    system_type_t system = args_dynmc->system;
    model_type_t model = args_dynmc->model;
    double mi = args_dynmc->mi;
    double mj = args_dynmc->mj;
    double Cij = args_dynmc->Cij;
    double ldsi = args_dynmc->OLSi;
    double ldsj = args_dynmc->OLSj;
    double b = args_model->b;
    double alpha_s = args_model->alpha_s;

    /* Overall minus times s_i/m_i². Meson adds s_j/m_j²; baryon subtracts it. */
    double coeff = ldsi / (mi * mi);
    if (system == SYSTEM_MESON) coeff += ldsj / (mj * mj);
    else if (system == SYSTEM_BARYON) coeff -= ldsj / (mj * mj);

    /* Coulomb plus string. -3/4 Cij makes a meson recover b, or b e^{-μr}. */
    double slope = -Cij * alpha_s / (r * r);
    if (model == MODEL_NRSTRING) {
        slope += -0.75 * Cij * b;
    }
    else if (model == MODEL_NRSCREEN) {
        slope += -0.75 * Cij * b * exp(-args_model->mu * r);
    }

    return -coeff / (2.0 * r) * slope;
}

double NRVtens(double r, void *ctx)
{
    NR_UNPACK(ctx);
    if (r == 0.0) return 0.0;

    double mi = args_dynmc->mi;
    double mj = args_dynmc->mj;
    double Cij = args_dynmc->Cij;
    double tens = args_dynmc->OTens;
    double alpha_s = args_model->alpha_s;

    /* Meson Cij = -4/3 recovers (4/3) αs / (m_i m_j r³) ⟨S_{12}⟩. */
    return -Cij * alpha_s * tens / (mi * mj * r * r * r);
}

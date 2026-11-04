/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Local Yukawa. Each term is (V/r) e^{-μ r} (W + B P_σ), P_σ=(1+σ1·σ2)/2.
 * V is MeV·fm, μ is fm^{-1}. Quadrature passes r in GeV^{-1}; kernels return GeV.
 */

#include <gemstore/model/ykwmodel.h>
#include <gemstore/basis/orbit.h>

#include <math.h>

#define YKW_UNPACK(ctx) \
    const argsYKWModel_t *args_model = ((const ykw_pot_ctx_t *)(ctx))->model; \
    const argsYKWModelDy_t *args_dynmc = ((const ykw_pot_ctx_t *)(ctx))->dyn; \
    (void)args_model; \
    (void)args_dynmc

double YKWVt(double p, void *ctx)
{
    YKW_UNPACK(ctx);
    double mi = args_dynmc->mi;
    double mj = args_dynmc->mj;
    double cent = args_dynmc->OCent;
    double mu;

    if (mi <= 0.0 || mj <= 0.0) {
        return 0.0;
    }
    mu = mi * mj / (mi + mj);

    return cent * (mi + mj + p * p / (2.0 * mu));
}

double YKWVc(double r, void *ctx)
{
    YKW_UNPACK(ctx);
    double cent = args_dynmc->OCent;
    double rfm;
    double vc = 0.0;
    int i;

    if (r < 1e-8) {
        r = 1e-8;
    }
    rfm = r / GEMSTORE_FM;
    for (i = 0; i < args_model->nterm; i++) {
        const argsYKWTerm_t *t = &args_model->term[i];
        double y;

        if (t->mu < 0.0) {
            continue;
        }
        /* MeV·fm / fm = MeV, then 10^{-3} → GeV. */
        y = 1.0e-3 * t->V * exp(-t->mu * rfm) / rfm;
        vc += (t->wigner + 0.5 * t->bartlett) * y;
    }

    return cent * vc;
}

double YKWVss(double r, void *ctx)
{
    YKW_UNPACK(ctx);
    double sds = args_dynmc->OSdS;
    double rfm;
    double vss = 0.0;
    int i;

    if (r < 1e-8) {
        r = 1e-8;
    }
    rfm = r / GEMSTORE_FM;
    for (i = 0; i < args_model->nterm; i++) {
        const argsYKWTerm_t *t = &args_model->term[i];
        double y;

        if (t->mu < 0.0) {
            continue;
        }
        y = 1.0e-3 * t->V * exp(-t->mu * rfm) / rfm;
        vss += 0.5 * t->bartlett * y;
    }

    return sds * vss;
}
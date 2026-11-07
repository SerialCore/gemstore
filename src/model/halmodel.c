/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Leading-order HAL QCD potential, Miyamoto thesis Eq. (3.11):
 *   V = V0(r) + Vσ(r) (σ1·σ2) + VT(r) S12.
 * Radial functions are tables. r below the first node keeps V(r_first).
 * r above the last node is zero. Values are MeV on input and GeV here.
 * L·S and coupled-channel matrices are not in this kernel.
 */

#include <gemstore/model/halmodel.h>
#include <gemstore/basis/orbit.h>

#include <math.h>

#define HAL_UNPACK(ctx) \
    const argsHALModel_t *args_model = ((const hal_pot_ctx_t *)(ctx))->model; \
    const argsHALModelDy_t *args_dynmc = ((const hal_pot_ctx_t *)(ctx))->dyn; \
    (void)args_model; \
    (void)args_dynmc

static double hal_interp(const argsHALPoint_t *pt, int n, double rfm)
{
    int lo, hi;

    if (n <= 0) {
        return 0.0;
    }
    if (rfm <= pt[0].r) {
        return pt[0].V;
    }
    if (rfm >= pt[n - 1].r) {
        return 0.0;
    }
    lo = 0;
    hi = n - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (pt[mid].r <= rfm) {
            lo = mid;
        }
        else {
            hi = mid;
        }
    }
    {
        double dr = pt[hi].r - pt[lo].r;
        double t = (dr > 0.0) ? (rfm - pt[lo].r) / dr : 0.0;
        return pt[lo].V + t * (pt[hi].V - pt[lo].V);
    }
}

double HALKin(double p, void *ctx)
{
    HAL_UNPACK(ctx);
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

double HALVc(double r, void *ctx)
{
    HAL_UNPACK(ctx);
    double rfm = r / GEMSTORE_FM;
    double v = hal_interp(args_model->V0, args_model->n0, rfm);

    return args_dynmc->OCent * (1.0e-3 * v);
}

double HALVss(double r, void *ctx)
{
    HAL_UNPACK(ctx);
    double rfm = r / GEMSTORE_FM;
    double v = hal_interp(args_model->Vsig, args_model->nsig, rfm);

    return args_dynmc->OSdS * (1.0e-3 * v);
}

double HALVten(double r, void *ctx)
{
    HAL_UNPACK(ctx);
    double rfm = r / GEMSTORE_FM;
    double v = hal_interp(args_model->VT, args_model->nt, rfm);

    return args_dynmc->OTens * (1.0e-3 * v);
}

double HALV1S0(const argsHALModel_t *model, double rfm)
{
    double v0 = hal_interp(model->V0, model->n0, rfm);
    double vs = hal_interp(model->Vsig, model->nsig, rfm);

    return 1.0e-3 * (v0 - 3.0 * vs);
}

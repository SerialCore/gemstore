/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Leading-order EFT contact of arXiv:2610.28021, Eq. (40):
 *   V(r) = -Ca exp(-(r/Rc)^2),
 * with Ca in MeV and Rc in fm. One-pion exchange is absent, and the
 * normalized Gaussian that replaces δ^(3)(r) is absorbed into Ca.
 * The potential is central: Vss = Vt = 0.
 */

#include <gemstore/model/eftmodel.h>
#include <gemstore/basis/orbit.h>

#include <math.h>

#define EFT_UNPACK(ctx) \
    const argsEFTModel_t *args_model = ((const eft_pot_ctx_t *)(ctx))->model; \
    const argsEFTModelDy_t *args_dynmc = ((const eft_pot_ctx_t *)(ctx))->dyn; \
    (void)args_model; \
    (void)args_dynmc

double EFTVt(double p, void *ctx)
{
    EFT_UNPACK(ctx);
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

double EFTVc(double r, void *ctx)
{
    EFT_UNPACK(ctx);
    double cent = args_dynmc->OCent;
    double rc = args_model->Rc;
    double rfm;
    double x;

    if (rc <= 0.0) {
        return 0.0;
    }
    rfm = r / GEMSTORE_FM;
    x = rfm / rc;
    /* MeV → GeV. V = -Ca exp(-(r/Rc)^2). */
    return cent * (-1.0e-3 * args_model->Ca * exp(-x * x));
}

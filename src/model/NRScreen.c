/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/model/NRScreen.h>
#include <gemstore/model/model.h>

#include <math.h>

const argsModel_t argsNRScreen_meson = {
    .mn = 0.606,
    .ms = 0.780,
    .mc = 1.984,
    .mb = 5.368,
    .alpha_s = 0.3930,
    .b1 = 0.2312,
    .mu = 0.069,
    .c = -1.1711,
    .sigma = 1.842
};

double NRScreen_T(double p, const argsModel_t *args_model, const argsModelDy_t *args_dynmc)
{
    double mi = args_dynmc->mi;
    double mj = args_dynmc->mj;
    double cent = args_dynmc->OCent;

    return cent * (mi + mj + p * p / (2.0 * mi) + p * p / (2.0 * mj));
}

double NRScreen_Vconf(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc)
{
    if (r == 0.0) {
        return 0.0;
    }
    
    double b1 = args_model->b1;
    double mu = args_model->mu;
    double c = args_model->c;
    double alpha_s = args_model->alpha_s;

    double Cij = args_dynmc->Cij;
    double cent = args_dynmc->OCent;

    double coul = alpha_s / r;
    double screen = b1 * (1.0 - exp(-mu * r)) / mu;
    
    return Cij * cent * coul - (3.0 / 4.0) * Cij * cent * screen - (3.0 / 4.0) * Cij * cent * c;
}

double NRScreen_Vcont(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc)
{
    double alpha_s = args_model->alpha_s;
    double sigma = args_model->sigma;

    double mi = args_dynmc->mi;
    double mj = args_dynmc->mj;
    double sds = args_dynmc->OSdS;

    double magnet = 32.0 * alpha_s * sds / (9.0 * mi * mj);
    double factor = pow(sigma, 3.0) / sqrt(M_PI);
    double exp_term = exp(-sigma * sigma * r * r);

    return magnet * factor * exp_term;
}

double NRScreen_Vsocm(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc)
{
    if (r == 0.0) {
        return 0.0;
    }

    double alpha_s = args_model->alpha_s;

    double mi = args_dynmc->mi;
    double mj = args_dynmc->mj;
    double Cij = args_dynmc->Cij;
    double ldsi = args_dynmc->OLSi;
    double ldsj = args_dynmc->OLSj;
    
    double tensor = alpha_s / pow(r, 3.0);
    double magnet = (ldsi / (mi * mi) + ldsj / (mj * mj) + (ldsi + ldsj) / (mi * mj));

    return -Cij * tensor * magnet;
}

double NRScreen_Vsotp(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc)
{
    if (r == 0.0) {
        return 0.0;
    }

    double alpha_s = args_model->alpha_s;
    double b1 = args_model->b1;
    double mu = args_model->mu;

    double mi = args_dynmc->mi;
    double mj = args_dynmc->mj;
    double Cij = args_dynmc->Cij;
    double ldsi = args_dynmc->OLSi;
    double ldsj = args_dynmc->OLSj;

    double dcoul_dr = -alpha_s / pow(r, 2.0);
    double dscreen_dr = b1 * exp(-mu * r);
    double dV_dr = Cij * dcoul_dr - (3.0 / 4.0) * Cij * dscreen_dr;
    double magnet = (ldsi / (mi * mi) + ldsj / (mj * mj));
    
    return -0.5 / r * dV_dr * magnet;
}

double NRScreen_Vtens(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc)
{
    if (r == 0.0) {
        return 0.0;
    }

    double alpha_s = args_model->alpha_s;

    double mi = args_dynmc->mi;
    double mj = args_dynmc->mj;
    double tens = args_dynmc->OTens;

    double tensor = alpha_s / pow(r, 3.0);
    double magnet = 4.0 * tens / (3.0 * mi * mj);

    return tensor * magnet;
}
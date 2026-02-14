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
    double m1 = args_dynmc->m1;
    double m2 = args_dynmc->m2;

    return m1 + m2 + p * p / (2.0 * m1) + p * p / (2.0 * m2);
}

double NRScreen_Vconf(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc)
{
    if (r == 0.0) {
        return 0.0;
    }
    
    double b1 = args_model->b1;
    double mu = args_model->mu;
    double c = args_model->c;

    double C12 = args_dynmc->C12;
    double alpha_s = args_model->alpha_s;

    double coul = alpha_s / r;
    double screen = b1 * (1.0 - exp(-mu * r)) / mu;
    
    return C12 * coul - (3.0 / 4.0) * C12 * screen - (3.0 / 4.0) * C12 * c;
}

double NRScreen_Vcont(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc)
{
    double alpha_s = args_model->alpha_s;
    double sigma = args_model->sigma;

    double m1 = args_dynmc->m1;
    double m2 = args_dynmc->m2;
    double sds = args_dynmc->OSdS;

    double magnet = 32.0 * alpha_s * sds / (9.0 * m1 * m2);
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

    double m1 = args_dynmc->m1;
    double m2 = args_dynmc->m2;
    double C12 = args_dynmc->C12;
    double lds1 = args_dynmc->OLS1;
    double lds2 = args_dynmc->OLS2;
    
    double tensor = alpha_s / pow(r, 3.0);
    double magnet = (lds1 / (m1 * m1) + lds2 / (m2 * m2) + (lds1 + lds2) / (m1 * m2));

    return -C12 * tensor * magnet;
}

double NRScreen_Vsotp(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc)
{
    if (r == 0.0) {
        return 0.0;
    }

    double alpha_s = args_model->alpha_s;
    double b1 = args_model->b1;
    double mu = args_model->mu;

    double m1 = args_dynmc->m1;
    double m2 = args_dynmc->m2;
    double C12 = args_dynmc->C12;
    double lds1 = args_dynmc->OLS1;
    double lds2 = args_dynmc->OLS2;

    double dcoul_dr = -alpha_s / pow(r, 2.0);
    double dscreen_dr = b1 * exp(-mu * r);
    double dV_dr = C12 * dcoul_dr - (3.0 / 4.0) * C12 * dscreen_dr;
    double magnet = (lds1 / (m1 * m1) + lds2 / (m2 * m2));
    
    return -0.5 / r * dV_dr * magnet;
}

double NRScreen_Vtens(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc)
{
    if (r == 0.0) {
        return 0.0;
    }

    double alpha_s = args_model->alpha_s;

    double m1 = args_dynmc->m1;
    double m2 = args_dynmc->m2;
    double tens = args_dynmc->OTens;

    double tensor = alpha_s / pow(r, 3.0);
    double magnet = 4.0 * tens / (3.0 * m1 * m2);

    return tensor * magnet;
}
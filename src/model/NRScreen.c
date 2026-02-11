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
    .mt = 172.57,
    .alpha_s = 0.3930,
    .b1 = 0.2312,
    .mu = 0.069,
    .c = -1.1711,
    .sigma = 1.842
};

double NRScreen_T(double p, const argsFlavor_t *args_flavor, const argsSOC_t *args_soc, const argsModel_t *args_model)
{
    double m1 = args_flavor->m1;
    double m2 = args_flavor->m2;

    return m1 + m2 + p * p / (2.0 * m1) + p * p / (2.0 * m2);
}

double NRScreen_Vconf(double r, const argsFlavor_t *args_flavor, const argsSOC_t *args_soc, const argsModel_t *args_model)
{
    if (r == 0.0) {
        return 0.0;
    }
    
    double Cij = -4.0 / 3.0;
    double alpha_s = args_model->alpha_s;
    double b1 = args_model->b1;
    double mu = args_model->mu;
    double c = args_model->c;

    double coul = alpha_s / r;
    double screen = b1 * (1.0 - exp(-mu * r)) / mu;
    
    return Cij * coul - (3.0 / 4.0) * Cij * screen - (3.0 / 4.0) * Cij * c;
}

double NRScreen_Vcont(double r, const argsFlavor_t *args_flavor, const argsSOC_t *args_soc, const argsModel_t *args_model)
{
    double m1 = args_flavor->m1;
    double m2 = args_flavor->m2;
    double sds = args_soc->OSdS;
    double alpha_s = args_model->alpha_s;
    double sigma = args_model->sigma;

    double magnet = 32.0 * alpha_s * sds / (9.0 * m1 * m2);
    double factor = pow(sigma, 3.0) / sqrt(M_PI);
    double exp_term = exp(-sigma * sigma * r * r);

    return magnet * factor * exp_term;
}

double NRScreen_Vsocm(double r, const argsFlavor_t *args_flavor, const argsSOC_t *args_soc, const argsModel_t *args_model)
{
    if (r == 0.0) {
        return 0.0;
    }

    double m1 = args_flavor->m1;
    double m2 = args_flavor->m2;
    double lds1 = args_soc->OLS1;
    double lds2 = args_soc->OLS2;
    double Cij = -4.0 / 3.0;
    double alpha_s = args_model->alpha_s;
    
    double tensor = alpha_s / pow(r, 3.0);
    double magnet = (lds1 / (m1 * m1) + lds2 / (m2 * m2) + (lds1 + lds2) / (m1 * m2));

    return -Cij * tensor * magnet;
}

double NRScreen_Vsotp(double r, const argsFlavor_t *args_flavor, const argsSOC_t *args_soc, const argsModel_t *args_model)
{
    if (r == 0.0) {
        return 0.0;
    }

    double m1 = args_flavor->m1;
    double m2 = args_flavor->m2;
    double lds1 = args_soc->OLS1;
    double lds2 = args_soc->OLS2;
    double Cij = -4.0 / 3.0;
    double alpha_s = args_model->alpha_s;
    double b1 = args_model->b1;
    double mu = args_model->mu;

    double dcoul_dr = -alpha_s / pow(r, 2.0);
    double dscreen_dr = b1 * exp(-mu * r);
    double dV_dr = Cij * dcoul_dr - (3.0 / 4.0) * Cij * dscreen_dr;
    double magnet = (lds1 / (m1 * m1) + lds2 / (m2 * m2));
    
    return -0.5 / r * dV_dr * magnet;
}

double NRScreen_Vtens(double r, const argsFlavor_t *args_flavor, const argsSOC_t *args_soc, const argsModel_t *args_model)
{
    if (r == 0.0) {
        return 0.0;
    }

    double m1 = args_flavor->m1;
    double m2 = args_flavor->m2;
    double tens = args_soc->OTens;
    double alpha_s = args_model->alpha_s;
    
    double tensor = alpha_s / pow(r, 3.0);
    double magnet = 4.0 * tens / (3.0 * m1 * m2);

    return tensor * magnet;
}
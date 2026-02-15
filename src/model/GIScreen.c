/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/model/GIScreen.h>
#include <gemstore/model/model.h>

#include <math.h>

const argsModel_t argsGIScreen_meson = {
    /* quark masses (GeV) */
    .mn = 0.220,
    .ms = 0.419,
    .mc = 1.628,
    .mb = 4.977,
    .mt = 172.57,

    /* GI screened string model parameters */
    .b1 = 0.18,
    .mu = 0.15,
    .c = -0.253,
    .sigma_0 = 1.8,
    .s = 1.55,
    
    /* GI smearing / contact / spin-orbit / tensor coefficients */
    .epsilon_Coul = 0.0,
    .epsilon_cont = -0.168,
    .epsilon_sonu = -0.035,
    .epsilon_sos = 0.055,
    .epsilon_tens = 0.025,
};

const double GI_ALPHA_K[3] = {0.25, 0.15, 0.20};
const double GI_GAMMA_K[3] = {0.5, 1.5811388300841898, 15.811388300841896};

double GIScreen_T(double p, const argsModel_t *args_model, const argsModelDy_t *args_dynmc)
{
    double mi = args_dynmc->mi;
    double mj = args_dynmc->mj;
    double cent = args_dynmc->OCent;

    return cent * sqrt(mi * mi + p * p) + cent * sqrt(mj * mj + p * p);
}

double GIScreen_betaij(double p, const argsModel_t *args_model, const argsModelDy_t *args_dynmc)
{
    double mi = args_dynmc->mi;
    double mj = args_dynmc->mj;
    double cent = args_dynmc->OCent;

    return cent * (1.0 + p * p / (sqrt(p * p + mi * mi) * sqrt(p * p + mj * mj)));
}

double GIScreen_deltaij(double p, const argsModel_t *args_model, const argsModelDy_t *args_dynmc)
{
    double mi = args_dynmc->mi;
    double mj = args_dynmc->mj;
    double cent = args_dynmc->OCent;

    return cent * mi * mj / (sqrt(p * p + mi * mi) * sqrt(p * p + mj * mj));
}

double GIScreen_deltaii(double p, const argsModel_t *args_model, const argsModelDy_t *args_dynmc)
{
    double mi = args_dynmc->mi;
    double cent = args_dynmc->OCent;

    return cent * mi * mi / (sqrt(p * p + mi * mi) * sqrt(p * p + mi * mi));
}

double GIScreen_deltajj(double p, const argsModel_t *args_model, const argsModelDy_t *args_dynmc)
{
    double mj = args_dynmc->mj;
    double cent = args_dynmc->OCent;

    return cent * mj * mj / (sqrt(p * p + mj * mj) * sqrt(p * p + mj * mj));
}

double GIScreen_Vcoul(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc)
{
    if (r == 0.0) {
        return 0.0;
    }

    double Cij = args_dynmc->Cij;
    double cent = args_dynmc->OCent;
    double sigmak[3] = {args_dynmc->Sigkij[0], args_dynmc->Sigkij[1], args_dynmc->Sigkij[2]};

    double sum = 0.0;
    for (int k = 0; k < 3; k++) {
        sum += GI_ALPHA_K[k] * erf(sigmak[k] * r);
    }

    return Cij * cent * sum / r;
}

double GIScreen_Vconf(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc)
{
    if (r == 0.0) {
        return 0.0;
    }

    double b1 = args_model->b1;
    double mu = args_model->mu;
    double c = args_model->c;
    double Cij = args_dynmc->Cij;
    double cent = args_dynmc->OCent;
    double sigmaij = args_dynmc->Sigij;

    double sig2 = sigmaij * sigmaij;
    double mu2_4sig2 = mu * mu / (4.0 * sig2);
    double mu_m_2rsig2 = mu - 2.0 * r * sig2;
    double mu_p_2rsig2 = mu + 2.0 * r * sig2;

    double pref = -3.0 * Cij * cent * b1 * exp(-r * mu) / (16.0 * r * mu * sig2);
    double inner1 = 4.0 * r * sig2 * exp(r * mu);
    double inner2 = mu_m_2rsig2 * exp(mu2_4sig2) * erfc(mu_m_2rsig2 / (2.0 * sigmaij));
    double inner3 = mu_p_2rsig2 * exp(mu2_4sig2 + 2 * r * mu) * erfc(mu_p_2rsig2 / (2.0 * sigmaij));

    return pref * (inner1 + inner2 - inner3) - 0.75 * Cij * cent * c;
}

double GIScreen_Vcont(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc)
{
    double mi = args_dynmc->mi;
    double mj = args_dynmc->mj;
    double Cij = args_dynmc->Cij;
    double sds = args_dynmc->OSdS;
    double sigmak[3] = {args_dynmc->Sigkij[0], args_dynmc->Sigkij[1], args_dynmc->Sigkij[2]};

    double pref = -4.0 * Cij * M_2_SQRTPI * sds / (3.0 * mi * mj);

    double sum = 0.0;
    for (int k = 0; k < 3; k++) {
        double sk = sigmak[k];
        sum += GI_ALPHA_K[k] * (sk * sk * sk) * exp(-sk * sk * r * r);
    }

    return pref * sum;
}

double GIScreen_Vsovi(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc)
{
    if (r == 0.0) {
        return 0.0;
    }

    double mi = args_dynmc->mi;
    double Cij = args_dynmc->Cij;
    double ldsi = args_dynmc->OLSi;
    double sigmak[3] = {args_dynmc->Sigkij[0], args_dynmc->Sigkij[1], args_dynmc->Sigkij[2]};

    double pref = -ldsi / (2 * r * mi * mi);

    double dVcoul_dr = 0.0;
    for (int k = 0; k < 3; k++) {
        dVcoul_dr += Cij * M_2_SQRTPI * exp(-r * r * sigmak[k] * sigmak[k]) * GI_ALPHA_K[k] * sigmak[k] / r
            - Cij * GI_ALPHA_K[k] * erf(r * sigmak[k]) / (r * r);
    }

    return pref * dVcoul_dr;
}

double GIScreen_Vsovj(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc)
{
    if (r == 0.0) {
        return 0.0;
    }

    double mj = args_dynmc->mj;
    double Cij = args_dynmc->Cij;
    double ldsj = args_dynmc->OLSj;
    double sigmak[3] = {args_dynmc->Sigkij[0], args_dynmc->Sigkij[1], args_dynmc->Sigkij[2]};

    double pref = ldsj / (2 * r * mj * mj);

    double dVcoul_dr = 0.0;
    for (int k = 0; k < 3; k++) {
        dVcoul_dr += Cij * M_2_SQRTPI * exp(-r * r * sigmak[k] * sigmak[k]) * GI_ALPHA_K[k] * sigmak[k] / r
            - Cij * GI_ALPHA_K[k] * erf(r * sigmak[k]) / (r * r);
    }

    return pref * dVcoul_dr;
}

double GIScreen_Vsovij(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc)
{
    if (r == 0.0) {
        return 0.0;
    }

    double mi = args_dynmc->mi;
    double mj = args_dynmc->mj;
    double Cij = args_dynmc->Cij;
    double ldsi = args_dynmc->OLSi;
    double ldsj = args_dynmc->OLSj;
    double sigmak[3] = {args_dynmc->Sigkij[0], args_dynmc->Sigkij[1], args_dynmc->Sigkij[2]};

    double pref = (ldsi + ldsj) / (2 * r * mi * mj);

    double dVcoul_dr = 0.0;
    for (int k = 0; k < 3; k++) {
        dVcoul_dr += Cij * M_2_SQRTPI * exp(-r * r * sigmak[k] * sigmak[k]) * GI_ALPHA_K[k] * sigmak[k] / r
            - Cij * GI_ALPHA_K[k] * erf(r * sigmak[k]) / (r * r);
    }

    return pref * dVcoul_dr;
}

double GIScreen_Vsosi(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc)
{
    if (r == 0.0) {
        return 0.0;
    }

    double mi = args_dynmc->mi;
    double Cij = args_dynmc->Cij;
    double ldsi = args_dynmc->OLSi;
    double sigmak[3] = {args_dynmc->Sigkij[0], args_dynmc->Sigkij[1], args_dynmc->Sigkij[2]};

    double pref = -ldsi / (2 * r * mi * mi);

    /* numerical differential */
    double h = 1e-6 * (r > 1e-8 ? r : 1.0);
    double dVconf_dr = (GIScreen_Vconf(r + h, args_model, args_dynmc) - GIScreen_Vconf(r - h, args_model, args_dynmc)) / (2.0 * h);

    return pref * dVconf_dr;
}

double GIScreen_Vsosj(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc)
{
    if (r == 0.0) {
        return 0.0;
    }

    double mj = args_dynmc->mj;
    double Cij = args_dynmc->Cij;
    double ldsj = args_dynmc->OLSj;
    double sigmak[3] = {args_dynmc->Sigkij[0], args_dynmc->Sigkij[1], args_dynmc->Sigkij[2]};

    double pref = -ldsj / (2 * r * mj * mj);

    /* numerical differential */
    double h = 1e-6 * (r > 1e-8 ? r : 1.0);
    double dVconf_dr = (GIScreen_Vconf(r + h, args_model, args_dynmc) - GIScreen_Vconf(r - h, args_model, args_dynmc)) / (2.0 * h);

    return pref * dVconf_dr;
}

double GIScreen_Vtens(double r, const argsModel_t *args_model, const argsModelDy_t *args_dynmc)
{
    if (r == 0.0) {
        return 0.0;
    }

    double mi = args_dynmc->mi;
    double mj = args_dynmc->mj;
    double Cij = args_dynmc->Cij;
    double tens = args_dynmc->OTens;
    double sigmak[3] = {args_dynmc->Sigkij[0], args_dynmc->Sigkij[1], args_dynmc->Sigkij[2]};

    double pref = tens / (3 * mi * mj);

    double sum = 0.0;
    double dVcoul_dr;
    double d2Vcoul_dr2;
    double exp_r2sigk2;
    double erf_rsigk;
    for (int k = 0; k < 3; k++) {
        exp_r2sigk2 = exp(-r * r * sigmak[k] * sigmak[k]);
        erf_rsigk = erf(r * sigmak[k]);

        dVcoul_dr = Cij * M_2_SQRTPI * exp_r2sigk2 * GI_ALPHA_K[k] * sigmak[k] / r
            - Cij * GI_ALPHA_K[k] * erf_rsigk / (r * r);

        d2Vcoul_dr2 = -2 * Cij * M_2_SQRTPI * exp_r2sigk2 * GI_ALPHA_K[k] * sigmak[k] / (r * r)
            - 2 * Cij * M_2_SQRTPI * exp_r2sigk2 * GI_ALPHA_K[k] * sigmak[k] * sigmak[k] * sigmak[k]
            + 2 * Cij * GI_ALPHA_K[k] * erf_rsigk / (r * r * r);

        sum += dVcoul_dr / r - d2Vcoul_dr2;
    }

    return pref * sum;
}
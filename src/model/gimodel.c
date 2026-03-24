/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/model/gimodel.h>
#include <gemstore/param/argset.h>

#include <math.h>

const double GI_ALPHA_K[3] = {0.25, 0.15, 0.20};
const double GI_GAMMA_K[3] = {0.5, 1.5811388300841898, 15.811388300841896};

double GIVt(double p, const argsGIModel_t *args_model, const argsGIModelDy_t *args_dynmc)
{
    double mi = args_dynmc->mi;
    double mj = args_dynmc->mj;
    double cent = args_dynmc->OCent;

    return cent * sqrt(mi * mi + p * p) + cent * sqrt(mj * mj + p * p);
}

double GIVbetaijcoul(double p, const argsGIModel_t *args_model, const argsGIModelDy_t *args_dynmc)
{
    double mi = args_dynmc->mi;
    double mj = args_dynmc->mj;
    double cent = args_dynmc->OCent;

    double betaij = cent * (1.0 + p * p / (sqrt(p * p + mi * mi) * sqrt(p * p + mj * mj)));

    return pow(betaij, 0.5);
}

double GIVdeltaijcont(double p, const argsGIModel_t *args_model, const argsGIModelDy_t *args_dynmc)
{
    double epsilon_cont = args_model->epsilon_cont;
    double mi = args_dynmc->mi;
    double mj = args_dynmc->mj;
    double cent = args_dynmc->OCent;

    double deltaij = cent * mi * mj / (sqrt(p * p + mi * mi) * sqrt(p * p + mj * mj));

    return pow(deltaij, 0.5 + epsilon_cont + epsilon_cont * epsilon_cont);
}

double GIVdeltaiisov(double p, const argsGIModel_t *args_model, const argsGIModelDy_t *args_dynmc)
{
    double epsilon_sov = args_model->epsilon_sov;
    double mi = args_dynmc->mi;
    double cent = args_dynmc->OCent;

    double deltaii = cent * mi * mi / (sqrt(p * p + mi * mi) * sqrt(p * p + mi * mi));

    return pow(deltaii, 0.5 + epsilon_sov + epsilon_sov * epsilon_sov);
}

double GIVdeltajjsov(double p, const argsGIModel_t *args_model, const argsGIModelDy_t *args_dynmc)
{
    double epsilon_sov = args_model->epsilon_sov;
    double mj = args_dynmc->mj;
    double cent = args_dynmc->OCent;

    double deltajj = cent * mj * mj / (sqrt(p * p + mj * mj) * sqrt(p * p + mj * mj));

    return pow(deltajj, 0.5 + epsilon_sov + epsilon_sov * epsilon_sov);
}

double GIVdeltaijsov(double p, const argsGIModel_t *args_model, const argsGIModelDy_t *args_dynmc)
{
    double epsilon_sov = args_model->epsilon_sov;
    double mi = args_dynmc->mi;
    double mj = args_dynmc->mj;
    double cent = args_dynmc->OCent;

    double deltaij = cent * mi * mj / (sqrt(p * p + mi * mi) * sqrt(p * p + mj * mj));

    return pow(deltaij, 0.5 + epsilon_sov + epsilon_sov * epsilon_sov);
}

double GIVdeltaiisos(double p, const argsGIModel_t *args_model, const argsGIModelDy_t *args_dynmc)
{
    double epsilon_sos = args_model->epsilon_sos;
    double mi = args_dynmc->mi;
    double cent = args_dynmc->OCent;

    double deltaii = cent * mi * mi / (sqrt(p * p + mi * mi) * sqrt(p * p + mi * mi));

    return pow(deltaii, 0.5 + epsilon_sos + epsilon_sos * epsilon_sos);
}

double GIVdeltajjsos(double p, const argsGIModel_t *args_model, const argsGIModelDy_t *args_dynmc)
{
    double epsilon_sos = args_model->epsilon_sos;
    double mj = args_dynmc->mj;
    double cent = args_dynmc->OCent;

    double deltajj = cent * mj * mj / (sqrt(p * p + mj * mj) * sqrt(p * p + mj * mj));

    return pow(deltajj, 0.5 + epsilon_sos + epsilon_sos * epsilon_sos);
}

double GIVdeltaijtens(double p, const argsGIModel_t *args_model, const argsGIModelDy_t *args_dynmc)
{
    double epsilon_tens = args_model->epsilon_tens;
    double mi = args_dynmc->mi;
    double mj = args_dynmc->mj;
    double cent = args_dynmc->OCent;

    double deltaij = cent * mi * mj / (sqrt(p * p + mi * mi) * sqrt(p * p + mj * mj));

    return pow(deltaij, 0.5 + epsilon_tens + epsilon_tens * epsilon_tens);
}

double GIVcoul(double r, const argsGIModel_t *args_model, const argsGIModelDy_t *args_dynmc)
{
    if (r == 0.0) return 0.0;

    double Cij = args_dynmc->Cij;
    double cent = args_dynmc->OCent;
    double sigmak[3] = {args_dynmc->Sigkij[0], args_dynmc->Sigkij[1], args_dynmc->Sigkij[2]};

    double sum = 0.0;
    for (int k = 0; k < 3; k++) {
        sum += GI_ALPHA_K[k] * erf(sigmak[k] * r);
    }

    return Cij * cent * sum / r;
}

double GIVconf(double r, const argsGIModel_t *args_model, const argsGIModelDy_t *args_dynmc)
{
    if (r == 0.0) return 0.0;

    model_type_t model = args_dynmc->model;
    double b1 = args_model->b1;
    double b2 = args_model->b2;
    double mu = args_model->mu;
    double c = args_model->c;
    double Cij = args_dynmc->Cij;
    double cent = args_dynmc->OCent;
    double sigmaij = args_dynmc->Sigij;

    if (model == MODEL_GI_STRING) {
        double rsig = r * sigmaij;

        double pref = -3.0 * Cij * cent * b1 / (8.0 * r * sigmaij * sigmaij);
        double inner1 = M_2_SQRTPI * r * sigmaij * exp(-rsig * rsig);
        double inner2 = (1.0 + 2.0 * rsig * rsig) * erf(rsig);
        
        return pref * (inner1 + inner2) - 0.75 * Cij * cent * c;
    }
    else if (model == MODEL_GI_SCREEN || model == MODEL_GI_QUADRA) {
        double sig2 = sigmaij * sigmaij;
        double mu2_4sig2 = mu * mu / (4.0 * sig2);
        double mu_m_2rsig2 = mu - 2.0 * r * sig2;
        double mu_p_2rsig2 = mu + 2.0 * r * sig2;

        double pref = -3.0 * Cij * cent * b1 * exp(-r * mu) / (16.0 * r * mu * sig2);
        double inner1 = 4.0 * r * sig2 * exp(r * mu);
        double inner2 = mu_m_2rsig2 * exp(mu2_4sig2) * erfc(mu_m_2rsig2 / (2.0 * sigmaij));
        double inner3 = mu_p_2rsig2 * exp(mu2_4sig2 + 2 * r * mu) * erfc(mu_p_2rsig2 / (2.0 * sigmaij));

        if (model == MODEL_GI_SCREEN) {
            return pref * (inner1 + inner2 - inner3) - 0.75 * Cij * cent * c;
        }
        if (model == MODEL_GI_QUADRA) {
            double pref_quadra = -3.0 * Cij * cent * b2 / (4.0 * mu);
            double inner_quadra = 1.0 - sig2 * sigmaij * exp(-mu * r * r * sig2 / (mu + sig2)) / pow(mu + sig2, 1.5);
            
            return pref * (inner1 + inner2 - inner3) + pref_quadra * inner_quadra - 0.75 * Cij * cent * c;
        }
    }

    return 0.0;
}

double GIVcont(double r, const argsGIModel_t *args_model, const argsGIModelDy_t *args_dynmc)
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

double GIVsovi(double r, const argsGIModel_t *args_model, const argsGIModelDy_t *args_dynmc)
{
    if (r == 0.0) return 0.0;

    double mi = args_dynmc->mi;
    double Cij = args_dynmc->Cij;
    double ldsi = args_dynmc->OLSi;
    double sigmak[3] = {args_dynmc->Sigkij[0], args_dynmc->Sigkij[1], args_dynmc->Sigkij[2]};

    double pref = ldsi / (2 * r * mi * mi);

    double dVcoul_dr = 0.0;
    for (int k = 0; k < 3; k++) {
        dVcoul_dr += Cij * M_2_SQRTPI * exp(-r * r * sigmak[k] * sigmak[k]) * GI_ALPHA_K[k] * sigmak[k] / r
            - Cij * GI_ALPHA_K[k] * erf(r * sigmak[k]) / (r * r);
    }

    return pref * dVcoul_dr;
}

double GIVsovj(double r, const argsGIModel_t *args_model, const argsGIModelDy_t *args_dynmc)
{
    if (r == 0.0) return 0.0;

    system_type_t system = args_dynmc->system;
    double mj = args_dynmc->mj;
    double Cij = args_dynmc->Cij;
    double ldsj = args_dynmc->OLSj;
    double sigmak[3] = {args_dynmc->Sigkij[0], args_dynmc->Sigkij[1], args_dynmc->Sigkij[2]};

    double sign = 0.0;
    if (system == SYSTEM_MESON) sign = 1.0;
    else if (system == SYSTEM_BARYON) sign = -1.0;
    double pref = sign * ldsj / (2 * r * mj * mj);

    double dVcoul_dr = 0.0;
    for (int k = 0; k < 3; k++) {
        dVcoul_dr += Cij * M_2_SQRTPI * exp(-r * r * sigmak[k] * sigmak[k]) * GI_ALPHA_K[k] * sigmak[k] / r
            - Cij * GI_ALPHA_K[k] * erf(r * sigmak[k]) / (r * r);
    }

    return pref * dVcoul_dr;
}

double GIVsovij(double r, const argsGIModel_t *args_model, const argsGIModelDy_t *args_dynmc)
{
    if (r == 0.0) return 0.0;

    system_type_t system = args_dynmc->system;
    double mi = args_dynmc->mi;
    double mj = args_dynmc->mj;
    double Cij = args_dynmc->Cij;
    double ldsi = args_dynmc->OLSi;
    double ldsj = args_dynmc->OLSj;
    double sigmak[3] = {args_dynmc->Sigkij[0], args_dynmc->Sigkij[1], args_dynmc->Sigkij[2]};

    double sign = 0.0;
    if (system == SYSTEM_MESON) sign = 1.0;
    else if (system == SYSTEM_BARYON) sign = -1.0;
    double pref = sign * (ldsi + sign * ldsj) / (r * mi * mj);

    double dVcoul_dr = 0.0;
    for (int k = 0; k < 3; k++) {
        dVcoul_dr += Cij * M_2_SQRTPI * exp(-r * r * sigmak[k] * sigmak[k]) * GI_ALPHA_K[k] * sigmak[k] / r
            - Cij * GI_ALPHA_K[k] * erf(r * sigmak[k]) / (r * r);
    }

    return pref * dVcoul_dr;
}

double GIVsosi(double r, const argsGIModel_t *args_model, const argsGIModelDy_t *args_dynmc)
{
    if (r == 0.0) return 0.0;

    double mi = args_dynmc->mi;
    double ldsi = args_dynmc->OLSi;

    double pref = -ldsi / (2 * r * mi * mi);

    /* numerical differential */
    double h = 1e-6 * (r > 1e-8 ? r : 1.0);
    double dVconf_dr = (GIVconf(r + h, args_model, args_dynmc) - GIVconf(r - h, args_model, args_dynmc)) / (2.0 * h);

    return pref * dVconf_dr;
}

double GIVsosj(double r, const argsGIModel_t *args_model, const argsGIModelDy_t *args_dynmc)
{
    if (r == 0.0) return 0.0;

    system_type_t system = args_dynmc->system;
    double mj = args_dynmc->mj;
    double ldsj = args_dynmc->OLSj;

    double sign = 0.0;
    if (system == SYSTEM_MESON) sign = -1.0;
    else if (system == SYSTEM_BARYON) sign = 1.0;
    double pref = sign * ldsj / (2 * r * mj * mj);

    /* numerical differential */
    double h = 1e-6 * (r > 1e-8 ? r : 1.0);
    double dVconf_dr = (GIVconf(r + h, args_model, args_dynmc) - GIVconf(r - h, args_model, args_dynmc)) / (2.0 * h);

    return pref * dVconf_dr;
}

double GIVtens(double r, const argsGIModel_t *args_model, const argsGIModelDy_t *args_dynmc)
{
    if (r == 0.0) return 0.0;

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
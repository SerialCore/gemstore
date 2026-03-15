/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/param/fitccbar.h>
#include <gemstore/param/fitting.h>
#include <gemstore/param/helper.h>

#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnPrint.h>

#include <vector>
#include <ctime>

static const std::vector<State> data = {
    // charmonium (c c-bar)
    {3, 3, 1, 0, 0, 0, 2984.1,   5},   // ηc(1S)
    {3, 3, 2, 0, 0, 0, 3637.8,   5},   // ηc(2S)
    {3, 3, 1, 1, 0, 1, 3096.9,   5},   // J/ψ(1S)
    {3, 3, 2, 1, 0, 1, 3686.1,   5},   // ψ(2S)
    {3, 3, 1, 0, 1, 1, 3525.4,   5},   // hc(1P)
    {3, 3, 1, 1, 1, 0, 3414.7,   5},   // χc0(1P)
    {3, 3, 1, 1, 1, 1, 3510.7,   5},   // χc1(1P)
    {3, 3, 1, 1, 1, 2, 3556.2,   5},   // χc2(1P)
};

void minuit2_ccbar_GIScreen(double *params_out)
{
    srand(time(0));
    DualStream dual("Fitting.out");

    /* set parameters */
    ROOT::Minuit2::MnUserParameters upar;
    //upar.Add(name, value, init_step, lower_limit, upper_limit);
    upar.Add("mn", 0.220, 0.01, 0.1, 0.5);
    upar.Add("ms", 0.419, 0.01, 0.3, 0.7);
    upar.Add("mc", 1.628, 0.01, 1.5, 2.0);
    upar.Add("mb", 4.977, 0.01, 4.5, 5.5);
    upar.Add("b1", 0.18, 0.01, 0.1, 0.3);
    upar.Add("mu", 0.15, 0.01, 0.1, 0.2);
    upar.Add("c", -0.253, 0.01, -2.0, 0.0);
    upar.Add("sig0", 1.8, 0.01, 1.0, 3.0);
    upar.Add("s", 1.55, 0.01, 1.0, 3.0);
    upar.Add("econt", -0.168, 0.01, -0.5, 0.0);
    upar.Add("esov", -0.035, 0.01, -1.0, 1.0);
    upar.Add("esos", 0.055, 0.01, -1.0, 1.0);
    upar.Add("etens", 0.025, 0.01, -1.0, 1.0);
    int N_PARAMS = upar.Params().size();
    int N_DATA = data.size();

    /* use of Migrad algorithm with strategy 2, high precision */
    Chi2Minimizer minuit_fit(data, MODEL_GI_SCREEN, N_PARAMS, N_DATA, 1.0);
    ROOT::Minuit2::MnMigrad migrad(minuit_fit, upar, 2);

    /* perform the fit */
    ROOT::Minuit2::FunctionMinimum min_result = migrad();
    compute_chi2(data, min_result.UserParameters().Params(), MODEL_GI_SCREEN, true);
    dual << min_result.UserParameters() << std::endl;

    auto params = min_result.UserParameters().Params();
    for (int i = 0; i < N_PARAMS; i++) {
        params_out[i] = params[i];
    }
}

void minuit2_ccbar_GIQuadra(double *params_out)
{
    srand(time(0));
    DualStream dual("Fitting.out");

    /* set parameters */
    ROOT::Minuit2::MnUserParameters upar;
    //upar.Add(name, value, init_step, lower_limit, upper_limit);
    upar.Add("mn", 0.220, 0.01, 0.1, 0.5);
    upar.Add("ms", 0.419, 0.01, 0.3, 0.7);
    upar.Add("mc", 1.628, 0.01, 1.5, 2.0);
    upar.Add("mb", 4.977, 0.01, 4.5, 5.5);
    upar.Add("b1", 0.18, 0.01, 0.1, 0.3);
    upar.Add("b2", 0.02, 0.01, 0.0, 0.1);
    upar.Add("mu", 0.15, 0.01, 0.1, 0.2);
    upar.Add("c", -0.253, 0.01, -2.0, 0.0);
    upar.Add("sig0", 1.8, 0.01, 1.0, 3.0);
    upar.Add("s", 1.55, 0.01, 1.0, 3.0);
    upar.Add("econt", -0.168, 0.01, -0.5, 0.0);
    upar.Add("esov", -0.035, 0.01, -1.0, 1.0);
    upar.Add("esos", 0.055, 0.01, -1.0, 1.0);
    upar.Add("etens", 0.025, 0.01, -1.0, 1.0);
    int N_PARAMS = upar.Params().size();
    int N_DATA = data.size();

    /* use of Migrad algorithm with strategy 2, high precision */
    Chi2Minimizer minuit_fit(data, MODEL_GI_QUADRA, N_PARAMS, N_DATA, 1.0);
    ROOT::Minuit2::MnMigrad migrad(minuit_fit, upar, 2);

    /* perform the fit */
    ROOT::Minuit2::FunctionMinimum min_result = migrad();
    compute_chi2(data, min_result.UserParameters().Params(), MODEL_GI_QUADRA, true);
    dual << min_result.UserParameters() << std::endl;

    auto params = min_result.UserParameters().Params();
    for (int i = 0; i < N_PARAMS; i++) {
        params_out[i] = params[i];
    }
}
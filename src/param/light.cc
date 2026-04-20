/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/param/light.h>
#include <gemstore/param/minuit.h>
#include <gemstore/types.h>

#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnPrint.h>

#include <vector>
#include <ctime>

static const std::vector<State> DATA_LIGHT = {
    // goldstone meson (c c-bar) 17
    {1, 1, 1, 1, 0, 1, 775.3,    5},   // rho(770)
    {1, 1, 1, 1, 1, 2, 1318.2,   5},   // a2(1320)
    {1, 2, 1, 0, 0, 0, 497.6,    5},   // K(1S)
    {1, 2, 1, 1, 0, 1, 895.6,    5},   // K*(1S)
    {1, 2, 1, 1, 1, 2, 1432.4,   5},   // K2*(1430)
    {2, 2, 1, 1, 0, 1, 1019.5,   5},   // phi(1S)
    {2, 2, 1, 1, 1, 2, 1517.3,   5},   // f2'(1525)
};

void minuit2_light_GIScreen(double *params_out)
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
    upar.Add("b", 0.18, 0.01, 0.1, 0.3);
    upar.Add("mu", 0.15, 0.01, 0.1, 0.2);
    upar.Add("c", -0.253, 0.01, -2.0, 0.0);
    upar.Add("sig0", 1.8, 0.01, 1.0, 3.0);
    upar.Add("s", 1.55, 0.01, 1.0, 3.0);
    upar.Add("econt", -0.168, 0.01, -0.5, 0.0);
    upar.Add("esov", -0.035, 0.01, -1.0, 1.0);
    upar.Add("esos", 0.055, 0.01, -1.0, 1.0);
    upar.Add("etens", 0.025, 0.01, -1.0, 1.0);
    upar.Fix("mc");
    upar.Fix("mb");
    int N_PARAMS = upar.Params().size();

    /* use of Migrad algorithm with strategy 2, high precision */
    Chi2Minimizer minuit_fit(DATA_LIGHT, MODEL_GI_SCREEN, N_PARAMS, 1.0);
    ROOT::Minuit2::MnMigrad migrad(minuit_fit, upar, 2);

    /* perform the fit */
    ROOT::Minuit2::FunctionMinimum min_result = migrad();
    compute_chi2(DATA_LIGHT, min_result.UserParameters().Params(), MODEL_GI_SCREEN, true);
    dual << min_result.UserParameters() << std::endl;

    auto params = min_result.UserParameters().Params();
    for (int i = 0; i < N_PARAMS; i++) {
        params_out[i] = params[i];
    }
}
/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/param/fbsbar.h>
#include <gemstore/param/minuit.h>
#include <gemstore/types.h>

#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnPrint.h>

#include <vector>
#include <ctime>

static const std::vector<State> DATA_BSBAR = {
    // Bs (s b-bar) 4
    {2, 4, 1, 0, 0, 0, 5366.9,   1},   // Bs(1^1S0)
    {2, 4, 1, 1, 0, 1, 5415.4,   1},   // Bs*(1^3S1)
    {2, 4, 1, 0, 1, 1, 5828.7,   1},   // Bs1(5830), mixed 1^1P1/1^3P1, stored here as dominant ^1P1-like state
    {2, 4, 1, 1, 1, 2, 5839.9,   1},   // Bs2*(5840), 1^3P2
};

void minuit2_bsbar_GIScreen(double *params_out)
{
    srand(time(0));
    DualStream dual("Fitting.out");

    /* set parameters */
    ROOT::Minuit2::MnUserParameters upar;
    //upar.Add(name, value, init_step, lower_limit, upper_limit);
    upar.Add("mn", 0.4560806209112, 0.01, 0.1, 0.5);
    upar.Add("ms", 0.6173440068792, 0.01, 0.3, 0.7);
    upar.Add("mc", 1.805387067165, 0.01, 1.5, 2.0);
    upar.Add("mb", 5.151269542382, 0.01, 4.5, 5.5);
    upar.Add("b", 0.2522010221331, 0.01, 0.1, 0.3);
    upar.Add("mu", 0.15600, 0.01, 0.1, 0.2);
    upar.Add("c", -0.6482214863381, 0.01, -2.0, 0.0);
    upar.Add("sig0", 1.770545357386, 0.01, 1.0, 3.0);
    upar.Add("s", 1.146340880694, 0.01, 1.0, 3.0);
    upar.Add("econt", -0.3199525941084, 0.01, -0.5, 0.0);
    upar.Add("esov", -0.3428467195716, 0.01, -1.0, 1.0);
    upar.Add("esos", 0.7905135472165, 0.01, -1.0, 1.0);
    upar.Add("etens", -0.5000469699086, 0.01, -1.0, 1.0);
    upar.Fix("mn");
    upar.Fix("ms");
    upar.Fix("mc");
    upar.Fix("mb");
    upar.Fix("mu");
    upar.Fix("c");
    upar.Fix("sig0");
    upar.Fix("s");
    int N_PARAMS = upar.Params().size();

    /* use of Migrad algorithm with strategy 2, high precision */
    Chi2Minimizer minuit_fit(DATA_BSBAR, MODEL_GISCREEN, N_PARAMS, 1.0);
    ROOT::Minuit2::MnMigrad migrad(minuit_fit, upar, 2);

    /* perform the fit */
    ROOT::Minuit2::FunctionMinimum min_result = migrad();
    compute_chi2(DATA_BSBAR, min_result.UserParameters().Params(), MODEL_GISCREEN, true);
    dual << min_result.UserParameters() << std::endl;

    auto params = min_result.UserParameters().Params();
    for (int i = 0; i < N_PARAMS; i++) {
        params_out[i] = params[i];
    }
}
/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/param/fscbar.h>
#include <gemstore/param/minuit.h>
#include <gemstore/types.h>

#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnPrint.h>

#include <vector>
#include <ctime>

static const std::vector<State> DATA_SCBAR = {
    // Ds (s c-bar) 10
    {2, 3, 1, 0, 0, 0, 1968.4,   1},    // Ds(1^1S0)
    {2, 3, 1, 1, 0, 1, 2112.2,   1},    // Ds*(1^3S1)
    {2, 3, 1, 1, 1, 0, 2317.8,   1},    // Ds0*(2317), candidate 1^3P0
    {2, 3, 1, 1, 1, 1, 2459.5,   1},    // Ds1(2460), mixed 1^1P1/1^3P1, stored here as dominant ^3P1-like state
    {2, 3, 1, 0, 1, 1, 2535.1,   1},    // Ds1(2536), 1^1P1
    {2, 3, 1, 1, 1, 2, 2569.1,   1},    // Ds2*(2573), 1^3P2
    {2, 3, 2, 1, 0, 1, 2714.0,   5},    // Ds1*(2700), mainly 2^3S1
                                        // sensitive: Song et al. favor 2^3S1-1^3D1 mixing, but this basis stores the dominant component
    {2, 3, 1, 1, 2, 1, 2859.0,  27},    // Ds1*(2860), 1^3D1
                                        // candidate identified after LHCb amplitude analysis
    {2, 3, 1, 1, 2, 3, 2860.5,   7},    // Ds3*(2860), 1^3D3
    {2, 3, 2, 1, 1, 1, 3044.0,  31},    // DsJ(3040), good candidate for 2P(1+)
                                        // sensitive: mixed 2^1P1/2^3P1 state; Song et al. favor this assignment, but alternative interpretations are discussed
};

void minuit2_scbar_GIScreen(double *params_out)
{
    srand(time(0));
    DualStream dual("Fitting.out");

    /* set parameters */
    ROOT::Minuit2::MnUserParameters upar;
    //upar.Add(name, value, init_step, lower_limit, upper_limit);
    upar.Add("mn", 0.3349266291038, 0.01, 0.1, 0.5);
    upar.Add("ms", 0.5283844975353, 0.01, 0.3, 0.7);
    upar.Add("mc", 1.747603574365, 0.01, 1.5, 2.0);
    upar.Add("mb", 5.095838715, 0.01, 4.5, 5.5);
    upar.Add("b", 0.248247135518, 0.01, 0.1, 0.3);
    upar.Add("mu", 0.153, 0.01, 0.14, 0.16);
    upar.Add("c", -0.5334999044266, 0.01, -2.0, 0.0);
    upar.Add("sig0", 1.56552865791, 0.01, 1.0, 3.0);
    upar.Add("s", 1.285723132711, 0.01, 1.0, 3.0);
    upar.Add("econt", -0.2864647624566, 0.01, -0.5, 0.0);
    upar.Add("esov", -0.349573212139, 0.01, -1.0, 1.0);
    upar.Add("esos", 0.7905135472165, 0.01, -1.0, 1.0);
    upar.Add("etens", -0.487322874302, 0.01, -1.0, 1.0);
    upar.Fix("mn");
    upar.Fix("ms");
    upar.Fix("mc");
    upar.Fix("mb");
    upar.Fix("c");
    upar.Fix("mu");
    int N_PARAMS = upar.Params().size();

    /* use of Migrad algorithm with strategy 2, high precision */
    Chi2Minimizer minuit_fit(DATA_SCBAR, MODEL_GISCREEN, N_PARAMS, 1.0);
    ROOT::Minuit2::MnMigrad migrad(minuit_fit, upar, 2);

    /* perform the fit */
    ROOT::Minuit2::FunctionMinimum min_result = migrad();
    compute_chi2(DATA_SCBAR, min_result.UserParameters().Params(), MODEL_GISCREEN, true);
    dual << min_result.UserParameters() << std::endl;

    auto params = min_result.UserParameters().Params();
    for (int i = 0; i < N_PARAMS; i++) {
        params_out[i] = params[i];
    }
}
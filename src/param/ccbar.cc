/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/param/ccbar.h>
#include <gemstore/param/minuit.h>
#include <gemstore/types.h>

#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnPrint.h>

#include <vector>
#include <ctime>

static const std::vector<State> DATA_CCBAR = {
    // charmonium (c c-bar) 17
    {3, 3, 1, 0, 0, 0, 2984.1,   1},   // ηc(1S)
    {3, 3, 2, 0, 0, 0, 3637.8,   1},   // ηc(2S)
    {3, 3, 1, 1, 0, 1, 3096.9,   1},   // ψ(1S)
    {3, 3, 2, 1, 0, 1, 3686.1,   1},   // ψ(2S)
    {3, 3, 3, 1, 0, 1, 4039.6,   4},   // ψ(3S)
    {3, 3, 4, 1, 0, 1, 4222.2,   3},   // ψ(4S)
    {3, 3, 5, 1, 0, 1, 4415.0,   5},   // ψ(5S)
    {3, 3, 1, 1, 2, 1, 3778.1,   1},   // ψ(1D)
    {3, 3, 1, 1, 2, 2, 3823.5,   1},   // ψ2(1D)
    {3, 3, 1, 1, 2, 3, 3842.7,   1},   // ψ3(1D)
    {3, 3, 2, 1, 2, 1, 4191.0,   5},   // ψ(2D)
    {3, 3, 3, 1, 2, 1, 4374.0,   7},   // ψ(3D)
    {3, 3, 1, 0, 1, 1, 3525.4,   1},   // hc(1P)
    {3, 3, 1, 1, 1, 0, 3414.7,   1},   // χc0(1P)
    {3, 3, 1, 1, 1, 1, 3510.7,   1},   // χc1(1P)
    {3, 3, 1, 1, 1, 2, 3556.2,   1},   // χc2(1P)
    {3, 3, 2, 1, 1, 2, 3922.5,   1},   // χc2(2P)
};

void minuit2_ccbar_GIScreen(double *params_out)
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
    upar.Add("mu", 0.1333931469096, 0.01, 0.1, 0.2);
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
    upar.Fix("sig0");
    upar.Fix("s");
    int N_PARAMS = upar.Params().size();

    /* use of Migrad algorithm with strategy 2, high precision */
    Chi2Minimizer minuit_fit(DATA_CCBAR, MODEL_GISCREEN, N_PARAMS, 1.0);
    ROOT::Minuit2::MnMigrad migrad(minuit_fit, upar, 2);

    /* perform the fit */
    ROOT::Minuit2::FunctionMinimum min_result = migrad();
    compute_chi2(DATA_CCBAR, min_result.UserParameters().Params(), MODEL_GISCREEN, true);
    dual << min_result.UserParameters() << std::endl;

    auto params = min_result.UserParameters().Params();
    for (int i = 0; i < N_PARAMS; i++) {
        params_out[i] = params[i];
    }
}
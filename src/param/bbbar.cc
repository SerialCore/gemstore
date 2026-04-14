/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/param/bbbar.h>
#include <gemstore/param/minuit.h>
#include <gemstore/types.h>

#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnPrint.h>

#include <vector>
#include <ctime>

static const std::vector<State> DATA_BBBAR = {
    // bottomonium (b b-bar) 19
    {4, 4, 1, 0, 0, 0, 9398.7,   2},   // ηb(1S)
    {4, 4, 2, 0, 0, 0, 9999.0,   4},   // ηb(2S)
    {4, 4, 1, 1, 0, 1, 9460.4,   1},   // Υ(1S)
    {4, 4, 2, 1, 0, 1, 10023.4,  1},   // Υ(2S)
    {4, 4, 3, 1, 0, 1, 10355.1,  1},   // Υ(3S)
    {4, 4, 4, 1, 0, 1, 10579.4,  2},   // Υ(4S)
    {4, 4, 5, 1, 0, 1, 10885.2,  3},   // Υ(5S)
    {4, 4, 6, 1, 0, 1, 11000.0,  4},   // Υ(6S)
    {4, 4, 1, 1, 2, 2, 10163.7,  2},   // Υ2(1D)
    {4, 4, 1, 0, 1, 1, 9899.3,   1},   // hb(1P)
    {4, 4, 2, 0, 1, 1, 10259.8,  2},   // hb(2P)
    {4, 4, 1, 1, 1, 0, 9859.4,   1},   // χb0(1P)
    {4, 4, 1, 1, 1, 1, 9892.8,   1},   // χb1(1P)
    {4, 4, 1, 1, 1, 2, 9912.2,   1},   // χb2(1P)
    {4, 4, 2, 1, 1, 0, 10232.5,  1},   // χb0(2P)
    {4, 4, 2, 1, 1, 1, 10255.5,  1},   // χb1(2P)
    {4, 4, 2, 1, 1, 2, 10268.7,  1},   // χb2(2P)
    {4, 4, 3, 1, 1, 1, 10513.4,  1},   // χb1(3P)
    {4, 4, 3, 1, 1, 2, 10524.0,  1}    // χb2(3P)
};

void minuit2_bbbar_GIScreen(double *params_out)
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
    upar.Add("b1", 0.248247135518, 0.01, 0.1, 0.3);
    upar.Add("mu", 0.1333931469096, 0.01, 0.05, 0.2);
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
    Chi2Minimizer minuit_fit(DATA_BBBAR, MODEL_GI_SCREEN, N_PARAMS, 1.0);
    ROOT::Minuit2::MnMigrad migrad(minuit_fit, upar, 2);

    /* perform the fit */
    ROOT::Minuit2::FunctionMinimum min_result = migrad();
    compute_chi2(DATA_BBBAR, min_result.UserParameters().Params(), MODEL_GI_SCREEN, true);
    dual << min_result.UserParameters() << std::endl;

    auto params = min_result.UserParameters().Params();
    for (int i = 0; i < N_PARAMS; i++) {
        params_out[i] = params[i];
    }
}

void minuit2_bbbar_GIQuadra(double *params_out)
{
    srand(time(0));
    DualStream dual("Fitting.out");

    /* set parameters */
    ROOT::Minuit2::MnUserParameters upar;
    //upar.Add(name, value, init_step, lower_limit, upper_limit);
    upar.Add("mn", 0.4024291451329, 0.01, 0.1, 0.5);
    upar.Add("ms", 0.5749451022621, 0.01, 0.3, 0.7);
    upar.Add("mc", 1.773779650685, 0.01, 1.5, 2.0);
    upar.Add("mb", 5.119870696908, 0.01, 4.5, 5.5);
    upar.Add("b1", 0.2020007923859, 0.01, 0.1, 0.3);
    upar.Add("b2", 0.01106920061884, 0.01, 0.0, 0.1);
    upar.Add("mu", 0.1079930218451, 0.01, 0.05, 0.2);
    upar.Add("c", -0.5485460827658, 0.01, -2.0, 0.0);
    upar.Add("sig0", 1.700037337287, 0.01, 1.0, 3.0);
    upar.Add("s", 1.288608951331, 0.01, 1.0, 3.0);
    upar.Add("econt", -0.273180413884, 0.01, -0.5, 0.0);
    upar.Add("esov", -0.7373105432325, 0.01, -1.0, 1.0);
    upar.Add("esos", 0.9890025702377, 0.01, -1.0, 1.0);
    upar.Add("etens", -0.5045346358157, 0.01, -1.0, 1.0);
    upar.Fix("mn");
    upar.Fix("ms");
    upar.Fix("mc");
    upar.Fix("mb");
    upar.Fix("c");
    upar.Fix("sig0");
    upar.Fix("s");
    int N_PARAMS = upar.Params().size();

    /* use of Migrad algorithm with strategy 2, high precision */
    Chi2Minimizer minuit_fit(DATA_BBBAR, MODEL_GI_QUADRA, N_PARAMS, 1.0);
    ROOT::Minuit2::MnMigrad migrad(minuit_fit, upar, 2);

    /* perform the fit */
    ROOT::Minuit2::FunctionMinimum min_result = migrad();
    compute_chi2(DATA_BBBAR, min_result.UserParameters().Params(), MODEL_GI_QUADRA, true);
    dual << min_result.UserParameters() << std::endl;

    auto params = min_result.UserParameters().Params();
    for (int i = 0; i < N_PARAMS; i++) {
        params_out[i] = params[i];
    }
}
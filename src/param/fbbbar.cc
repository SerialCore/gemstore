/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/param/fbbbar.h>
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
    upar.Add("mn", 0.4560806209112, 0.01, 0.1, 0.5);
    upar.Add("ms", 0.6173440068792, 0.01, 0.3, 0.7);
    upar.Add("mc", 1.805387067165, 0.01, 1.5, 2.0);
    upar.Add("mb", 5.151269542382, 0.01, 4.5, 5.5);
    upar.Add("b", 0.2522010221331, 0.01, 0.1, 0.3);
    upar.Add("mu", 0.1345502083599, 0.01, 0.1, 0.2);
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
    upar.Fix("c");
    upar.Fix("sig0");
    upar.Fix("s");
    int N_PARAMS = upar.Params().size();

    /* use of Migrad algorithm with strategy 2, high precision */
    Chi2Minimizer minuit_fit(DATA_BBBAR, MODEL_GISCREEN, N_PARAMS, 1.0);
    ROOT::Minuit2::MnMigrad migrad(minuit_fit, upar, 2);

    /* perform the fit */
    ROOT::Minuit2::FunctionMinimum min_result = migrad();
    compute_chi2(DATA_BBBAR, min_result.UserParameters().Params(), MODEL_GISCREEN, true);
    dual << min_result.UserParameters() << std::endl;

    auto params = min_result.UserParameters().Params();
    for (int i = 0; i < N_PARAMS; i++) {
        params_out[i] = params[i];
    }
}
/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/param/ccbar.h>
#include <gemstore/param/helper.h>
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
    upar.Add("mn", 0.4433275191676, 0.01, 0.1, 0.5);
    upar.Add("ms", 0.606437654085, 0.01, 0.3, 0.7);
    upar.Add("mc", 1.797694082834, 0.01, 1.5, 2.0);
    upar.Add("mb", 5.142643566529, 0.01, 4.5, 5.5);
    upar.Add("b1", 0.2530091124005, 0.01, 0.1, 0.3);
    upar.Add("mu", 0.1401631404922, 0.01, 0.1, 0.2);
    upar.Add("c", -0.6300631479611, 0.01, -2.0, 0.0);
    upar.Add("sig0", 1.776663048762, 0.01, 1.0, 3.0);
    upar.Add("s", 1.180712032887, 0.01, 1.0, 3.0);
    upar.Add("econt", -0.3068965865533, 0.01, -0.5, 0.0);
    upar.Add("esov", -0.3715004921572, 0.01, -1.0, 1.0);
    upar.Add("esos", 0.9276232568028, 0.01, -1.0, 1.0);
    upar.Add("etens", 0.5056764482414, 0.01, -1.0, 1.0);
    upar.Fix("mn");
    upar.Fix("ms");
    upar.Fix("mc");
    upar.Fix("mb");
    upar.Fix("sig0");
    upar.Fix("s");
    int N_PARAMS = upar.Params().size();

    /* use of Migrad algorithm with strategy 2, high precision */
    Chi2Minimizer minuit_fit(DATA_CCBAR, MODEL_GI_SCREEN, N_PARAMS, 1.0);
    ROOT::Minuit2::MnMigrad migrad(minuit_fit, upar, 2);

    /* perform the fit */
    ROOT::Minuit2::FunctionMinimum min_result = migrad();
    compute_chi2(DATA_CCBAR, min_result.UserParameters().Params(), MODEL_GI_SCREEN, true);
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
    upar.Add("mn", 0.2722163527103, 0.01, 0.1, 0.5);
    upar.Add("ms", 0.4947009852984, 0.01, 0.3, 0.7);
    upar.Add("mc", 1.731560039229, 0.01, 1.5, 2.0);
    upar.Add("mb", 5.079977591962, 0.01, 4.5, 5.5);
    upar.Add("b1", 0.2047731657768, 0.01, 0.1, 0.3);
    upar.Add("b2", 0.0116897321522, 0.01, 0.0, 0.1);
    upar.Add("mu", 0.1149494963463, 0.01, 0.1, 0.2);
    upar.Add("c", -0.473796639106, 0.01, -2.0, 0.0);
    upar.Add("sig0", 1.405806263295, 0.01, 1.0, 3.0);
    upar.Add("s", 1.355912150877, 0.01, 1.0, 3.0);
    upar.Add("econt", -0.3097205997153, 0.01, -0.5, 0.0);
    upar.Add("esov", -0.309578055483, 0.01, -1.0, 1.0);
    upar.Add("esos", 0.9476943481252, 0.01, -1.0, 1.0);
    upar.Add("etens", 0.5403018589044, 0.01, -1.0, 1.0);
    upar.Fix("mn");
    upar.Fix("ms");
    upar.Fix("mc");
    upar.Fix("mb");
    upar.Fix("sig0");
    upar.Fix("s");
    int N_PARAMS = upar.Params().size();

    /* use of Migrad algorithm with strategy 2, high precision */
    Chi2Minimizer minuit_fit(DATA_CCBAR, MODEL_GI_QUADRA, N_PARAMS, 1.0);
    ROOT::Minuit2::MnMigrad migrad(minuit_fit, upar, 2);

    /* perform the fit */
    ROOT::Minuit2::FunctionMinimum min_result = migrad();
    compute_chi2(DATA_CCBAR, min_result.UserParameters().Params(), MODEL_GI_QUADRA, true);
    dual << min_result.UserParameters() << std::endl;

    auto params = min_result.UserParameters().Params();
    for (int i = 0; i < N_PARAMS; i++) {
        params_out[i] = params[i];
    }
}
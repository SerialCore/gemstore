/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/param/argset.h>

#include <gemstore/types.h>
#include <gemstore/parse.h>

const argsGIModel_t argsGIString_meson = {
    .mn = 0.220,
    .ms = 0.419,
    .mc = 1.628,
    .mb = 4.977,
    .b = 0.18,
    .c = -0.253,
    .sigma_0 = 1.8,
    .s = 1.55,
    .epsilon_cont = -0.168,
    .epsilon_sov = -0.035,
    .epsilon_sos = 0.055,
    .epsilon_tens = 0.025
};

const argsGIModel_t argsGIScreen_meson = {
    .mn = 0.3349266291038,
    .ms = 0.5283844975353,
    .mc = 1.747603574365,
    .mb = 5.095838715,
    .b = 0.248247135518,
    .mu = 0.1333931469096,
    .c = -0.5334999044266,
    .sigma_0 = 1.56552865791,
    .s = 1.285723132711,
    .epsilon_cont = -0.2864647624566,
    .epsilon_sov = -0.349573212139,
    .epsilon_sos = 0.7905135472165,
    .epsilon_tens = -0.487322874302,
};

const argsGIModel_t argsGIScreen_bbbar = {
    .mn = 0.3349266291038,
    .ms = 0.5283844975353,
    .mc = 1.747603574365,
    .mb = 5.095838715,
    .b = 0.2443022498196,
    .mu = 0.1241517438531,
    .c = -0.5334999044266,
    .sigma_0 = 1.56552865791,
    .s = 1.285723132711,
    .epsilon_cont = -0.4998459653183,
    .epsilon_sov = -0.8599140847974,
    .epsilon_sos = -0.49970641935,
    .epsilon_tens = -0.7579113279148,
};

const argsGIModel_t argsGIScreen_ccbar = {
    .mn = 0.3349266291038,
    .ms = 0.5283844975353,
    .mc = 1.747603574365,
    .mb = 5.095838715,
    .b = 0.2518030121301,
    .mu = 0.1417429469983,
    .c = -0.5334999044266,
    .sigma_0 = 1.56552865791,
    .s = 1.285723132711,
    .epsilon_cont = -0.2888523192093,
    .epsilon_sov = -0.3389117883675,
    .epsilon_sos = 0.9999999547283,
    .epsilon_tens = -0.499991780808,
};

argsGIModel_t argsGIModel_from(const argsInput_t *input)
{
    argsGIModel_t args_model = {0};

    if (input->param == PARAM_GISTRING_MESON) args_model = argsGIString_meson;
    else if (input->param == PARAM_GISCREEN_MESON) args_model = argsGIScreen_meson;
    else if (input->param == PARAM_GISCREEN_BBBAR) args_model = argsGIScreen_bbbar;
    else if (input->param == PARAM_GISCREEN_CCBAR) args_model = argsGIScreen_ccbar;
    else if (input->param == PARAM_GISTRING_CUSTOM) {
        parse_param_GISTRING(input->param_file, &args_model);
    }
    else if (input->param == PARAM_GISCREEN_CUSTOM) {
        parse_param_GISCREEN(input->param_file, &args_model);
    }

    return args_model;
}
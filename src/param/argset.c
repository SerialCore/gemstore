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
    .mn = 0.4560806209112,
    .ms = 0.6173440068792,
    .mc = 1.805387067165,
    .mb = 5.151269542382,
    .b = 0.2522010221331,
    .mu = 0.1345502083599,
    .c = -0.6482214863381,
    .sigma_0 = 1.770545357386,
    .s = 1.146340880694,
    .epsilon_cont = -0.3199525941084,
    .epsilon_sov = -0.3428467195716,
    .epsilon_sos = 0.9999999905835,
    .epsilon_tens = -0.5000469699086,
};

const argsGIModel_t argsGIScreen_bbbar = {
    .mn = 0.4560806209112,
    .ms = 0.6173440068792,
    .mc = 1.805387067165,
    .mb = 5.151269542382,
    .b = 0.2469561138537,
    .mu = 0.1221612749279,
    .c = -0.6482214863381,
    .sigma_0 = 1.770545357386,
    .s = 1.146340880694,
    .epsilon_cont = -0.4982927251991,
    .epsilon_sov = -0.1461227889302,
    .epsilon_sos = -0.5018806829405,
    .epsilon_tens = -0.8456374804674,
};

const argsGIModel_t argsGIScreen_bcbar = {
    .mn = 0.4560806209112,
    .ms = 0.6173440068792,
    .mc = 1.805387067165,
    .mb = 5.151269542382,
    .b = 0.2554804581,
    .mu = 0.13694,
    .c = -0.6482214863381,
    .sigma_0 = 1.770545357386,
    .s = 1.146340880694,
    .epsilon_cont = -0.499130628,
    .epsilon_sov = -0.978788257,
    .epsilon_sos = -0.7682287286,
    .epsilon_tens = -0.2737721829,
};

const argsGIModel_t argsGIScreen_bsbar = {
    .mn = 0.4560806209112,
    .ms = 0.6173440068792,
    .mc = 1.805387067165,
    .mb = 5.151269542382,
    .b = 0.2562164758,
    .mu = 0.14610,
    .c = -0.6482214863381,
    .sigma_0 = 1.770545357386,
    .s = 1.146340880694,
    .epsilon_cont = -0.292816925,
    .epsilon_sov = 0.223940353,
    .epsilon_sos = 0.7325063004,
    .epsilon_tens = -0.5195779148,
};

const argsGIModel_t argsGIScreen_ccbar = {
    .mn = 0.4560806209112,
    .ms = 0.6173440068792,
    .mc = 1.805387067165,
    .mb = 5.151269542382,
    .b = 0.256309879544,
    .mu = 0.1440535603044,
    .c = -0.6482214863381,
    .sigma_0 = 1.770545357386,
    .s = 1.146340880694,
    .epsilon_cont = -0.3585942884847,
    .epsilon_sov = -0.4985927372035,
    .epsilon_sos = 0.9999999949382,
    .epsilon_tens = -0.4999999166133,
};

const argsGIModel_t argsGIScreen_csbar = {
    .mn = 0.4560806209112,
    .ms = 0.6173440068792,
    .mc = 1.805387067165,
    .mb = 5.151269542382,
    .b = 0.252102902,
    .mu = 0.14725,
    .c = -0.6482214863381,
    .sigma_0 = 1.770545357386,
    .s = 1.146340880694,
    .epsilon_cont = -0.2629809578,
    .epsilon_sov = -0.499998198,
    .epsilon_sos = 0.9999999652,
    .epsilon_tens = -0.4999997792,
};

argsGIModel_t argsGIModel_from(const argsInput_t *input)
{
    argsGIModel_t args_model = {0};

    if (input->param == PARAM_GISTRING_MESON) args_model = argsGIString_meson;
    else if (input->param == PARAM_GISCREEN_MESON) args_model = argsGIScreen_meson;
    else if (input->param == PARAM_GISCREEN_BBBAR) args_model = argsGIScreen_bbbar;
    else if (input->param == PARAM_GISCREEN_BCBAR) args_model = argsGIScreen_bcbar;
    else if (input->param == PARAM_GISCREEN_BSBAR) args_model = argsGIScreen_bsbar;
    else if (input->param == PARAM_GISCREEN_CCBAR) args_model = argsGIScreen_ccbar;
    else if (input->param == PARAM_GISCREEN_CSBAR) args_model = argsGIScreen_csbar;
    else if (input->param == PARAM_GISTRING_CUSTOM) {
        parse_param_GISTRING(input->param_file, &args_model);
    }
    else if (input->param == PARAM_GISCREEN_CUSTOM) {
        parse_param_GISCREEN(input->param_file, &args_model);
    }

    return args_model;
}
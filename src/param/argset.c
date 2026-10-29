/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/param/argset.h>

#include <gemstore/types.h>
#include <gemstore/parse.h>

#include <stdio.h>
#include <stdlib.h>

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

const argsGIModel_t argsGIString_baryon = {
    .mn = 0.220,
    .ms = 0.419,
    .mc = 1.628,
    .mb = 4.977,
    .b = 0.141,
    .c = -0.204,
    .sigma_0 = 1.889,
    .s = 1.422,
    .epsilon_cont = -0.156,
    .epsilon_tens = 0.379,
    .epsilon_sov = 0.006,
    .epsilon_sos = 0.449
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
    .b = 0.2549474924,
    .mu = 0.13561,
    .c = -0.6482214863381,
    .sigma_0 = 1.770545357386,
    .s = 1.146340880694,
    .epsilon_cont = -0.4993633033,
    .epsilon_sov = -0.3428460826,
    .epsilon_sos = 0.7914990349,
    .epsilon_tens = -0.500062988,
};

const argsGIModel_t argsGIScreen_bsbar = {
    .mn = 0.4560806209112,
    .ms = 0.6173440068792,
    .mc = 1.805387067165,
    .mb = 5.151269542382,
    .b = 0.2594958766,
    .mu = 0.15600,
    .c = -0.6482214863381,
    .sigma_0 = 1.770545357386,
    .s = 1.146340880694,
    .epsilon_cont = -0.3202159303,
    .epsilon_sov = -0.2947420493,
    .epsilon_sos = -0.1691792532,
    .epsilon_tens = -0.9995011891,
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
    .b = 0.2575524695,
    .mu = 0.16255,
    .c = -0.6482214863381,
    .sigma_0 = 1.770545357386,
    .s = 1.146340880694,
    .epsilon_cont = -0.2795867762,
    .epsilon_sov = -0.4999906411,
    .epsilon_sos = 0.9999999437,
    .epsilon_tens = -0.5000014157,
};

const argsNRModel_t argsNRString_meson = {
    .mn = 0.606,
    .ms = 0.780,
    .mc = 1.984,
    .mb = 5.368,
    .b = 0.2312,
    .c = -1.1711,
    .alpha_s = 0.3930,
    .sigma = 1.842
};

const argsNRModel_t argsNRScreen_meson = {
    .mn = 0.606,
    .ms = 0.780,
    .mc = 1.984,
    .mb = 5.368,
    .b = 0.2312,
    .mu = 0.0690,
    .c = -1.1711,
    .alpha_s = 0.3930,
    .sigma = 1.842
};

const argsBAGModel_t argsBAGMIT_default = {
    .mn = 0.000,
    .ms = 0.279,
    .mc = 1.641,
    .mb = 5.093,
    .Bcs = -0.0246,
    .Bbs = -0.0319,
    .Bcc = -0.0772,
    .Bbc = -0.1013,
    .Bbb = -0.1279,
    .B = 0.000442051,
    .Z = 1.83
};

argsGIModel_t argsGIModel_from(const argsInput_t *input)
{
    argsGIModel_t args_model = {0};

    if (input->param == PARAM_GISTRING_MESON) args_model = argsGIString_meson;
    else if (input->param == PARAM_GISTRING_BARYON) args_model = argsGIString_baryon;
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

argsNRModel_t argsNRModel_from(const argsInput_t *input)
{
    argsNRModel_t args_model = {0};

    if (input->param == PARAM_NRSTRING_MESON) args_model = argsNRString_meson;
    else if (input->param == PARAM_NRSCREEN_MESON) args_model = argsNRScreen_meson;
    else if (input->param == PARAM_NRSTRING_CUSTOM) {
        parse_param_NRSTRING(input->param_file, &args_model);
    }
    else if (input->param == PARAM_NRSCREEN_CUSTOM) {
        parse_param_NRSCREEN(input->param_file, &args_model);
    }

    return args_model;
}

argsBAGModel_t argsBAGModel_from(const argsInput_t *input)
{
    argsBAGModel_t args_model = {0};

    if (input->param == PARAM_BAGMIT_DEFAULT) args_model = argsBAGMIT_default;
    else if (input->param == PARAM_BAGMIT_CUSTOM) {
        parse_param_BAGMIT(input->param_file, &args_model);
    }

    return args_model;
}

argsOBEModel_t argsOBEModel_from(const argsInput_t *input)
{
    argsOBEModel_t args_model = {0};

    parse_param_OBESTAT(input->param_file, &args_model);
    args_model.type1 = input->type1;
    args_model.type2 = input->type2;
    args_model.cc1 = input->cc1;
    args_model.cc2 = input->cc2;
    args_model.eta = input->eta;
    args_model.iso_pi = input->iso_pi;
    args_model.iso_rho = input->iso_rho;

    return args_model;
}
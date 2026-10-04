/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/param/fitting.h>

static const fit_state_t k_states[] = {
    // Bc (c b-bar) 2
    {"Bc(1S)", 3, 4, 1, 0, 0, 0, 6274.5, 1},        /* Bc(1S) */
    {"Bc*(1^3S1)", 3, 4, 1, 1, 0, 1, 6339.0, 2},    /* Bc*(1^3S1) */
    {"Bc(2S)", 3, 4, 2, 0, 0, 0, 6871.2, 1},        /* Bc(2S) */
};

static const fit_param_t k_params[] = {
    {"mn", 0.4560806209112, 0.01, 0.1, 0.5, 1},
    {"ms", 0.6173440068792, 0.01, 0.3, 0.7, 1},
    {"mc", 1.805387067165, 0.01, 1.5, 2.0, 1},
    {"mb", 5.151269542382, 0.01, 4.5, 5.5, 1},
    {"b", 0.2522010221331, 0.01, 0.1, 0.3, 0},
    {"mu", 0.13561, 0.01, 0.1, 0.2, 1},
    {"c", -0.6482214863381, 0.01, -2.0, 0.0, 1},
    {"sigma_0", 1.770545357386, 0.01, 1.0, 3.0, 1},
    {"s", 1.146340880694, 0.01, 1.0, 3.0, 1},
    {"epsilon_cont", -0.3199525941084, 0.01, -0.5, 0.0, 0},
    {"epsilon_sov", -0.3428467195716, 0.01, -1.0, 1.0, 0},
    {"epsilon_sos", 0.7905135472165, 0.01, -1.0, 1.0, 0},
    {"epsilon_tens", -0.5000469699086, 0.01, -1.0, 1.0, 0},
};

const fit_target_t fit_target_giscreen_bcbar = {
    MODEL_GISCREEN,
    k_states,
    (int)(sizeof k_states / sizeof k_states[0]),
    k_params,
    (int)(sizeof k_params / sizeof k_params[0])
};

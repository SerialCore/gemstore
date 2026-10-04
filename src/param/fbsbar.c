/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/param/fitting.h>

static const fit_state_t k_states[] = {
    // Bs (s b-bar) 4
    {"Bs(1^1S0)", 2, 4, 1, 0, 0, 0, 5366.9, 1},     /* Bs(1^1S0) */
    {"Bs*(1^3S1)", 2, 4, 1, 1, 0, 1, 5415.4, 1},    /* Bs*(1^3S1) */
    {"Bs1(5830)", 2, 4, 1, 0, 1, 1, 5828.7, 1},     /* Bs1(5830), mixed 1^1P1/1^3P1, stored here as dominant ^1P1-like state */
    {"Bs2*(5840)", 2, 4, 1, 1, 1, 2, 5839.9, 1},    /* Bs2*(5840), 1^3P2 */
};

static const fit_param_t k_params[] = {
    {"mn", 0.4560806209112, 0.01, 0.1, 0.5, 1},
    {"ms", 0.6173440068792, 0.01, 0.3, 0.7, 1},
    {"mc", 1.805387067165, 0.01, 1.5, 2.0, 1},
    {"mb", 5.151269542382, 0.01, 4.5, 5.5, 1},
    {"b", 0.2522010221331, 0.01, 0.1, 0.3, 0},
    {"mu", 0.15600, 0.01, 0.1, 0.2, 1},
    {"c", -0.6482214863381, 0.01, -2.0, 0.0, 1},
    {"sigma_0", 1.770545357386, 0.01, 1.0, 3.0, 1},
    {"s", 1.146340880694, 0.01, 1.0, 3.0, 1},
    {"epsilon_cont", -0.3199525941084, 0.01, -0.5, 0.0, 0},
    {"epsilon_sov", -0.3428467195716, 0.01, -1.0, 1.0, 0},
    {"epsilon_sos", 0.7905135472165, 0.01, -1.0, 1.0, 0},
    {"epsilon_tens", -0.5000469699086, 0.01, -1.0, 1.0, 0},
};

const fit_target_t fit_target_giscreen_bsbar = {
    MODEL_GISCREEN,
    k_states,
    (int)(sizeof k_states / sizeof k_states[0]),
    k_params,
    (int)(sizeof k_params / sizeof k_params[0])
};

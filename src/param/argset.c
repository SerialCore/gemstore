/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/param/argset.h>

const argsGIModel_t argsGIString_meson = {
    .mn = 0.220,
    .ms = 0.419,
    .mc = 1.628,
    .mb = 4.977,
    .b1 = 0.18,
    .c = -0.253,
    .sigma_0 = 1.8,
    .s = 1.55,
    .epsilon_cont = -0.168,
    .epsilon_sov = -0.035,
    .epsilon_sos = 0.055,
    .epsilon_tens = 0.025
};

const argsGIModel_t argsGIScreen_meson = {
    .mn = 0.4713455847642,
    .ms = 0.6283121820133,
    .mc = 1.810505119204,
    .mb = 5.156014766761,
    .b1 = 0.2532976556232,
    .mu = 0.1356286583032,
    .c = -0.658943240626,
    .sigma_0 = 1.884145499156,
    .s = 1.113514380624,
    .epsilon_cont = -0.2754617951752,
    .epsilon_sov = -0.6472631695864,
    .epsilon_sos = 0.9425803439318,
    .epsilon_tens = -0.4494901823012,
};

const argsGIModel_t argsGIScreen_bbbar = {
    .mn = 0.4713455847642,
    .ms = 0.6283121820133,
    .mc = 1.810505119204,
    .mb = 5.156014766761,
    .b1 = 0.2481528961448,
    .mu = 0.1236411839524,
    .c = -0.658943240626,
    .sigma_0 = 1.884145499156,
    .s = 1.113514380624,
    .epsilon_cont = -0.4995058655201,
    .epsilon_sov = -0.8495903230387,
    .epsilon_sos = -0.498981474843,
    .epsilon_tens = -0.1856496411742,
};

const argsGIModel_t argsGIScreen_ccbar = {
    .mn = 0.4713455847642,
    .ms = 0.6283121820133,
    .mc = 1.810505119204,
    .mb = 5.156014766761,
    .b1 = 0.2575467075473,
    .mu = 0.1453562021339,
    .c = -0.658943240626,
    .sigma_0 = 1.884145499156,
    .s = 1.113514380624,
    .epsilon_cont = -0.32452949845,
    .epsilon_sov = -0.5404734834836,
    .epsilon_sos = 0.9999999508829,
    .epsilon_tens = -0.4999502878773,
};

const argsGIModel_t argsGIScreen_light = {
    .mn = 0.1699938791233,
    .ms = 0.3590841512991,
    .mc = 1.6280000000000,
    .mb = 4.9770000000000,
    .b1 = 0.199309678258,
    .mu = 0.157365986975,
    .c = -0.0845975807305,
    .sigma_0 = 2.899140501817,
    .s = 2.922503884794,
    .epsilon_cont = -0.1031116447922,
    .epsilon_sov = -0.5013015725126,
    .epsilon_sos = 0.3188288566508,
    .epsilon_tens = 0.3714409488222
};

const argsGIModel_t argsGIQuadra_meson = {
    .mn = 0.2177513490091,
    .ms = 0.4643758833827,
    .mc = 1.714570728467,
    .mb = 5.065203991177,
    .b1 = 0.2465976978894,
    .b2 = 1.198634858368e-07,
    .mu = 0.1302357633854,
    .c = -0.472129834672,
    .sigma_0 = 1.402122058075,
    .s = 1.343892159583,
    .epsilon_cont = -0.2563034069373,
    .epsilon_sov = -0.2401174575969,
    .epsilon_sos = 0.9992485642149,
    .epsilon_tens = -0.4928359200016
};

const argsGIModel_t argsGIQuadra_bbbar = {
    .mn = 0.2177513490091,
    .ms = 0.4643758833827,
    .mc = 1.714570728467,
    .mb = 5.065203991177,
    .b1 = 0.2407181813395,
    .b2 = 2.853350333787e-11,
    .mu = 0.1159569166466,
    .c = -0.472129834672,
    .sigma_0 = 1.402122058075,
    .s = 1.343892159583,
    .epsilon_cont = -0.3848890250076,
    .epsilon_sov = -0.1200714395148,
    .epsilon_sos = -0.4977643882449,
    .epsilon_tens = -0.05506924478218
};

const argsGIModel_t argsGIQuadra_ccbar = {
    .mn = 0.2177513490091,
    .ms = 0.4643758833827,
    .mc = 1.714570728467,
    .mb = 5.065203991177,
    .b1 = 0.2274304999277,
    .b2 = 0.01059832207369,
    .mu = 0.1336694609948,
    .c = -0.472129834672,
    .sigma_0 = 1.402122058075,
    .s = 1.343892159583,
    .epsilon_cont = -0.2035767996312,
    .epsilon_sov = -0.249087491019,
    .epsilon_sos = 0.9999995908872,
    .epsilon_tens = -0.5002223930929
};

const argsGIModel_t argsGIQuadra_light = {
    .mn = 0.1976715749049,
    .ms = 0.3807811992895,
    .mc = 1.6280000000000,
    .mb = 4.9770000000000,
    .b1 = 0.1663064200075,
    .b2 = 0.008609382200722,
    .mu = 0.1239587117965,
    .c = -0.1055333840203,
    .sigma_0 = 2.785790175651,
    .s = 1.363344797476,
    .epsilon_cont = -0.1192800474937,
    .epsilon_sov = -0.4983175742134,
    .epsilon_sos = 0.2206097555963,
    .epsilon_tens = -0.1877007219321
};
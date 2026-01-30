/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/basis/color.h>
#include <gemstore/basis/intrin.h>

#include <math.h>
#include <string.h>

static void *CM8(intrin_wfn_t *CM8)
{
    for (int i = 0; i < 8; i++) {
        CM8[i] = intrin_wfn_init(2);
    }

    /* CM8[0] */
    intrin_wfn_push(&CM8[0], 1.0, "bR");

    /* CM8[1] */
    intrin_wfn_push(&CM8[1], 1.0, "bG");

    /* CM8[2] */
    intrin_wfn_push(&CM8[2], -1.0, "gR");

    /* CM8[3] */
    intrin_wfn_push(&CM8[3], 1.0 / sqrt(2.0), "rR");
    intrin_wfn_push(&CM8[3], -1.0 / sqrt(2.0), "gG");

    /* CM8[4] */
    intrin_wfn_push(&CM8[4], 1.0, "rG");

    /* CM8[5] */
    intrin_wfn_push(&CM8[5], 2.0 / sqrt(6.0), "bB");
    intrin_wfn_push(&CM8[5], -1.0 / sqrt(6.0), "rR");
    intrin_wfn_push(&CM8[5], -1.0 / sqrt(6.0), "gG");

    /* CM8[6] */
    intrin_wfn_push(&CM8[6], -1.0, "gB");

    /* CM8[7] */
    intrin_wfn_push(&CM8[7], 1.0, "rB");
}

static void *CB38(intrin_wfn_t *CB38)
{
    for (int i = 0; i < 8; i++) {
        CB38[i] = intrin_wfn_init(3);
    }

    /* CB38[0] */
    intrin_wfn_push(&CB38[0], 1.0 / sqrt(2.0), "rgr");
    intrin_wfn_push(&CB38[0], -1.0 / sqrt(2.0), "grr");

    /* CB38[1] */
    intrin_wfn_push(&CB38[1], 1.0 / sqrt(2.0), "rgg");
    intrin_wfn_push(&CB38[1], -1.0 / sqrt(2.0), "grg");

    /* CB38[2] */
    intrin_wfn_push(&CB38[2], 1.0 / sqrt(2.0), "rbr");
    intrin_wfn_push(&CB38[2], -1.0 / sqrt(2.0), "brr");

    /* CB38[3] */
    intrin_wfn_push(&CB38[3], 1.0 / 2.0, "rbg");
    intrin_wfn_push(&CB38[3], -1.0 / 2.0, "brg");
    intrin_wfn_push(&CB38[3], 1.0 / 2.0, "gbr");
    intrin_wfn_push(&CB38[3], -1.0 / 2.0, "bgr");

    /* CB38[4] */
    intrin_wfn_push(&CB38[4], 1.0 / sqrt(2.0), "gbg");
    intrin_wfn_push(&CB38[4], -1.0 / sqrt(2.0), "bgg");

    /* CB38[5] */
    intrin_wfn_push(&CB38[5], 1.0 / sqrt(3.0), "rgb");
    intrin_wfn_push(&CB38[5], -1.0 / sqrt(3.0), "grb");
    intrin_wfn_push(&CB38[5], -1.0 / sqrt(12.0), "gbr");
    intrin_wfn_push(&CB38[5], 1.0 / sqrt(12.0), "bgr");
    intrin_wfn_push(&CB38[5], -1.0 / sqrt(12.0), "brg");
    intrin_wfn_push(&CB38[5], 1.0 / sqrt(12.0), "rbg");

    /* CB38[6] */
    intrin_wfn_push(&CB38[6], 1.0 / sqrt(2.0), "rbb");
    intrin_wfn_push(&CB38[6], -1.0 / sqrt(2.0), "brb");

    /* CB38[7] */
    intrin_wfn_push(&CB38[7], 1.0 / sqrt(2.0), "gbb");
    intrin_wfn_push(&CB38[7], -1.0 / sqrt(2.0), "bgb");
}

static void *CB68(intrin_wfn_t *CB68)
{
    for (int i = 0; i < 8; i++) {
        CB68[i] = intrin_wfn_init(3);
    }

    /* CB68[0] */
    intrin_wfn_push(&CB68[0], 2.0 / sqrt(6.0), "rrg");
    intrin_wfn_push(&CB68[0], -1.0 / sqrt(6.0), "rgr");
    intrin_wfn_push(&CB68[0], -1.0 / sqrt(6.0), "grr");

    /* CB68[1] */
    intrin_wfn_push(&CB68[1], 1.0 / sqrt(6.0), "rgg");
    intrin_wfn_push(&CB68[1], 1.0 / sqrt(6.0), "grg");
    intrin_wfn_push(&CB68[1], -2.0 / sqrt(6.0), "ggr");

    /* CB68[2] */
    intrin_wfn_push(&CB68[2], 2.0 / sqrt(6.0), "rrb");
    intrin_wfn_push(&CB68[2], -1.0 / sqrt(6.0), "rbr");
    intrin_wfn_push(&CB68[2], -1.0 / sqrt(6.0), "brr");

    /* CB68[3] */
    intrin_wfn_push(&CB68[3], 1.0 / sqrt(3.0), "rgb");
    intrin_wfn_push(&CB68[3], 1.0 / sqrt(3.0), "grb");
    intrin_wfn_push(&CB68[3], -1.0 / sqrt(12.0), "gbr");
    intrin_wfn_push(&CB68[3], -1.0 / sqrt(12.0), "rbg");
    intrin_wfn_push(&CB68[3], -1.0 / sqrt(12.0), "brg");
    intrin_wfn_push(&CB68[3], -1.0 / sqrt(12.0), "bgr");

    /* CB68[4] */
    intrin_wfn_push(&CB68[4], 2.0 / sqrt(6.0), "ggb");
    intrin_wfn_push(&CB68[4], -1.0 / sqrt(6.0), "gbg");
    intrin_wfn_push(&CB68[4], -1.0 / sqrt(6.0), "bgg");

    /* CB68[5] */
    intrin_wfn_push(&CB68[5], 1.0 / 2.0, "rbg");
    intrin_wfn_push(&CB68[5], -1.0 / 2.0, "gbr");
    intrin_wfn_push(&CB68[5], 1.0 / 2.0, "brg");
    intrin_wfn_push(&CB68[5], -1.0 / 2.0, "bgr");

    /* CB68[6] */
    intrin_wfn_push(&CB68[6], 1.0 / sqrt(6.0), "rbb");
    intrin_wfn_push(&CB68[6], 1.0 / sqrt(6.0), "brb");
    intrin_wfn_push(&CB68[6], -2.0 / sqrt(6.0), "bbr");

    /* CB68[7] */
    intrin_wfn_push(&CB68[7], 1.0 / sqrt(6.0), "gbb");
    intrin_wfn_push(&CB68[7], 1.0 / sqrt(6.0), "bgb");
    intrin_wfn_push(&CB68[7], -2.0 / sqrt(6.0), "bbg");
}

intrin_wfn_t color_wfn_meson()
{
    intrin_wfn_t cwf = intrin_wfn_init(2);

    intrin_wfn_push(&cwf, 1.0 / sqrt(3.0), "rR");
    intrin_wfn_push(&cwf, 1.0 / sqrt(3.0), "gG");
    intrin_wfn_push(&cwf, 1.0 / sqrt(3.0), "bB");

    intrin_wfn_trim(&cwf);
    return cwf;
}

intrin_wfn_t color_wfn_baryon()
{
    intrin_wfn_t cwf = intrin_wfn_init(3);

    intrin_wfn_push(&cwf, 1.0 / sqrt(6.0), "rgb");
    intrin_wfn_push(&cwf, -1.0 / sqrt(6.0), "grb");
    intrin_wfn_push(&cwf, 1.0 / sqrt(6.0), "gbr");
    intrin_wfn_push(&cwf, -1.0 / sqrt(6.0), "bgr");
    intrin_wfn_push(&cwf, 1.0 / sqrt(6.0), "brg");
    intrin_wfn_push(&cwf, -1.0 / sqrt(6.0), "rbg");

    intrin_wfn_trim(&cwf);
    return cwf;
}

intrin_wfn_t color_wfn_tetra1()
{
    intrin_wfn_t cwf = intrin_wfn_init(4);
    intrin_wfn_t part1T = color_wfn_meson();
    intrin_wfn_t part2T = color_wfn_meson();

    intrin_wfn_product(&part1T, &part2T, &cwf, 1.0);
    intrin_wfn_free(&part1T);
    intrin_wfn_free(&part2T);

    intrin_wfn_trim(&cwf);
    return cwf;
}

intrin_wfn_t color_wfn_tetra8()
{
    intrin_wfn_t cwf = intrin_wfn_init(4);
    intrin_wfn_t part1T[8];
    intrin_wfn_t part2T[8];
    CM8(part1T);
    CM8(part2T);

    for (int i = 0; i < 8; i++) {
        intrin_wfn_product(&part1T[i], &part2T[i], &cwf, 1.0 / sqrt(8.0));
    }
    for (int i = 0; i < 8; i++) {
        intrin_wfn_free(&part1T[i]);
        intrin_wfn_free(&part2T[i]);
    }

    intrin_wfn_trim(&cwf);
    return cwf;
}

intrin_wfn_t color_wfn_penta1()
{
    intrin_wfn_t cwf = intrin_wfn_init(5);
    intrin_wfn_t part1P = color_wfn_baryon();
    intrin_wfn_t part2P = color_wfn_meson();

    intrin_wfn_product(&part1P, &part2P, &cwf, 1.0);
    intrin_wfn_free(&part1P);
    intrin_wfn_free(&part2P);

    intrin_wfn_trim(&cwf);
    return cwf;
}

intrin_wfn_t color_wfn_penta38()
{
    intrin_wfn_t cwf = intrin_wfn_init(5);
    intrin_wfn_t part1P[8];
    intrin_wfn_t part2P[8];
    CB38(part1P);
    CM8(part2P);

    for (int i = 0; i < 8; i++) {
        intrin_wfn_product(&part1P[i], &part2P[i], &cwf, 1.0 / sqrt(8.0));
    }
    for (int i = 0; i < 8; i++) {
        intrin_wfn_free(&part1P[i]);
        intrin_wfn_free(&part2P[i]);
    }

    intrin_wfn_trim(&cwf);
    return cwf;
}

intrin_wfn_t color_wfn_penta68()
{
    intrin_wfn_t cwf = intrin_wfn_init(5);
    intrin_wfn_t part1P[8];
    intrin_wfn_t part2P[8];
    CB68(part1P);
    CM8(part2P);

    for (int i = 0; i < 8; i++) {
        intrin_wfn_product(&part1P[i], &part2P[i], &cwf, 1.0 / sqrt(8.0));
    }
    for (int i = 0; i < 8; i++) {
        intrin_wfn_free(&part1P[i]);
        intrin_wfn_free(&part2P[i]);
    }

    intrin_wfn_trim(&cwf);
    return cwf;
}

double OrthogonalColor(const intrin_wfn_t *wfn, const intrin_wfn_t *ref)
{
    double overlap = 0;

    for (int i = 0; i < wfn->num_terms; i++) {
        for (int j = 0; j < ref->num_terms; j++) {
            if (strcmp(wfn->configs[i], ref->configs[j]) == 0) {
                overlap += wfn->coeffs[i] * ref->coeffs[j];
            }
        }
    }

    return overlap;
}
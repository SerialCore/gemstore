/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/basis/color.h>
#include <gemstore/basis/intrin.h>

#include <math.h>

intrin_wfn_t color_wfn_meson()
{
    intrin_wfn_t cwf = intrin_wfn_init(2);

    double s3 = 1.0 / sqrt(3.0);

    intrin_wfn_push(&cwf, s3, "rR");
    intrin_wfn_push(&cwf, s3, "gG");
    intrin_wfn_push(&cwf, s3, "bB");

    return cwf;
}

intrin_wfn_t color_wfn_baryon()
{
    intrin_wfn_t cwf = intrin_wfn_init(3);

    double s6 = 1.0 / sqrt(6.0);

    intrin_wfn_push(&cwf, s6, "rgb");
    intrin_wfn_push(&cwf, -s6, "grb");
    intrin_wfn_push(&cwf, s6, "gbr");
    intrin_wfn_push(&cwf, -s6, "bgr");
    intrin_wfn_push(&cwf, s6, "brg");
    intrin_wfn_push(&cwf, -s6, "rbg");

    return cwf;
}

intrin_wfn_t color_wfn_tetra1()
{
    intrin_wfn_t cwf = intrin_wfn_init(4);

    double s9 = 1.0 / 3.0;

    intrin_wfn_push(&cwf, s9, "rRrR");
    intrin_wfn_push(&cwf, s9, "bBbB");
    intrin_wfn_push(&cwf, s9, "gGgG");
    intrin_wfn_push(&cwf, s9, "rRbB");
    intrin_wfn_push(&cwf, s9, "rRgG");
    intrin_wfn_push(&cwf, s9, "bBrR");
    intrin_wfn_push(&cwf, s9, "bBgG");
    intrin_wfn_push(&cwf, s9, "gGrR");
    intrin_wfn_push(&cwf, s9, "gGbB");

    return cwf;
}

intrin_wfn_t color_wfn_tetra8()
{
    intrin_wfn_t cwf = intrin_wfn_init(4);

    double s8 = 1.0 / sqrt(8.0);    /* 1/(2√2) = 1/√8 */
    double s18 = 1.0 / sqrt(18.0);  /* 1/(3√2) = 1/√18 */
    double s72 = 1.0 / sqrt(72.0);  /* 1/(6√2) = 1/√72 */

    intrin_wfn_push(&cwf, s8, "bRrB");
    intrin_wfn_push(&cwf, s8, "rBbR");
    intrin_wfn_push(&cwf, s8, "gRrG");
    intrin_wfn_push(&cwf, s8, "rGgR");
    intrin_wfn_push(&cwf, s8, "gBbG");
    intrin_wfn_push(&cwf, s8, "bGgB");
    intrin_wfn_push(&cwf, s18, "rRrR");
    intrin_wfn_push(&cwf, s18, "gGgG");
    intrin_wfn_push(&cwf, s18, "bBbB");
    intrin_wfn_push(&cwf, -s72, "rRgG");
    intrin_wfn_push(&cwf, -s72, "gGrR");
    intrin_wfn_push(&cwf, -s72, "bBgG");
    intrin_wfn_push(&cwf, -s72, "gGbB");
    intrin_wfn_push(&cwf, -s72, "bBrR");
    intrin_wfn_push(&cwf, -s72, "rRbB");

    return cwf;
}

intrin_wfn_t color_wfn_tetra3()
{
    intrin_wfn_t cwf = intrin_wfn_init(4);

    double s12 = 1.0 / sqrt(12.0);  /* 1/(2√3) = 1/√12 */

    intrin_wfn_push(&cwf,  s12, "rbBR");
    intrin_wfn_push(&cwf, -s12, "brBR");
    intrin_wfn_push(&cwf, -s12, "grGR");
    intrin_wfn_push(&cwf,  s12, "rgGR");
    intrin_wfn_push(&cwf,  s12, "gbBG");
    intrin_wfn_push(&cwf, -s12, "bgBG");
    intrin_wfn_push(&cwf,  s12, "grRG");
    intrin_wfn_push(&cwf, -s12, "rgRG");
    intrin_wfn_push(&cwf, -s12, "gbGB");
    intrin_wfn_push(&cwf,  s12, "bgGB");
    intrin_wfn_push(&cwf, -s12, "rbRB");
    intrin_wfn_push(&cwf,  s12, "brRB");

    return cwf;
}

intrin_wfn_t color_wfn_tetra6()
{
    intrin_wfn_t cwf = intrin_wfn_init(4);

    double s6 = 1.0 / sqrt(6.0);
    double s24 = 1.0 / sqrt(24.0);  /* 1/(2√6) = 1/√24 */

    intrin_wfn_push(&cwf, s6, "rrRR");
    intrin_wfn_push(&cwf, s6, "ggGG");
    intrin_wfn_push(&cwf, s6, "bbBB");
    intrin_wfn_push(&cwf, s24, "rbBR");
    intrin_wfn_push(&cwf, s24, "brBR");
    intrin_wfn_push(&cwf, s24, "grGR");
    intrin_wfn_push(&cwf, s24, "rgGR");
    intrin_wfn_push(&cwf, s24, "gbBG");
    intrin_wfn_push(&cwf, s24, "bgBG");
    intrin_wfn_push(&cwf, s24, "grRG");
    intrin_wfn_push(&cwf, s24, "rgRG");
    intrin_wfn_push(&cwf, s24, "gbGB");
    intrin_wfn_push(&cwf, s24, "bgGB");
    intrin_wfn_push(&cwf, s24, "rbRB");
    intrin_wfn_push(&cwf, s24, "brRB");

    return cwf;
}
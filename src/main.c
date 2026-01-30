/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/basis/color.h>
#include <gemstore/basis/intrin.h>

#include <stdio.h>

void print_help();

int main(int argc, char **argv)
{
    intrin_wfn_t wf, ref_wf;

    wf = color_wfn_penta38();
    printf("ColorWFPenta38:\n");
    intrin_wfn_print(&wf);

    ref_wf = color_wfn_penta68();
    printf("ColorWFPenta68:\n");
    intrin_wfn_print(&ref_wf);

    double ortho = OrthogonalColor(&wf, &ref_wf);
    printf("Orthogonal degree: %f\n", ortho);

    double normal = OrthogonalColor(&ref_wf, &ref_wf);
    printf("Normalized degree: %f\n", normal);

    intrin_wfn_free(&wf);
    intrin_wfn_free(&ref_wf);

    return 0;
}

void print_help()
{
    printf("This is the help page.\n");
}
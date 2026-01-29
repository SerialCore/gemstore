/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/basis/isospin.h>
#include <gemstore/basis/intrin.h>

#include <stdio.h>

void print_help();

int main(int argc, char **argv)
{
    intrin_wfn_t wf;

    wf = isospin_basis(0.5, 'q');
    printf("IsospinBasis[1/2]:\n");
    intrin_wfn_print(&wf);
    intrin_wfn_free(&wf);

    wf = isospin_wfn_meson(1.0, 0.0);
    printf("IsospinWFMeson[1][0]:\n");
    intrin_wfn_print(&wf);
    intrin_wfn_free(&wf);

    wf = isospin_wfn_baryon(1.0, 1.5, 0.5);
    printf("IsospinWFBaryon[1,3/2][1/2]:\n");
    intrin_wfn_print(&wf);
    intrin_wfn_free(&wf);

    wf = isospin_wfn_tetra(1.0, 1.0, 1.0, 0.0);
    printf("IsospinWFTetra[1,1,1][0]:\n");
    intrin_wfn_print(&wf);
    intrin_wfn_free(&wf);

    wf = isospin_wfn_penta(1.0, 0.5, 1.0, 1.5, 0.5);
    printf("IsospinWFPenta[1,1/2,1,3/2][1/2]:\n");
    intrin_wfn_print(&wf);
    intrin_wfn_free(&wf);

    wf = isospin_wfn_hexa(1.0, 0.5, 1.0, 1.5, 1.0, 0.0);
    printf("IsospinWFHexa[1,1/2,1,3/2,1][0]:\n");
    intrin_wfn_print(&wf);
    intrin_wfn_free(&wf);

    return 0;
}

void print_help()
{
    printf("This is the help page.\n");
}
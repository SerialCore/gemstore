/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/entry.h>
#include <gemstore/debug.h>

#include <gemstore/numerical/matrix.h>
#include <gemstore/numerical/eigen.h>

#include <stdio.h>

void print_help();
void print_help()
{
    printf("This is the help page.\n");
}

int main(int argc, char **argv)
{
    //debug_soc_operator();
    //debug_matrix_element();
    //debug_eigen_system();
    //call_spectra_meson_NRScreen(3, 3, 1, 1, 0, 20, 10.0, 0.1);
    call_spectra_meson_GIScreen(3, 3, 0, 0, 0, 20, 10.0, 0.1);
    //call_minuit2_chi2();
}
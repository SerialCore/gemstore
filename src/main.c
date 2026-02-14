/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/entry.h>
#include <gemstore/debug.h>

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
    call_spectra_meson_NRScreen(3, 3, 0, 0, 0, 20, 10.0, 0.1);

    return 0;
}
/*
 * Copyright (C) 2025 Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <cg.h>

#include <stdio.h>

void print_help();

int main(int argc, char **argv)
{
    //double cg = clebsch_gordan(1, 1, 1, -1, 1, 0);
    double cg = cg_3J_symbol(1, 1, 1, -1, 1, 0);
    //double cg = cg_6J_symbol(1, 1, 1, 1, 2, 1);
    //double cg = cg_9J_symbol(1, 1, 2, 2, 1, 1, 2, 2, 3);
    printf("result: %lf", cg);

    return 0;
}

void print_help()
{
    printf("This is the help page.\n");
}
/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/basis/orbit.h>

#include <stdio.h>

void print_help();

int main(int argc, char **argv)
{
    double result = NormalizedSp();
	printf("%f\n", result);

    return 0;
}

void print_help()
{
    printf("This is the help page.\n");
}
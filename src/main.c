/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/entry.h>
#include <gemstore/debug.h>

#include <stdio.h>

void print_help();

int main(int argc, char **argv)
{
    debug_orbit_wfn();

    return 0;
}

void print_help()
{
    printf("This is the help page.\n");
}
/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/entry.h>

#include <stdio.h>
#include <string.h>
#include <getopt.h>

void print_help()
{
    printf("gemstore: hadron spectroscopy tools using Gaussian Expanding Method, Godfrey-Isgur models and more.\n\n");
    printf("Usage: gemstore [--input FILE] [--fitting TARGET] [--print ITEM] [--debug UNIT]\n\n");
    printf("Arguments:\n");
    printf("  -i, --input           input FILE that constains full instructions\n");
    printf("  -f, --fitting         fit TARGET such as GIScreen_meson, GIScreen_ccbar, GIScreen_bbbar, GIQuadra_light\n");
    printf("  -d, --debug           debug UNIT such as su3_product, soc_operator, casimir_operator, \n");
    printf("                        color_wfn, spin_wfn, isospin_wfn, orbit_wfn, eigen_system\n");
    printf("  -p, --print           print ITEM such as potential, wavefunction\n");
    printf("  -h,--help             show this help\n");
    printf("  -v,--version          show version\n\n");
}

int main(int argc, char **argv)
{
    if (argc == 1) {
        print_help();
        return 0;
    }

    int option;
    static struct option long_options[] = {
        {"input",   required_argument, 0, 'i'},
        {"fitting", required_argument, 0, 'f'},
        {"debug",   required_argument, 0, 'd'},
        {"print",   required_argument, 0, 'p'},
        {"help",    no_argument, 0, 'h'},
        {"version", no_argument, 0, 'v'},
        {0, 0, 0, 0}
    };
    while ((option = getopt_long(argc, argv, "i:f:d:p:hv", long_options, NULL)) != -1) {
        switch (option) {
            case 'i':
                entry_compute(optarg);
                break;
            case 'f':
                entry_fitting(optarg);
                break;
            case 'd':
                entry_debug(optarg);
                break;
            case 'p':
                entry_print(optarg);
                break;
            case 'h':
                print_help();
                return 0;
            case 'v':
                printf("gemstore version 0.1.3\n");
                return 0;
        }
    }

    return 0;
}
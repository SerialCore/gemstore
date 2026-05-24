/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/entry.h>
#include <gemstore/print.h>

#include <stdio.h>
#include <string.h>
#include <getopt.h>

int main(int argc, char **argv)
{
    print_logo();
    if (argc == 1) {
        print_help();
        return 0;
    }

    int option;
    static struct option long_options[] = {
        {"compute", required_argument, 0, 'c'},
        {"fitting", required_argument, 0, 'f'},
        {"debug",   required_argument, 0, 'd'},
        {"help",    no_argument, 0, 'h'},
        {"version", no_argument, 0, 'v'},
        {0, 0, 0, 0}
    };
    while ((option = getopt_long(argc, argv, ":c:f:d:hv", long_options, NULL)) != -1) {
        switch (option) {
            case 'c':
                entry_compute(optarg);
                break;
            case 'f':
                entry_fitting(optarg);
                break;
            case 'd':
                entry_debug(optarg);
                break;
            case 'h':
                print_help();
                return 0;
            case 'v':
                printf("gemstore version 1.3\n");
                return 0;
            default:
                print_help();
                return 0;
        }
    }

    print_copyright();

    return 0;
}
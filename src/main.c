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
                printf("gemstore version 0.1.5\n");
                return 0;
        }
    }

    print_copyright();

    return 0;
}
/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/fileio.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int fileio_write_spectra(const char *name, const double *mass, const double *rmsradius,
                         const double *eigenvectors, int nmax, int dim)
{
    FILE *pf;
    int state = 0;

    char path[256];
    sprintf(path, "%s%s", name, ".out");
    pf = fopen(path, "w");

    fprintf(pf, "%s\t%s\t%s\t\t%s\n", "Radial", "Mass", "RMS Radius", "Eigenvectors");
    if (pf != NULL) {
        for (int n = 0; n < nmax; n++) {
            fprintf(pf, "%d\t%10.6f\t%10.6f\t\t", n, mass[n], rmsradius[n]);
            for (int d = 0; d < dim; d++) {
                fprintf(pf, "%10.10f ", eigenvectors[n * dim + d]);
            }
            fprintf(pf, "\n");
        }
        state = fclose(pf) == 0 ? 1 : 0;
    }

    return state;
}

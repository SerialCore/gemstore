/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/fileio.h>
#include <gemstore/math/matrix.h>
#include <gemstore/param/argset.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int write_meson_spectra(const argsInput_t *input, const array_t *mass, const array_t *radius, const matrix_t *vector, int len)
{
    int nmax = input->nmax;

    FILE *pf;
    int state = 0;

    char path[260];
    sprintf(path, "%s%s", input->project, ".out");
    pf = fopen(path, "w");

    fprintf(pf, "%s\t%s\t%s\t%s\n", "Radial", "Mass", "RMSRadius", "Eigenvectors");
    if (pf != NULL) {
        for (int n = 0; n < len; n++) {
            fprintf(pf, "%d\t%10.6f\t%10.6f\t", n + 1, mass->value[n], radius->value[n]);
            for (int d = 0; d < nmax; d++) {
                fprintf(pf, "%10.10f ", vector->value[n][d]);
            }
            fprintf(pf, "\n");
        }
        state = fclose(pf) == 0 ? 1 : 0;
    }

    return state;
}

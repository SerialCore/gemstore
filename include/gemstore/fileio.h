/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef GEMSTORE_FILEIO
#define GEMSTORE_FILEIO

int fileio_write_spectra(const char *path, const double *mass, const double *rmsradius,
                         const double *eigenvectors, int nmax, int dim);

#endif

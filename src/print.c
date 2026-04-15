/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/print.h>
#include <gemstore/math/matrix.h>
#include <gemstore/param/argset.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void print_logo()
{
    printf("\n");
    printf("  ██████╗  ███████╗  ███╗   ███╗  ███████╗  ████████╗   ██████╗  ██████╗  ███████╗\n");
    printf(" ██╔════╝  ██╔════╝  ████╗ ████║  ██╔════╝  ╚══██╔══╝  ██╔═══██║ ██╔══██╗ ██╔════╝\n");
    printf(" ██║  ███╗ █████╗    ██╔████╔██║  ███████╗     ██║     ██║   ██║ ██████╔╝ █████╗  \n");
    printf(" ██║   ██║ ██╔══╝    ██║╚██╔╝██║    ╚══██║     ██║     ██║   ██║ ██╔══██╗ ██╔══╝  \n");
    printf(" ╚██████╔╝ ███████╗  ██║ ╚═╝ ██║  ███████║     ██║     ╚██████╔╝ ██║  ██║ ███████╗\n");
    printf("  ╚═════╝  ╚══════╝  ╚═╝     ╚═╝  ╚══════╝     ╚═╝      ╚═════╝  ╚═╝  ╚═╝ ╚══════╝\n");
    printf("\n");
}

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

void print_copyright()
{
    printf("\n");
    printf("================================================================================\n");
    printf("                              COPYRIGHT NOTICE                                \n");
    printf("================================================================================\n");
    printf("\n");
    printf("GEMSTORE - Hadron Spectroscopy Simulation Tools\n");
    printf("Version 1.0\n");
    printf("\n");
    printf("Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>\n");
    printf("\n");
    printf("This program is free software: you can redistribute it and/or modify\n");
    printf("it under the terms of the GNU General Public License as published by\n");
    printf("the Free Software Foundation, either version 3 of the License, or\n");
    printf("(at your option) any later version.\n");
    printf("\n");
    printf("This program is distributed in the hope that it will be useful,\n");
    printf("but WITHOUT ANY WARRANTY; without even the implied warranty of\n");
    printf("MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the\n");
    printf("GNU General Public License for more details.\n");
    printf("\n");
    printf("You should have received a copy of the GNU General Public License\n");
    printf("along with this program. If not, see <https://www.gnu.org/licenses/>.\n");
    printf("\n");
    printf("SPDX-License-Identifier: GPL-3.0-or-later\n");
    printf("\n");
    printf("================================================================================\n");
    printf("\n");
}

void print_input_parameters(const argsInput_t *input)
{
    if (input == NULL) {
        return;
    }

    printf("\n");
    printf("================================================================================\n");
    printf("                          INPUT PARAMETERS SUMMARY                            \n");
    printf("================================================================================\n");
    printf("\n");

    /* Task and Model Configuration */
    printf("Configuration Settings:\n");
    printf("  Project Name:         %-50s\n", input->project);
    printf("  Task Type:            %-50d\n", input->task);
    printf("  Model Type:           %-50d\n", input->model);
    printf("  System Type:          %-50d\n", input->system);
    printf("\n");

    /* Quark Flavor Configuration */
    printf("Quark Flavor Configuration:\n");
    printf("  Flavor 1 (f1):        %-50d\n", input->f1);
    printf("  Flavor 2 (f2):        %-50d\n", input->f2);
    printf("  Flavor 3 (f3):        %-50d\n", input->f3);
    printf("  Flavor 4 (f4):        %-50d\n", input->f4);
    printf("\n");

    /* Angular Momentum Quantum Numbers */
    printf("Angular Momentum Quantum Numbers:\n");
    printf("  Spin Momentum (S):    %-50.6f\n", input->S);
    printf("  Orbital Momentum (L): %-50.6f\n", input->L);
    printf("  jl (Orbital jl):      %-50.6f\n", input->jl);
    printf("  Total Momentum (J):   %-50.6f\n", input->J);
    printf("\n");

    /* Gaussian Basis Parameters */
    printf("Gaussian Basis Parameters:\n");
    printf("  Number of Gaussians (nmax): %-46d\n", input->nmax);
    printf("  Minimum Range (rmin):       %-46.6f fm\n", input->rmin);
    printf("  Maximum Range (rmax):       %-46.6f fm\n", input->rmax);
    printf("\n");

    /* Model Parameters */
    printf("Model Parameters:\n");
    printf("  Quark Masses:\n");
    printf("    Light Quark (n)   (mn):  %-46.6f GeV\n", input->params.mn);
    printf("    Strange Quark (s) (ms):  %-46.6f GeV\n", input->params.ms);
    printf("    Charm Quark (c)   (mc):  %-46.6f GeV\n", input->params.mc);
    printf("    Bottom Quark (b)  (mb):  %-46.6f GeV\n", input->params.mb);
    printf("\n");

    printf("  Potential Parameters:\n");
    printf("    String Tension (b1):     %-46.6f\n", input->params.b1);
    printf("    Surface Tension (b2):    %-46.6f\n", input->params.b2);
    printf("    Screen Length (mu):      %-46.6f\n", input->params.mu);
    printf("    Constant Potential (c):  %-46.6f\n", input->params.c);
    printf("\n");

    printf("  Gaunov-Isgur Smearing Parameters:\n");
    printf("    sigma_0 (Center):        %-46.6f\n", input->params.sigma_0);
    printf("    sigma (Spin-Spin):       %-46.6f\n", input->params.s);
    printf("    epsilon_cont (Contact):  %-46.6f\n", input->params.epsilon_cont);
    printf("    epsilon_sov (Spin-Orbit):%-46.6f\n", input->params.epsilon_sov);
    printf("    epsilon_sos (Thomas):    %-46.6f\n", input->params.epsilon_sos);
    printf("    epsilon_tens (Tensor):   %-46.6f\n", input->params.epsilon_tens);
    printf("\n");
}

void print_debug_results(const array_t *eigenvalue, const array_t *rmsradius, const matrix_t *eigenvector, int nmax)
{
    if (eigenvalue == NULL || rmsradius == NULL || eigenvector == NULL) {
        return;
    }

    printf("\n");
    printf("================================================================================\n");
    printf("                            DEBUG RESULTS SUMMARY                             \n");
    printf("================================================================================\n");
    printf("\n");

    for (int n = 0; n < nmax; n++) {
        double norm = 0.0;
        double maxc = 0.0;
        for (int i = 0; i < nmax; i++) {
            double c = fabs(eigenvector->value[n][i]);
            norm += c * c;
            if (c > maxc) maxc = c;
        }
        printf("State %2d:  mass=%9.6f  RMS=%6.3f  max|c|=%6.3f  ||c||^2=%13.10f\n", 
            n+1, eigenvalue->value[n], rmsradius->value[n], maxc, norm);
    }

    printf("\n");
}

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

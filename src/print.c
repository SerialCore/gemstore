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
#include <time.h>

void print_logo()
{
    printf("\n");
    printf("  ██████╗  ███████╗  ███╗   ███╗  ███████╗  ████████╗   ██████╗  ██████╗  ███████╗\n");
    printf(" ██╔════╝  ██╔════╝  ████╗ ████║  ██╔════╝  ╚══██╔══╝  ██╔═══██║ ██╔══██╗ ██╔════╝\n");
    printf(" ██║  ███╗ █████╗    ██╔████╔██║  ███████╗     ██║     ██║   ██║ ██████╔╝ █████╗  \n");
    printf(" ██║   ██║ ██╔══╝    ██║╚██╔╝██║  ╚════██║     ██║     ██║   ██║ ██╔══██╗ ██╔══╝  \n");
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

    /* Compute global statistics */
    double min_mass = eigenvalue->value[0];
    double max_mass = eigenvalue->value[0];
    double min_rms = rmsradius->value[0];
    double max_rms = rmsradius->value[0];
    double total_norm = 0.0;
    int num_states = 0;

    for (int n = 0; n < nmax; n++) {
        if (n < eigenvalue->len) {
            if (eigenvalue->value[n] < min_mass) min_mass = eigenvalue->value[n];
            if (eigenvalue->value[n] > max_mass) max_mass = eigenvalue->value[n];
        }
        if (n < rmsradius->len) {
            if (rmsradius->value[n] < min_rms) min_rms = rmsradius->value[n];
            if (rmsradius->value[n] > max_rms) max_rms = rmsradius->value[n];
        }
        num_states++;
    }

    /* Print header row */
    printf("%-6s%-12s%-12s%-12s%-12s%-15s\n", 
        "State", "Mass(GeV)", "RMS(fm)", "max|c|", "||c||^2", "norm_check");
    printf("------+----------+----------+----------+----------+---------------\n");

    /* Print each state's information */
    for (int n = 0; n < nmax; n++) {
        double norm = 0.0;
        double maxc = 0.0;
        
        for (int i = 0; i < nmax; i++) {
            double c = fabs(eigenvector->value[n][i]);
            norm += c * c;
            if (c > maxc) maxc = c;
        }
        
        total_norm += norm;
        
        /* Determine normalization status */
        const char *norm_status = "";
        if (fabs(norm - 1.0) < 1e-6) {
            norm_status = "✓ OK";
        } else if (norm > 1.0) {
            norm_status = "⚠ OVER";
        } else {
            norm_status = "⚠ UNDER";
        }
        
        printf("%-6d%-12.6f%-12.6f%-12.6f%-12.10f%-15s\n", 
            n+1, eigenvalue->value[n], rmsradius->value[n], maxc, norm, norm_status);
    }

    printf("\n");
    printf("GLOBAL STATISTICS:\n");
    printf("  Number of states:    %d\n", num_states);
    printf("  Mass range:          %.6f - %.6f GeV (Δ=%.6f GeV)\n", 
        min_mass, max_mass, max_mass - min_mass);
    printf("  RMS radius range:    %.6f - %.6f fm (Δ=%.6f fm)\n", 
        min_rms, max_rms, max_rms - min_rms);
    printf("  Total norm sum:      %.10f (should ≈ %d)\n", total_norm, nmax);
    printf("  Average norm per state: %.10f\n", total_norm / num_states);
    printf("\n");
}

int write_meson_spectra(const argsInput_t *input, const array_t *mass, const array_t *radius, const matrix_t *vector, int len)
{
    if (input == NULL || mass == NULL || radius == NULL || vector == NULL || len <= 0) {
        fprintf(stderr, "Error: Invalid input parameters to write_meson_spectra()\n");
        return 0;
    }

    int nmax = input->nmax;
    FILE *pf;
    int state = 0;

    char path[260];
    sprintf(path, "%s%s", input->project, ".out");
    pf = fopen(path, "w");

    if (pf == NULL) {
        fprintf(stderr, "Error: Cannot open file %s for writing\n", path);
        return 0;
    }

    /* ========== HEADER SECTION ========== */
    fprintf(pf, "================================================================================\n");
    fprintf(pf, "                    HADRON SPECTROSCOPY RESULTS SUMMARY                       \n");
    fprintf(pf, "================================================================================\n\n");

    /* Timestamp and metadata */
    time_t now = time(NULL);
    struct tm *timeinfo = localtime(&now);
    char time_str[80];
    strftime(time_str, sizeof(time_str), "%Y-%m-%d %H:%M:%S", timeinfo);
    fprintf(pf, "Generated:   %s\n", time_str);
    fprintf(pf, "Project:     %s\n\n", input->project);

    /* ========== INPUT CONFIGURATION ========== */
    fprintf(pf, "INPUT CONFIGURATION:\n");
    fprintf(pf, "  Task:                %d\n", input->task);
    fprintf(pf, "  Model:               %d\n", input->model);
    fprintf(pf, "  System:              %d\n", input->system);
    fprintf(pf, "  Quark Flavors:       f1=%d  f2=%d  f3=%d  f4=%d\n", 
        input->f1, input->f2, input->f3, input->f4);
    fprintf(pf, "  Angular Momentum:    S=%.1f  L=%.1f  J=%.1f  jl=%.1f\n",
        input->S, input->L, input->J, input->jl);
    fprintf(pf, "  Gaussian Basis:      nmax=%d  rmin=%.6f fm  rmax=%.6f fm\n",
        input->nmax, input->rmin, input->rmax);
    fprintf(pf, "\n");

    /* ========== MODEL PARAMETERS ========== */
    fprintf(pf, "MODEL PARAMETERS:\n");
    fprintf(pf, "  Quark Masses:        mn=%.6f  ms=%.6f  mc=%.6f  mb=%.6f GeV\n",
        input->params.mn, input->params.ms, input->params.mc, input->params.mb);
    fprintf(pf, "  Potential:           b1=%.6f  b2=%.6f  mu=%.6f  c=%.6f\n",
        input->params.b1, input->params.b2, input->params.mu, input->params.c);
    fprintf(pf, "  GI Parameters:       sigma_0=%.6f  s=%.6f\n",
        input->params.sigma_0, input->params.s);
    fprintf(pf, "  Smearing:            epsilon_cont=%.6f  epsilon_sov=%.6f  epsilon_sos=%.6f  epsilon_tens=%.6f\n",
        input->params.epsilon_cont, input->params.epsilon_sov, input->params.epsilon_sos, input->params.epsilon_tens);
    fprintf(pf, "\n");

    /* ========== COMPUTE STATISTICS ========== */
    double min_mass = mass->value[0];
    double max_mass = mass->value[0];
    double sum_mass = 0.0;
    double min_radius = radius->value[0];
    double max_radius = radius->value[0];
    double sum_radius = 0.0;

    for (int n = 0; n < len; n++) {
        sum_mass += mass->value[n];
        sum_radius += radius->value[n];
        if (mass->value[n] < min_mass) min_mass = mass->value[n];
        if (mass->value[n] > max_mass) max_mass = mass->value[n];
        if (radius->value[n] < min_radius) min_radius = radius->value[n];
        if (radius->value[n] > max_radius) max_radius = radius->value[n];
    }

    double mean_mass = sum_mass / len;
    double mean_radius = sum_radius / len;

    /* Calculate standard deviation */
    double var_mass = 0.0;
    double var_radius = 0.0;
    for (int n = 0; n < len; n++) {
        var_mass += (mass->value[n] - mean_mass) * (mass->value[n] - mean_mass);
        var_radius += (radius->value[n] - mean_radius) * (radius->value[n] - mean_radius);
    }
    double std_mass = sqrt(var_mass / len);
    double std_radius = sqrt(var_radius / len);

    /* ========== DETAILED SPECTRAL DATA ========== */
    fprintf(pf, "SPECTRAL DATA:\n");
    fprintf(pf, "%-6s%-12s%-12s%-12s%-15s%-15s%-15s\n",
        "State", "Mass(GeV)", "RMS(fm)", "Δmass", "max|coeff|", "||coeff||^2", "norm_stat");
    fprintf(pf, "------+----------+----------+----------+---------------+---------------+---------------\n");

    for (int n = 0; n < len; n++) {
        double norm = 0.0;
        double maxc = 0.0;

        for (int i = 0; i < nmax; i++) {
            double c = fabs(vector->value[n][i]);
            norm += c * c;
            if (c > maxc) maxc = c;
        }

        /* Normalization status */
        const char *norm_stat = "";
        if (fabs(norm - 1.0) < 1e-6) {
            norm_stat = "✓ OK";
        } else if (norm > 1.0) {
            norm_stat = "⚠ OVER";
        } else {
            norm_stat = "⚠ UNDER";
        }

        /* Mass difference from mean */
        double delta_mass = mass->value[n] - mean_mass;

        fprintf(pf, "%-6d%-12.6f%-12.6f%-12.6f%-15.8f%-15.10f%-15s\n",
            n+1, mass->value[n], radius->value[n], delta_mass, maxc, norm, norm_stat);
    }
    fprintf(pf, "\n");

    /* ========== EIGENVECTOR COMPONENTS ========== */
    fprintf(pf, "EIGENVECTOR COMPONENTS:\n");
    fprintf(pf, "(Each row represents one state; columns are eigenvector coefficients)\n\n");

    /* Print column headers for eigenvector components */
    fprintf(pf, "%-6s", "State");
    for (int i = 0; i < nmax; i++) {
        fprintf(pf, "%-14s", "");
        fprintf(pf, "c[%d]", i);
    }
    fprintf(pf, "\n");

    /* Print component separators */
    fprintf(pf, "------");
    for (int i = 0; i < nmax; i++) {
        fprintf(pf, "+----------+----");
    }
    fprintf(pf, "\n");

    /* Print eigenvector components */
    for (int n = 0; n < len; n++) {
        fprintf(pf, "%-6d", n+1);
        for (int d = 0; d < nmax; d++) {
            fprintf(pf, "%-14.10f ", vector->value[n][d]);
        }
        fprintf(pf, "\n");
    }
    fprintf(pf, "\n");

    /* ========== STATISTICAL SUMMARY ========== */
    fprintf(pf, "STATISTICAL SUMMARY:\n");
    fprintf(pf, "  Number of states:    %d\n", len);
    fprintf(pf, "  Gaussian basis size: %d\n", nmax);
    fprintf(pf, "\n");
    fprintf(pf, "  Mass Statistics (GeV):\n");
    fprintf(pf, "    Min:               %.6f\n", min_mass);
    fprintf(pf, "    Max:               %.6f\n", max_mass);
    fprintf(pf, "    Mean:              %.6f\n", mean_mass);
    fprintf(pf, "    Std Dev:           %.6f\n", std_mass);
    fprintf(pf, "    Range:             %.6f\n", max_mass - min_mass);
    fprintf(pf, "\n");
    fprintf(pf, "  RMS Radius Statistics (fm):\n");
    fprintf(pf, "    Min:               %.6f\n", min_radius);
    fprintf(pf, "    Max:               %.6f\n", max_radius);
    fprintf(pf, "    Mean:              %.6f\n", mean_radius);
    fprintf(pf, "    Std Dev:           %.6f\n", std_radius);
    fprintf(pf, "    Range:             %.6f\n", max_radius - min_radius);
    fprintf(pf, "\n");

    /* Eigenvector statistics */
    double min_norm = 1.0;
    double max_norm = 1.0;
    double sum_norm = 0.0;
    double min_maxc = 1.0;
    double max_maxc = 0.0;
    double sum_maxc = 0.0;

    for (int n = 0; n < len; n++) {
        double norm = 0.0;
        double maxc = 0.0;
        for (int i = 0; i < nmax; i++) {
            double c = fabs(vector->value[n][i]);
            norm += c * c;
            if (c > maxc) maxc = c;
        }
        norm = sqrt(norm);
        sum_norm += norm;
        sum_maxc += maxc;
        if (norm < min_norm) min_norm = norm;
        if (norm > max_norm) max_norm = norm;
        if (maxc < min_maxc) min_maxc = maxc;
        if (maxc > max_maxc) max_maxc = maxc;
    }

    fprintf(pf, "  Eigenvector Normalization:\n");
    fprintf(pf, "    Min ||coeff||:     %.10f\n", min_norm);
    fprintf(pf, "    Max ||coeff||:     %.10f\n", max_norm);
    fprintf(pf, "    Mean ||coeff||:    %.10f\n", sum_norm / len);
    fprintf(pf, "\n");
    fprintf(pf, "  Maximum Coefficient Statistics:\n");
    fprintf(pf, "    Min max|c|:        %.10f\n", min_maxc);
    fprintf(pf, "    Max max|c|:        %.10f\n", max_maxc);
    fprintf(pf, "    Mean max|c|:       %.10f\n", sum_maxc / len);
    fprintf(pf, "\n");

    /* ========== FOOTER ========== */
    fprintf(pf, "================================================================================\n");
    fprintf(pf, "End of Hadron Spectroscopy Results\n");
    fprintf(pf, "================================================================================\n");

    state = fclose(pf) == 0 ? 1 : 0;
    if (state == 1) {
        printf("Successfully written hadron spectra to: %s\n", path);
    } else {
        fprintf(stderr, "Error: Failed to close file %s\n", path);
    }

    return state;
}

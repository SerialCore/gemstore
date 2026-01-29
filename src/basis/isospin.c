/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/basis/isospin.h>
#include <gemstore/basis/intrin.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

intrin_wfn_t isospin_basis(double it, const char config)
{
    intrin_wfn_t wf = intrin_wfn_init(1);

    /* q for quark, Q for antiquark */
    if (it == 0.5 && config == 'q') {
        intrin_wfn_push(&wf, 1.0, "u");
    }
    else if (it == -0.5 && config == 'q') {
        intrin_wfn_push(&wf, 1.0, "d");
    }
    else if (it == 0.5 && config == 'Q') {
        intrin_wfn_push(&wf, 1.0, "D");
    }
    else if (it == -0.5 && config == 'Q') {
        intrin_wfn_push(&wf, -1.0, "U");
    }

    return wf;
}

intrin_wfn_t isospin_wfn_meson(double it, double it3)
{
    intrin_wfn_t iwM = intrin_wfn_init(2);
    cg_table_t coupleM = CallCGTable(0.5, 0.5, it, it3);

    for (int M = 0; M < coupleM.num; M++) {
        double ms1 = coupleM.tuples[M].ms1;
        double ms2 = coupleM.tuples[M].ms2;
        double cg = coupleM.tuples[M].cg;
        intrin_wfn_t part1M = isospin_basis(ms1, 'q');
        intrin_wfn_t part2M = isospin_basis(ms2, 'Q');

        double coeff = cg * part1M.coeffs[0] * part2M.coeffs[0];
        char config[3];
        sprintf(config, "%s%s", part1M.configs[0], part2M.configs[0]);
        intrin_wfn_push(&iwM, coeff, config);
        intrin_wfn_free(&part1M);
        intrin_wfn_free(&part2M);
    }
    cg_table_free(&coupleM);

    intrin_wfn_t trimmed = intrin_wfn_trim(&iwM);
    intrin_wfn_free(&iwM);

    return trimmed;
}

intrin_wfn_t isospin_wfn_diquark(double it, double it3)
{
    intrin_wfn_t iwM = intrin_wfn_init(2);
    cg_table_t coupleM = CallCGTable(0.5, 0.5, it, it3);

    for (int M = 0; M < coupleM.num; M++) {
        double ms1 = coupleM.tuples[M].ms1;
        double ms2 = coupleM.tuples[M].ms2;
        double cg = coupleM.tuples[M].cg;
        intrin_wfn_t part1M = isospin_basis(ms1, 'q');
        intrin_wfn_t part2M = isospin_basis(ms2, 'q');

        double coeff = cg * part1M.coeffs[0] * part2M.coeffs[0];
        char config[3];
        sprintf(config, "%s%s", part1M.configs[0], part2M.configs[0]);
        intrin_wfn_push(&iwM, coeff, config);
        intrin_wfn_free(&part1M);
        intrin_wfn_free(&part2M);
    }
    cg_table_free(&coupleM);

    intrin_wfn_t trimmed = intrin_wfn_trim(&iwM);
    intrin_wfn_free(&iwM);

    return trimmed;
}

intrin_wfn_t isospin_wfn_baryon(double i12, double it, double it3)
{
    intrin_wfn_t iwB = intrin_wfn_init(3);

    if (i12 == 0.0) {
        intrin_wfn_t part1B = isospin_wfn_diquark(0.0, 0.0);
        intrin_wfn_t part2B = isospin_basis(it3, 'q');

        for (int p1B = 0; p1B < part1B.num_terms; p1B++) {
            for (int p2B = 0; p2B < part2B.num_terms; p2B++) {
                double coeff = part1B.coeffs[p1B] * part2B.coeffs[p2B];
                char config[4];
                sprintf(config, "%s%s", part1B.configs[p1B], part2B.configs[p2B]);
                intrin_wfn_push(&iwB, coeff, config);
            }
        }
        intrin_wfn_free(&part1B);
        intrin_wfn_free(&part2B);
    }
    if (i12 == 1.0) {
        cg_table_t coupleB = CallCGTable(i12, 0.5, it, it3);

        for (int B = 0; B < coupleB.num; B++) {
            double ms12 = coupleB.tuples[B].ms1;
            double ms3 = coupleB.tuples[B].ms2;
            double cg = coupleB.tuples[B].cg;
            intrin_wfn_t part1B = isospin_wfn_diquark(1.0, ms12);
            intrin_wfn_t part2B = isospin_basis(ms3, 'q');

            for (int p1B = 0; p1B < part1B.num_terms; p1B++) {
                for (int p2B = 0; p2B < part2B.num_terms; p2B++) {
                    double coeff = cg * part1B.coeffs[p1B] * part2B.coeffs[p2B];
                    char config[4];
                    sprintf(config, "%s%s", part1B.configs[p1B], part2B.configs[p2B]);
                    intrin_wfn_push(&iwB, coeff, config);
                }
            }
            intrin_wfn_free(&part1B);
            intrin_wfn_free(&part2B);
        }
        cg_table_free(&coupleB);
    }

    intrin_wfn_t trimmed = intrin_wfn_trim(&iwB);
    intrin_wfn_free(&iwB);

    return trimmed;
}

intrin_wfn_t isospin_wfn_tetra(double i12, double i34, double it, double it3)
{
    intrin_wfn_t iwT = intrin_wfn_init(4);

    if (i12 == 0.0 && i34 == 0.0) {
        intrin_wfn_t part1T = isospin_wfn_meson(0.0, 0.0);
        intrin_wfn_t part2T = isospin_wfn_meson(0.0, 0.0);

        for (int p1T = 0; p1T < part1T.num_terms; p1T++) {
            for (int p2T = 0; p2T < part2T.num_terms; p2T++) {
                double coeff = part1T.coeffs[p1T] * part2T.coeffs[p2T];
                char config[5];
                sprintf(config, "%s%s", part1T.configs[p1T], part2T.configs[p2T]);
                intrin_wfn_push(&iwT, coeff, config);
            }
        }
        intrin_wfn_free(&part1T);
        intrin_wfn_free(&part2T);
    }
    if (i12 == 0.0 && i34 == 1.0) {
        intrin_wfn_t part1T = isospin_wfn_meson(0.0, 0.0);
        intrin_wfn_t part2T = isospin_wfn_meson(1.0, it3);

        for (int p1T = 0; p1T < part1T.num_terms; p1T++) {
            for (int p2T = 0; p2T < part2T.num_terms; p2T++) {
                double coeff = part1T.coeffs[p1T] * part2T.coeffs[p2T];
                char config[5];
                sprintf(config, "%s%s", part1T.configs[p1T], part2T.configs[p2T]);
                intrin_wfn_push(&iwT, coeff, config);
            }
        }
        intrin_wfn_free(&part1T);
        intrin_wfn_free(&part2T);
    }
    if (i12 == 1.0 && i34 == 0.0) {
        intrin_wfn_t part1T = isospin_wfn_meson(1.0, it3);
        intrin_wfn_t part2T = isospin_wfn_meson(0.0, 0.0);

        for (int p1T = 0; p1T < part1T.num_terms; p1T++) {
            for (int p2T = 0; p2T < part2T.num_terms; p2T++) {
                double coeff = part1T.coeffs[p1T] * part2T.coeffs[p2T];
                char config[5];
                sprintf(config, "%s%s", part1T.configs[p1T], part2T.configs[p2T]);
                intrin_wfn_push(&iwT, coeff, config);
            }
        }
        intrin_wfn_free(&part1T);
        intrin_wfn_free(&part2T);
    }
    if (i12 == 1.0 && i34 == 1.0) {
        cg_table_t coupleT = CallCGTable(i12, i34, it, it3);

        for (int T = 0; T < coupleT.num; T++) {
            double ms12 = coupleT.tuples[T].ms1;
            double ms34 = coupleT.tuples[T].ms2;
            double cg = coupleT.tuples[T].cg;
            intrin_wfn_t part1T = isospin_wfn_meson(1.0, ms12);
            intrin_wfn_t part2T = isospin_wfn_meson(1.0, ms34);

            for (int p1T = 0; p1T < part1T.num_terms; p1T++) {
                for (int p2T = 0; p2T < part2T.num_terms; p2T++) {
                    double coeff = cg * part1T.coeffs[p1T] * part2T.coeffs[p2T];
                    char config[5];
                    sprintf(config, "%s%s", part1T.configs[p1T], part2T.configs[p2T]);
                    intrin_wfn_push(&iwT, coeff, config);
                }
            }
            intrin_wfn_free(&part1T);
            intrin_wfn_free(&part2T);
        }
        cg_table_free(&coupleT);
    }

    intrin_wfn_t trimmed = intrin_wfn_trim(&iwT);
    intrin_wfn_free(&iwT);

    return trimmed;
}

intrin_wfn_t isospin_wfn_penta(double i12, double i123, double i45, double it, double it3)
{   
    intrin_wfn_t iwP = intrin_wfn_init(5);

    if (i45 == 0.0) {
        intrin_wfn_t part1P = isospin_wfn_baryon(i12, i123, it3);
        intrin_wfn_t part2P = isospin_wfn_meson(0.0, 0.0);

        for (int p1P = 0; p1P < part1P.num_terms; p1P++) {
            for (int p2P = 0; p2P < part2P.num_terms; p2P++) {
                double coeff = part1P.coeffs[p1P] * part2P.coeffs[p2P];
                char config[6];
                sprintf(config, "%s%s", part1P.configs[p1P], part2P.configs[p2P]);
                intrin_wfn_push(&iwP, coeff, config);
            }
        }
        intrin_wfn_free(&part1P);
        intrin_wfn_free(&part2P);
    }
    if (i45 == 1.0) {
        cg_table_t coupleP = CallCGTable(i123, i45, it, it3);

        for (int P = 0; P < coupleP.num; P++) {
            double ms123 = coupleP.tuples[P].ms1;
            double ms45 = coupleP.tuples[P].ms2;
            double cg = coupleP.tuples[P].cg;
            intrin_wfn_t part1P = isospin_wfn_baryon(i12, i123, ms123);
            intrin_wfn_t part2P = isospin_wfn_meson(1.0, ms45);

            for (int p1P = 0; p1P < part1P.num_terms; p1P++) {
                for (int p2P = 0; p2P < part2P.num_terms; p2P++) {
                    double coeff = cg * part1P.coeffs[p1P] * part2P.coeffs[p2P];
                    char config[6];
                    sprintf(config, "%s%s", part1P.configs[p1P], part2P.configs[p2P]);
                    intrin_wfn_push(&iwP, coeff, config);
                }
            }
            intrin_wfn_free(&part1P);
            intrin_wfn_free(&part2P);
        }
        cg_table_free(&coupleP);
    }

    intrin_wfn_t trimmed = intrin_wfn_trim(&iwP);
    intrin_wfn_free(&iwP);

    return trimmed;
}

intrin_wfn_t isospin_wfn_hexa(double i12, double i123, double i45, double i456, double it, double it3)
{
    intrin_wfn_t iwH = intrin_wfn_init(6);
    cg_table_t coupleH = CallCGTable(i123, i456, it, it3);

    for (int H = 0; H < coupleH.num; H++) {
        double ms123 = coupleH.tuples[H].ms1;
        double ms456 = coupleH.tuples[H].ms2;
        double cg = coupleH.tuples[H].cg;
        intrin_wfn_t part1H = isospin_wfn_baryon(i12, i123, ms123);
        intrin_wfn_t part2H = isospin_wfn_baryon(i45, i456, ms456);

        for (int p1H = 0; p1H < part1H.num_terms; p1H++) {
            for (int p2H = 0; p2H < part2H.num_terms; p2H++) {
                double coeff = cg * part1H.coeffs[p1H] * part2H.coeffs[p2H];
                char config[7];
                sprintf(config, "%s%s", part1H.configs[p1H], part2H.configs[p2H]);
                intrin_wfn_push(&iwH, coeff, config);
            }
        }
        intrin_wfn_free(&part1H);
        intrin_wfn_free(&part2H);
    }
    cg_table_free(&coupleH);

    intrin_wfn_t trimmed = intrin_wfn_trim(&iwH);
    intrin_wfn_free(&iwH);

    return trimmed;
}
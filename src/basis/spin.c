/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/basis/spin.h>
#include <gemstore/basis/soc.h>
#include <gemstore/basis/intrin.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct cg_tuple
{
    double ms1;             /* the third component of s1 */
    double ms2;             /* the third component of s2 */
    double cg;              /* Clebsch-Gordan coefficient */
} cg_tuple_t;

typedef struct cg_table
{
    int num;                /* number of entries */
    cg_tuple_t *tuples;     /* Clebsch-Gordan table */
} cg_table_t;

static void cg_table_free(cg_table_t *t)
{
    free(t->tuples);
    t->num = 0;
    t->tuples = NULL;
}

static cg_table_t CallCGTable(double s1, double s2, double st, double st3)
{
    cg_table_t couple;
    couple.num = 0;
    couple.tuples = NULL;

    for (double ms1 = -s1; ms1 <= s1 + 1e-10; ms1 += 1.0) {
        for (double ms2 = -s2; ms2 <= s2 + 1e-10; ms2 += 1.0) {
            if (ms1 + ms2 - st3 == 0.0) {
                double cg = clebsch_gordan(s1, ms1, s2, ms2, st, st3);
                if (cg != 0.0) {
                    couple.num++;
                    couple.tuples = (cg_tuple_t *)realloc(couple.tuples, couple.num * sizeof(cg_tuple_t));
                    couple.tuples[couple.num - 1].ms1 = ms1;
                    couple.tuples[couple.num - 1].ms2 = ms2;
                    couple.tuples[couple.num - 1].cg = cg;
                }
            }
        }
    }

    return couple;
}

intrin_wfn_t spin_basis(double st)
{
    intrin_wfn_t wf = intrin_wfn_init(1);

    if (st == 0.5) {
        intrin_wfn_push(&wf, 1.0, "1");
    }
    else if (st == -0.5) {
        intrin_wfn_push(&wf, 1.0, "0");
    }

    return wf;
}

intrin_wfn_t spin_wfn_meson(double st, double st3)
{
    intrin_wfn_t swM = intrin_wfn_init(2);
    cg_table_t coupleM = CallCGTable(0.5, 0.5, st, st3);

    for (int M = 0; M < coupleM.num; M++) {
        double ms1 = coupleM.tuples[M].ms1;
        double ms2 = coupleM.tuples[M].ms2;
        double cg = coupleM.tuples[M].cg;
        intrin_wfn_t part1M = spin_basis(ms1);
        intrin_wfn_t part2M = spin_basis(ms2);

        double coeff = cg * part1M.coeffs[0] * part2M.coeffs[0];
        char config[3];
        sprintf(config, "%s%s", part1M.configs[0], part2M.configs[0]);
        intrin_wfn_push(&swM, coeff, config);
        intrin_wfn_free(&part1M);
        intrin_wfn_free(&part2M);
    }
    cg_table_free(&coupleM);

    intrin_wfn_t trimmed = intrin_wfn_trim(&swM);
    intrin_wfn_free(&swM);
    
    return trimmed;
}

intrin_wfn_t spin_wfn_baryon(double s12, double st, double st3)
{
    intrin_wfn_t swB = intrin_wfn_init(3);

    if (s12 == 0.0) {
        intrin_wfn_t part1B = spin_wfn_meson(0.0, 0.0);
        intrin_wfn_t part2B = spin_basis(st3);

        for (int p1B = 0; p1B < part1B.num_terms; p1B++) {
            for (int p2B = 0; p2B < part2B.num_terms; p2B++) {
                double coeff = part1B.coeffs[p1B] * part2B.coeffs[p2B];
                char config[4];
                sprintf(config, "%s%s", part1B.configs[p1B], part2B.configs[p2B]);
                intrin_wfn_push(&swB, coeff, config);
            }
        }
        intrin_wfn_free(&part1B);
        intrin_wfn_free(&part2B);
    }
    if (s12 == 1.0) {
        cg_table_t coupleB = CallCGTable(s12, 0.5, st, st3);

        for (int B = 0; B < coupleB.num; B++) {
            double ms12 = coupleB.tuples[B].ms1;
            double ms3 = coupleB.tuples[B].ms2;
            double cg = coupleB.tuples[B].cg;
            intrin_wfn_t part1B = spin_wfn_meson(1.0, ms12);
            intrin_wfn_t part2B = spin_basis(ms3);

            for (int p1B = 0; p1B < part1B.num_terms; p1B++) {
                for (int p2B = 0; p2B < part2B.num_terms; p2B++) {
                    double coeff = cg * part1B.coeffs[p1B] * part2B.coeffs[p2B];
                    char config[4];
                    sprintf(config, "%s%s", part1B.configs[p1B], part2B.configs[p2B]);
                    intrin_wfn_push(&swB, coeff, config);
                }
            }
            intrin_wfn_free(&part1B);
            intrin_wfn_free(&part2B);
        }
        cg_table_free(&coupleB);
    }

    intrin_wfn_t trimmed = intrin_wfn_trim(&swB);
    intrin_wfn_free(&swB);
    
    return trimmed;
}

intrin_wfn_t spin_wfn_tetra(double s12, double s34, double st, double st3)
{
    intrin_wfn_t swT = intrin_wfn_init(4);

    if (s12 == 0.0 && s34 == 0.0) {
        intrin_wfn_t part1T = spin_wfn_meson(0.0, 0.0);
        intrin_wfn_t part2T = spin_wfn_meson(0.0, 0.0);

        for (int p1T = 0; p1T < part1T.num_terms; p1T++) {
            for (int p2T = 0; p2T < part2T.num_terms; p2T++) {
                double coeff = part1T.coeffs[p1T] * part2T.coeffs[p2T];
                char config[5];
                sprintf(config, "%s%s", part1T.configs[p1T], part2T.configs[p2T]);
                intrin_wfn_push(&swT, coeff, config);
            }
        }
        intrin_wfn_free(&part1T);
        intrin_wfn_free(&part2T);
    }
    if (s12 == 0.0 && s34 == 1.0) {
        intrin_wfn_t part1T = spin_wfn_meson(0.0, 0.0);
        intrin_wfn_t part2T = spin_wfn_meson(1.0, st3);

        for (int p1T = 0; p1T < part1T.num_terms; p1T++) {
            for (int p2T = 0; p2T < part2T.num_terms; p2T++) {
                double coeff = part1T.coeffs[p1T] * part2T.coeffs[p2T];
                char config[5];
                sprintf(config, "%s%s", part1T.configs[p1T], part2T.configs[p2T]);
                intrin_wfn_push(&swT, coeff, config);
            }
        }
        intrin_wfn_free(&part1T);
        intrin_wfn_free(&part2T);
    }
    if (s12 == 1.0 && s34 == 0.0) {
        intrin_wfn_t part1T = spin_wfn_meson(1.0, st3);
        intrin_wfn_t part2T = spin_wfn_meson(0.0, 0.0);

        for (int p1T = 0; p1T < part1T.num_terms; p1T++) {
            for (int p2T = 0; p2T < part2T.num_terms; p2T++) {
                double coeff = part1T.coeffs[p1T] * part2T.coeffs[p2T];
                char config[5];
                sprintf(config, "%s%s", part1T.configs[p1T], part2T.configs[p2T]);
                intrin_wfn_push(&swT, coeff, config);
            }
        }
        intrin_wfn_free(&part1T);
        intrin_wfn_free(&part2T);
    }
    if (s12 == 1.0 && s34 == 1.0) {
        cg_table_t coupleT = CallCGTable(s12, s34, st, st3);

        for (int T = 0; T < coupleT.num; T++) {
            double ms12 = coupleT.tuples[T].ms1;
            double ms34 = coupleT.tuples[T].ms2;
            double cg = coupleT.tuples[T].cg;
            intrin_wfn_t part1T = spin_wfn_meson(1.0, ms12);
            intrin_wfn_t part2T = spin_wfn_meson(1.0, ms34);

            for (int p1T = 0; p1T < part1T.num_terms; p1T++) {
                for (int p2T = 0; p2T < part2T.num_terms; p2T++) {
                    double coeff = cg * part1T.coeffs[p1T] * part2T.coeffs[p2T];
                    char config[5];
                    sprintf(config, "%s%s", part1T.configs[p1T], part2T.configs[p2T]);
                    intrin_wfn_push(&swT, coeff, config);
                }
            }
            intrin_wfn_free(&part1T);
            intrin_wfn_free(&part2T);
        }
        cg_table_free(&coupleT);
    }

    intrin_wfn_t trimmed = intrin_wfn_trim(&swT);
    intrin_wfn_free(&swT);
    
    return trimmed;
}

intrin_wfn_t spin_wfn_penta(double s12, double s123, double s45, double st, double st3)
{   
    intrin_wfn_t swP = intrin_wfn_init(5);

    if (s45 == 0.0) {
        intrin_wfn_t part1P = spin_wfn_baryon(s12, s123, st3);
        intrin_wfn_t part2P = spin_wfn_meson(0.0, 0.0);

        for (int p1P = 0; p1P < part1P.num_terms; p1P++) {
            for (int p2P = 0; p2P < part2P.num_terms; p2P++) {
                double coeff = part1P.coeffs[p1P] * part2P.coeffs[p2P];
                char config[6];
                sprintf(config, "%s%s", part1P.configs[p1P], part2P.configs[p2P]);
                intrin_wfn_push(&swP, coeff, config);
            }
        }
        intrin_wfn_free(&part1P);
        intrin_wfn_free(&part2P);
    }
    if (s45 == 1.0) {
        cg_table_t coupleP = CallCGTable(s123, s45, st, st3);

        for (int P = 0; P < coupleP.num; P++) {
            double ms123 = coupleP.tuples[P].ms1;
            double ms45 = coupleP.tuples[P].ms2;
            double cg = coupleP.tuples[P].cg;
            intrin_wfn_t part1P = spin_wfn_baryon(s12, s123, ms123);
            intrin_wfn_t part2P = spin_wfn_meson(1.0, ms45);

            for (int p1P = 0; p1P < part1P.num_terms; p1P++) {
                for (int p2P = 0; p2P < part2P.num_terms; p2P++) {
                    double coeff = cg * part1P.coeffs[p1P] * part2P.coeffs[p2P];
                    char config[6];
                    sprintf(config, "%s%s", part1P.configs[p1P], part2P.configs[p2P]);
                    intrin_wfn_push(&swP, coeff, config);
                }
            }
            intrin_wfn_free(&part1P);
            intrin_wfn_free(&part2P);
        }
        cg_table_free(&coupleP);
    }

    intrin_wfn_t trimmed = intrin_wfn_trim(&swP);
    intrin_wfn_free(&swP);
    
    return trimmed;
}

intrin_wfn_t spin_wfn_hexa(double s12, double s123, double s45, double s456, double st, double st3)
{
    intrin_wfn_t swH = intrin_wfn_init(6);
    cg_table_t coupleH = CallCGTable(s123, s456, st, st3);

    for (int H = 0; H < coupleH.num; H++) {
        double ms123 = coupleH.tuples[H].ms1;
        double ms456 = coupleH.tuples[H].ms2;
        double cg = coupleH.tuples[H].cg;
        intrin_wfn_t part1H = spin_wfn_baryon(s12, s123, ms123);
        intrin_wfn_t part2H = spin_wfn_baryon(s45, s456, ms456);

        for (int p1H = 0; p1H < part1H.num_terms; p1H++) {
            for (int p2H = 0; p2H < part2H.num_terms; p2H++) {
                double coeff = cg * part1H.coeffs[p1H] * part2H.coeffs[p2H];
                char config[7];
                sprintf(config, "%s%s", part1H.configs[p1H], part2H.configs[p2H]);
                intrin_wfn_push(&swH, coeff, config);
            }
        }
        intrin_wfn_free(&part1H);
        intrin_wfn_free(&part2H);
    }
    cg_table_free(&coupleH);

    intrin_wfn_t trimmed = intrin_wfn_trim(&swH);
    intrin_wfn_free(&swH);
    
    return trimmed;
}
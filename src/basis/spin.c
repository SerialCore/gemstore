/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/basis/spin.h>
#include <gemstore/basis/intrin.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

intrin_wfn_t spin_basis(double st)
{
    intrin_wfn_t wf = intrin_wfn_init(1);

    /* 1 for spin-up, 0 for spin-down */
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

        intrin_wfn_product(&part1M, &part2M, &swM, cg);
        intrin_wfn_free(&part1M);
        intrin_wfn_free(&part2M);
    }
    cg_table_free(&coupleM);

    intrin_wfn_trim(&swM);    
    return swM;
}

intrin_wfn_t spin_wfn_baryon(double s12, double st, double st3)
{
    intrin_wfn_t swB = intrin_wfn_init(3);

    if (s12 == 0.0) {
        intrin_wfn_t part1B = spin_wfn_meson(0.0, 0.0);
        intrin_wfn_t part2B = spin_basis(st3);

        intrin_wfn_product(&part1B, &part2B, &swB, 1.0);
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

            intrin_wfn_product(&part1B, &part2B, &swB, cg);
            intrin_wfn_free(&part1B);
            intrin_wfn_free(&part2B);
        }
        cg_table_free(&coupleB);
    }

    intrin_wfn_trim(&swB);  
    return swB;
}

intrin_wfn_t spin_wfn_tetra(double s12, double s34, double st, double st3)
{
    intrin_wfn_t swT = intrin_wfn_init(4);

    if (s12 == 0.0 && s34 == 0.0) {
        intrin_wfn_t part1T = spin_wfn_meson(0.0, 0.0);
        intrin_wfn_t part2T = spin_wfn_meson(0.0, 0.0);

        intrin_wfn_product(&part1T, &part2T, &swT, 1.0);
        intrin_wfn_free(&part1T);
        intrin_wfn_free(&part2T);
    }
    if (s12 == 0.0 && s34 == 1.0) {
        intrin_wfn_t part1T = spin_wfn_meson(0.0, 0.0);
        intrin_wfn_t part2T = spin_wfn_meson(1.0, st3);

        intrin_wfn_product(&part1T, &part2T, &swT, 1.0);
        intrin_wfn_free(&part1T);
        intrin_wfn_free(&part2T);
    }
    if (s12 == 1.0 && s34 == 0.0) {
        intrin_wfn_t part1T = spin_wfn_meson(1.0, st3);
        intrin_wfn_t part2T = spin_wfn_meson(0.0, 0.0);

        intrin_wfn_product(&part1T, &part2T, &swT, 1.0);
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

            intrin_wfn_product(&part1T, &part2T, &swT, cg);
            intrin_wfn_free(&part1T);
            intrin_wfn_free(&part2T);
        }
        cg_table_free(&coupleT);
    }

    intrin_wfn_trim(&swT);
    return swT;
}

intrin_wfn_t spin_wfn_penta(double s12, double s123, double s45, double st, double st3)
{   
    intrin_wfn_t swP = intrin_wfn_init(5);

    if (s45 == 0.0) {
        intrin_wfn_t part1P = spin_wfn_baryon(s12, s123, st3);
        intrin_wfn_t part2P = spin_wfn_meson(0.0, 0.0);

        intrin_wfn_product(&part1P, &part2P, &swP, 1.0);
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

            intrin_wfn_product(&part1P, &part2P, &swP, cg);
            intrin_wfn_free(&part1P);
            intrin_wfn_free(&part2P);
        }
        cg_table_free(&coupleP);
    }

    intrin_wfn_trim(&swP);
    return swP;
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

        intrin_wfn_product(&part1H, &part2H, &swH, cg);
        intrin_wfn_free(&part1H);
        intrin_wfn_free(&part2H);
    }
    cg_table_free(&coupleH);

    intrin_wfn_trim(&swH);
    return swH;
}
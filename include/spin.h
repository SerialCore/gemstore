/*
 * Copyright (C) 2025 Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef SPIN_H
#define SPIN_H

/* return Clebsch-Gordan coefficient for given angular momenta:
 * <j1m1j2m2|jm>
 */
double clebsch_gordan(double j1, double m1, double j2, double m2, double j, double m);

/* return 3J symbol for given angular momenta:
 * j1, j2, j
 * m1, m2, m
 */
double cg_3J_symbol(double j1, double m1, double j2, double m2, double j, double m);

/* return 6J symbol for given angular momenta:
 * j1, j2, j12
 * j3,  j, j23
 */
double cg_6J_symbol(double j1, double j2, double j12, double j3, double j, double j23);

/* return 9J symbol for given angular momenta:
 *  j1,  j2, j12
 *  j3,  j4, j34
 * j13, j24,   j
*/
double cg_9J_symbol(double j1, double j2, double j12, double j3, double j4, double j34, double j13, double j24, double j);

#endif
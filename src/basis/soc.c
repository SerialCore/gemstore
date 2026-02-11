/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <gemstore/basis/soc.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Precomputed factorials (d/2)! indexed by d */
static const double factorial_2[128] = {
	1.00000000000000E+00, 8.86226925452758E-01, 1.00000000000000E+00, 1.32934038817914E+00,
    2.00000000000000E+00, 3.32335097044784E+00, 6.00000000000000E+00, 1.16317283965675E+01,
    2.40000000000000E+01, 5.23427777845535E+01, 1.20000000000000E+02, 2.87885277815044E+02, 
    7.20000000000000E+02, 1.87125430579779E+03, 5.04000000000000E+03, 1.40344072934834E+04,
    4.03200000000000E+04, 1.19292461994609E+05, 3.62880000000000E+05, 1.13327838894879E+06,
    3.62880000000000E+06, 1.18994230839623E+07, 3.99168000000000E+07, 1.36843365465566E+08,
    4.79001600000000E+08, 1.71054206831957E+09, 6.22702080000000E+09, 2.30923179223142E+10,
    8.71782912000000E+10, 3.34838609873556E+11, 1.30767436800000E+12, 5.18999845304013E+12,
    2.09227898880000E+13, 8.56349744751621E+13, 3.55687428096000E+14, 1.49861205331534E+15,
    6.40237370572800E+15, 2.77243229863337E+16, 1.21645100408832E+17, 5.40624298233507E+17,
    2.43290200817664E+18, 1.10827981137869E+19, 5.10909421717094E+19, 2.38280159446418E+20,
    1.12400072777761E+21, 5.36130358754442E+21, 2.58520167388850E+22, 1.25990634307294E+23,
    6.20448401733239E+23, 3.08677054052870E+24, 1.55112100433310E+25, 7.87126487834818E+25,
    4.03291461126606E+26, 2.08588519276227E+27, 1.08888694504183E+28, 5.73618428009623E+28,
    3.04888344611714E+29, 1.63481251982743E+30, 8.84176199373970E+30, 4.82269693349091E+31,
    2.65252859812191E+32, 1.47092256471473E+33, 8.22283865417792E+33, 4.63340607885139E+34,
    2.63130836933694E+35, 1.50585697562670E+36, 8.68331761881189E+36, 5.04462086834945E+37,
    2.95232799039604E+38, 1.74039419958056E+39, 1.03331479663861E+40, 6.17839940851099E+40,
    3.71993326789901E+41, 2.25511578410651E+42, 1.37637530912264E+43, 8.45668419039942E+43,
    5.23022617466601E+44, 3.25582341330378E+45, 2.03978820811974E+46, 1.28605024825499E+47,
    8.15915283247898E+47, 5.20850350543272E+48, 3.34525266131638E+49, 2.16152895475458E+50,
    1.40500611775288E+51, 9.18649805770695E+51, 6.04152630633738E+52, 3.99612665510252E+53,
    2.65827157478844E+54, 1.77827636152062E+55, 1.19622220865480E+56, 8.09115744491884E+56,
    5.50262215981209E+57, 3.76238821188726E+58, 2.58623241511168E+59, 1.78713440064645E+60,
    1.24139155925361E+61, 8.66760184313527E+61, 6.08281864034267E+62, 4.29046291235196E+63,
    3.04140932017134E+64, 2.16668377073774E+65, 1.55111875328738E+66, 1.11584214192994E+67,
    8.06581751709439E+67, 5.85817124513216E+68, 4.27488328406003E+69, 3.13412161614571E+70,
    2.30843697339241E+71, 1.70809628079941E+72, 1.26964033536583E+73, 9.47993435843673E+73,
    7.10998587804864E+74, 5.35616291251675E+75, 4.05269195048772E+76, 3.07979367469713E+77,
    2.35056133128288E+78, 1.80167929969782E+79, 1.38683118545690E+80, 1.07199918332020E+81,
    8.32098711274139E+81, 6.48559505908724E+82, 5.07580213877225E+83, 3.98864096133865E+84,
    3.14699732603879E+85, 2.49290060083666E+86, 1.98260831540444E+87, 1.58299188153128E+88
};

/* Precomputed inverse factorials 1/(d/2)! indexed by d */
static const double factorial_2_inverse[128] = {
	1.00000000000000E+00, 1.12837916709551E+00, 1.00000000000000E+00, 7.52252778063675E-01,
    5.00000000000000E-01, 3.00901111225470E-01, 1.66666666666667E-01, 8.59717460644200E-02,
    4.16666666666667E-02, 1.91048324587600E-02, 8.33333333333333E-03, 3.47360590159273E-03,
    1.38888888888889E-03, 5.34400907937343E-04, 1.98412698412698E-04, 7.12534543916457E-05,
    2.48015873015873E-05, 8.38275934019361E-06, 2.75573192239859E-06, 8.82395720020380E-07,
    2.75573192239859E-07, 8.40376876209886E-08, 2.50521083854417E-08, 7.30762501052075E-09,
    2.08767569878681E-09, 5.84610000841660E-10, 1.60590438368216E-10, 4.33044445067892E-11,
    1.14707455977297E-11, 2.98651341426135E-12, 7.64716373181982E-13, 1.92678284791055E-13,
    4.77947733238738E-14, 1.16774718055185E-14, 2.81145725434552E-15, 6.67284103172485E-16,
    1.56192069685862E-16, 3.60694109822965E-17, 8.22063524662433E-18, 1.84971338370751E-18,
    4.11031762331217E-19, 9.02299211564640E-20, 1.95729410633913E-20, 4.19674051890530E-21,
    8.89679139245058E-22, 1.86521800840236E-22, 3.86817017063068E-23, 7.93709790809513E-24,
    1.61173757109612E-24, 3.23963179922250E-25, 6.44695028438447E-26, 1.27044384283235E-26,
    2.47959626322480E-27, 4.79412770880134E-28, 9.18368986379555E-29, 1.74331916683685E-29,
    3.27988923706984E-30, 6.11690935732228E-31, 1.13099628864477E-31, 2.07352859570247E-32,
    3.76998762881591E-33, 6.79845441213924E-34, 1.21612504155352E-34, 2.15823949591722E-35,
    3.80039075485474E-36, 6.64073691051452E-37, 1.15163356207720E-37, 1.98230952552672E-38,
    3.38715753552116E-39, 5.74582471167166E-40, 9.67759295863189E-41, 1.61854217230188E-41,
    2.68822026628664E-42, 4.43436211589555E-43, 7.26546017915307E-44, 1.18249656423881E-44,
    1.91196320504028E-45, 3.07141964737354E-46, 4.90246975651354E-47, 7.77574594271782E-48,
    1.22561743912839E-48, 1.91993726980687E-49, 2.98931082714241E-50, 4.62635486700451E-51,
    7.11740673129144E-52, 1.08855408635400E-52, 1.65521086774220E-53, 2.50242318702069E-54,
    3.76184288123226E-55, 5.62342289218133E-56, 8.35965084718281E-57, 1.23591711916073E-57,
    1.81731540156148E-58, 2.65788627776502E-59, 3.86662851396059E-60, 5.59555005845266E-61,
    8.05547607075124E-62, 1.15372166153663E-62, 1.64397470831658E-63, 2.33075083138714E-64,
    3.28794941663316E-65, 4.61534818096463E-66, 6.44695964045717E-67, 8.96184112808665E-68,
    1.23979993085715E-68, 1.70701735773079E-69, 2.33924515256066E-70, 3.19068664996409E-71,
    4.33193546770492E-72, 5.85447091736531E-73, 7.87624630491804E-74, 1.05485962475051E-74,
    1.40647255444965E-75, 1.86700818539913E-76, 2.46749570956079E-77, 3.24697075721587E-78,
    4.25430294751860E-79, 5.55037736276217E-80, 7.21068296189594E-81, 9.32836531556668E-82,
    1.20178049364932E-82, 1.54187856455648E-83, 1.97013195680217E-84, 2.50711961716500E-85,
    3.17763218839060E-86, 4.01139138746400E-87, 5.04386061649301E-88, 6.31715179128189E-89
};

/* scale the input by 2 to index precomputed arrays */
static inline double fact(double x);
static inline double fact(double x)
{
    return factorial_2[(int)(2*x)];
}

/* scale the input by 2 to index precomputed arrays */
static inline double factin(double x);
static inline double factin(double x)
{
    return factorial_2_inverse[(int)(2*x)];
}

double clebsch_gordan(double j1, double m1, double j2, double m2, double j, double m)
{
	int x, xmin, xmax;
    double factor, sum, factsum;
	
    /* constrain the inputs by mathematical conditions */
	if (j1 < 0 || j2 < 0 || j < 0)
        return 0;
	if (j1+j2 < j || abs(j1-j2) > j || (int)(2*(j1+j2+j))%2 != 0)
        return 0;
	if (abs(m1) > j1 || abs(m2) > j2 || abs(m) > j || (int)(2*(j1+m1))%2 != 0 || (int)(2*(j2+m2))%2 != 0 || (int)(2*(j+m))%2 != 0)
        return 0;
	if (m1+m2 != m)
        return 0;

    /* use scientific notation to compute factorials and the factor */
	factor = fact(j1+j2-j) * fact(j1-j2+j) * fact(-j1+j2+j) * fact(j1+m1) * fact(j1-m1) * fact(j2+m2) * fact(j2-m2) * fact(j+m) * fact(j-m) * factin(j1+j2+j+1);
	factor = sqrt((2*j+1)*factor);
	
    xmin = -j+j2-m1 > 0 ? -j+j2-m1 : 0;
    xmin = -j+j1+m2 > xmin ? -j+j1+m2 : xmin;
    xmax = j1-m1 < j1+j2-j ? j1-m1 : j1+j2-j;
    xmax = j2+m2 < xmax ? j2+m2 : xmax;
    if (xmax-xmin <= -1)
        return 0;
    
    /* compute sums */
    sum = 0;
    for (x = xmin; x <= xmax; x += 1) {
        factsum = factin(x) * factin(j1+j2-j-x) * factin(j1-m1-x) * factin(j2+m2-x) * factin(j-j2+m1+x) * factin(j-j1-m2+x);
		if (x%2 != 0)
			factsum = -factsum;
		sum += factsum;
	}

	return factor * sum;
}

double threeJ_symbol(double j1, double m1, double j2, double m2, double j, double m)
{
	return pow(-1, j1-j2-m) / sqrt(2*j+1) * clebsch_gordan(j1, m1, j2, m2, j, -m);
}

double sixJ_symbol(double j1, double j2, double j12, double j3, double j, double j23)
{
	int x, xmin, xmax;
    double factor, sum, factsum;
	
    /* constrain the inputs by mathematical conditions */
	if (j1 < 0 || j2 < 0 || j12 < 0 || j3 < 0 || j < 0 || j23 < 0)
        return 0;
	if (j1+j2 < j12 || abs(j1-j2) > j12 || (int)(2*(j1+j2+j12))%2 != 0)
        return 0;
    if (j2+j3 < j23 || abs(j2-j3) > j23 || (int)(2*(j2+j3+j23))%2 != 0)
        return 0;
    if (j12+j3 < j || abs(j12-j3) > j || (int)(2*(j12+j3+j))%2 != 0)
        return 0;
    if (j1+j23 < j || abs(j1-j23) > j || (int)(2*(j1+j23+j))%2 != 0)
        return 0;

    /* use scientific notation to compute factorials and the factor */
	factor = fact(j1+j2-j12) * fact(j1-j2+j12) * fact(-j1+j2+j12) * factin(j1+j2+j12+1) * fact(j2+j3-j23) * fact(j2-j3+j23) * fact(-j2+j3+j23) * factin(j2+j3+j23+1)
        * fact(j12+j3-j) * fact(j12-j3+j) * fact(-j12+j3+j) * factin(j12+j3+j+1) * fact(j1+j23-j) * fact(j1-j23+j) * fact(-j1+j23+j) * factin(j1+j23+j+1);
	factor = sqrt(factor);

    xmin = j12+j3+j > j1+j2+j12 ? j12+j3+j : j1+j2+j12;
    xmin = j1+j23+j > xmin ? j1+j23+j : xmin;
    xmin = j2+j3+j23 > xmin ? j2+j3+j23 : xmin;
    xmax = j1+j12+j3+j23 < j1+j2+j3+j ? j1+j12+j3+j23 : j1+j2+j3+j;
    xmax = j2+j12+j+j23 < xmax ? j2+j12+j+j23 : xmax;
	
    /* compute sums */
	sum = 0;
	for (x = xmin; x <= xmax; x += 1) {
        factsum = fact(x+1) * factin(x-j1-j2-j12) * factin(x-j12-j3-j) * factin(x-j1-j23-j) * factin(x-j2-j3-j23) 
            * factin(j1+j2+j3+j-x) * factin(j1+j12+j3+j23-x) * factin(j2+j12+j+j23-x);
		if (x%2 != 0)
			factsum = -factsum;   
		sum += factsum;
	}
	
	return factor * sum;
}

double nineJ_symbol(double j1, double j2, double j12, double j3, double j4, double j34, double j13, double j24, double j)
{
	int x, xmin, xmax;
	double sum;

    /* constrain the inputs by mathematical conditions */
	if (j1 < 0 || j2 < 0 || j12 < 0 || j3 < 0 || j4 < 0 || j34 < 0 || j13 < 0 || j24 < 0 || j < 0)
        return 0;
	if (j1+j2 < j12 || abs(j1-j2) > j12 || (int)(2*(j1+j2+j12))%2 != 0)
        return 0;
	if (j3+j4 < j34 || abs(j3-j4) > j34 || (int)(2*(j3+j4+j34))%2 != 0)
        return 0;
	if (j1+j3 < j13 || abs(j1-j3) > j13 || (int)(2*(j1+j3+j13))%2 != 0)
        return 0;
	if (j2+j4 < j24 || abs(j2-j4) > j24 || (int)(2*(j2+j4+j24))%2 != 0)
        return 0;
	if (j12+j34 < j || abs(j12-j34) > j || (int)(2*(j12+j34+j))%2 != 0)
        return 0;
	if (j13+j24 < j || abs(j13-j24) > j || (int)(2*(j13+j24+j))%2 != 0)
        return 0;

	xmin = abs((int)(j2-j34)) > abs((int)(j3-j24)) ? abs((int)(2*(j2-j34))) : abs((int)(2*(j3-j24)));
	xmax = abs((int)(j2+j34)) < abs((int)(j3+j24)) ? abs((int)(2*(j2+j34))) : abs((int)(2*(j3+j24)));

    /* compute sums */
	sum = 0;
	for (x = xmin; x <= xmax; x += 1) {
		sum += pow(-1, x) * (x+1) 
            * sixJ_symbol(j1, j2, j12, j34, j, 0.5*x)
            * sixJ_symbol(j3, j4, j34, j2, 0.5*x, j24)
            * sixJ_symbol(j13, j24, j, 0.5*x, j1, j3);
	}

	return sum;
}

double operator_center_sl(double s1, double s2, double s, double l, double s1p, double s2p, double sp, double lp, double j)
{
    if (s1 != s1p || s2 != s2p || s != sp || l != lp) {
        return 0.0;
    }
    else {
        return 1.0;
    }
}

double operator_sdots_sl(double s1, double s2, double s, double l, double s1p, double s2p, double sp, double lp, double j)
{
    if (s1 != s1p || s2 != s2p || s != sp || l != lp) {
        return 0.0;
    }

    double pre_factor = (s*(s+1)-s1*(s1+1)-s2*(s2+1));

    return 0.5 * pre_factor;
}

double operator_ldots1_sl(double s1, double s2, double s, double l, double s1p, double s2p, double sp, double lp, double j)
{
    if (s1 != s1p || s2 != s2p || l != lp) {
        return 0.0;
    }

    double phase = pow(-1, s+lp+j) * pow(-1, s1p+s2p+s+1);
    double sqrt_term = sqrt((2*sp+1)*(2*s+1)) * sqrt((2*s1+1)*(s1+1)*s1) * sqrt((2*l+1)*(l+1)*l);
    double symbol_term = sixJ_symbol(sp, s, 1, l, lp, j) * sixJ_symbol(s1, s1p, 1, sp, s, s2p);

    return phase * sqrt_term * symbol_term;
}

double operator_ldots2_sl(double s1, double s2, double s, double l, double s1p, double s2p, double sp, double lp, double j)
{
    if (s1 != s1p || s2 != s2p || l != lp) {
        return 0.0;
    }

    double phase = pow(-1, s+lp+j) * pow(-1, s1+s2+sp+1);
    double sqrt_term = sqrt((2*sp+1)*(2*s+1)) * sqrt((2*s2+1)*(s2+1)*s2) * sqrt((2*l+1)*(l+1)*l);
    double symbol_term = sixJ_symbol(sp, s, 1, l, lp, j) * sixJ_symbol(s2, s2p, 1, sp, s, s1p);

    return phase * sqrt_term * symbol_term;
}

double operator_tensor_sl(double s1, double s2, double s, double l, double s1p, double s2p, double sp, double lp, double j)
{
    if (s1 != s1p || s2 != s2p) {
        return 0.0;
    }

    double phase = pow(-1, s+lp+j) * pow(-1, lp);
    double sqrt_term = sqrt(4.8*M_PI) * sqrt(5*(2*sp+1)*(2*s+1)) * sqrt((2*s1+1)*(s1+1)*s1) * sqrt((2*s2+1)*(s2+1)*s2) * sqrt(1.25*M_1_PI*(2*lp+1)*(2*l+1));
    double symbol_term = sixJ_symbol(sp, s, 2, l, lp, j) * nineJ_symbol(s1p, s1, 1, s2p, s2, 1, sp, s, 2) * threeJ_symbol(lp, 0, 2, 0, l, 0);

    return phase * sqrt_term * symbol_term;
}

/* transform sl coupling into jj coupling */
static double trans_sl_jj(double s1, double s2, double L, double jl, double s1p, double s2p, double Lp, double jlp, double J, operator_sl osl);
static double trans_sl_jj(double s1, double s2, double L, double jl, double s1p, double s2p, double Lp, double jlp, double J, operator_sl osl)
{
    double result = 0.0;

    /* S loop for bra */
    double S_min = fabs(s1 - s2);
    double S_max = s1 + s2;
    for (double S = S_min; S <= S_max + 1e-10; S += 1.0) {
        /* Sp loop for ket */
        double Sp_min = fabs(s1p - s2p);
        double Sp_max = s1p + s2p;
        for (double Sp = Sp_min; Sp <= Sp_max + 1e-10; Sp += 1.0) {
            /* phase factor for bra */
            double phase_bra = pow(-1.0, s2 + L + S + jl);
            double factor_bra = sqrt((2.0 * S + 1.0) * (2.0 * jl + 1.0));
            double sixj_bra   = sixJ_symbol(s2, s1, S, L, J, jl);

            /* phase factor for ket */
            double phase_ket = pow(-1.0, s2p + Lp + Sp + jlp);
            double factor_ket = sqrt((2.0 * Sp + 1.0) * (2.0 * jlp + 1.0));
            double sixj_ket   = sixJ_symbol(s2p, s1p, Sp, Lp, J, jlp);

            /* operator value in sl coupling */
            double sl_value = osl(s1, s2, S, L, s1p, s2p, Sp, Lp, J);

            result += phase_bra * factor_bra * sixj_bra *
                      phase_ket * factor_ket * sixj_ket *
                      sl_value;
        }
    }

    return result;
}

double operator_center_jj(double s1, double s2, double l, double jl, double s1p, double s2p, double lp, double jlp, double j)
{
    return trans_sl_jj(s1, s2, l, jl, s1p, s2p, lp, jlp, j, operator_center_sl);
}

double operator_sdots_jj(double s1, double s2, double l, double jl, double s1p, double s2p, double lp, double jlp, double j)
{
    return trans_sl_jj(s1, s2, l, jl, s1p, s2p, lp, jlp, j, operator_sdots_sl);
}

double operator_ldots1_jj(double s1, double s2, double l, double jl, double s1p, double s2p, double lp, double jlp, double j)
{
    return trans_sl_jj(s1, s2, l, jl, s1p, s2p, lp, jlp, j, operator_ldots1_sl);
}

double operator_ldots2_jj(double s1, double s2, double l, double jl, double s1p, double s2p, double lp, double jlp, double j)
{
    return trans_sl_jj(s1, s2, l, jl, s1p, s2p, lp, jlp, j, operator_ldots2_sl);
}

double operator_tensor_jj(double s1, double s2, double l, double jl, double s1p, double s2p, double lp, double jlp, double j)
{
    return trans_sl_jj(s1, s2, l, jl, s1p, s2p, lp, jlp, j, operator_tensor_sl);
}
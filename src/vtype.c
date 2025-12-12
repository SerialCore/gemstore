/*
 * Copyright (C) 2025, Wen-Xuan Zhang <serialcore@outlook.com>
 * Copyright (C) 2025, Si-Qiang Luo <luosq15@lzu.edu.cn>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <vtype.h>
#include <spin.h>
#include <basis.h>
#include <matrix.h>
#include <sumckdk.h>

#include <math.h>
#include <unistd.h>
#include <pthread.h>

static double factorial[171] = {
	1.00000000000000000000E+000, 1.00000000000000000000E+000, 2.00000000000000000000E+000, 6.00000000000000000000E+000, 2.40000000000000000000E+001,
	1.20000000000000000000E+002, 7.20000000000000000000E+002, 5.04000000000000000000E+003, 4.03200000000000000000E+004, 3.62880000000000000000E+005,
	3.62880000000000000000E+006, 3.99168000000000000000E+007, 4.79001600000000000000E+008, 6.22702080000000000000E+009, 8.71782912000000000000E+010,
	1.30767436800000000000E+012, 2.09227898880000000000E+013, 3.55687428096000000000E+014, 6.40237370572800000000E+015, 1.21645100408832000000E+017,
	2.43290200817664000000E+018, 5.10909421717094400000E+019, 1.12400072777760768000E+021, 2.58520167388849766400E+022, 6.20448401733239439360E+023,
	1.55112100433309859840E+025, 4.03291461126605635584E+026, 1.08888694504183521607E+028, 3.04888344611713860501E+029, 8.84176199373970195454E+030,
	2.65252859812191058636E+032, 8.22283865417792281772E+033, 2.63130836933693530167E+035, 8.68331761881188649551E+036, 2.95232799039604140847E+038,
	1.03331479663861449296E+040, 3.71993326789901217468E+041, 1.37637530912263450463E+043, 5.23022617466601111760E+044, 2.03978820811974433586E+046,
	8.15915283247897734345E+047, 3.34525266131638071081E+049, 1.40500611775287989854E+051, 6.04152630633738356373E+052, 2.65827157478844876804E+054,
	1.19622220865480194561E+056, 5.50262215981208894985E+057, 2.58623241511168180642E+059, 1.24139155925360726708E+061, 6.08281864034267560872E+062,
	3.04140932017133780436E+064, 1.55111875328738228022E+066, 8.06581751709438785716E+067, 4.27488328406002556429E+069, 2.30843697339241380472E+071,
	1.26964033536582759259E+073, 7.10998587804863451854E+074, 4.05269195048772167556E+076, 2.35056133128287857182E+078, 1.38683118545689835737E+080,
	8.32098711274139014427E+081, 5.07580213877224798800E+083, 3.14699732603879375256E+085, 1.98260831540444006411E+087, 1.26886932185884164103E+089,
	8.24765059208247066672E+090, 5.44344939077443064003E+092, 3.64711109181886852882E+094, 2.48003554243683059960E+096, 1.71122452428141311372E+098,
	1.19785716699698917960E+100, 8.50478588567862317521E+101, 6.12344583768860868615E+103, 4.47011546151268434089E+105, 3.30788544151938641225E+107,
	2.48091408113953980919E+109, 1.88549470166605025498E+111, 1.45183092028285869634E+113, 1.13242811782062978314E+115, 8.94618213078297528685E+116,
	7.15694570462638022948E+118, 5.79712602074736798588E+120, 4.75364333701284174842E+122, 3.94552396972065865119E+124, 3.31424013456535326699E+126,
	2.81710411438055027694E+128, 2.42270953836727323817E+130, 2.10775729837952771721E+132, 1.85482642257398439114E+134, 1.65079551609084610812E+136,
	1.48571596448176149730E+138, 1.35200152767840296255E+140, 1.24384140546413072554E+142, 1.15677250708164157475E+144, 1.08736615665674308027E+146,
	1.03299784882390592626E+148, 9.91677934870949689209E+149, 9.61927596824821198533E+151, 9.42689044888324774562E+153, 9.33262154439441526817E+155,
	9.33262154439441526817E+157, 9.42594775983835942085E+159, 9.61446671503512660926E+161, 9.90290071648618040754E+163, 1.02990167451456276238E+166,
	1.08139675824029090050E+168, 1.14628056373470835453E+170, 1.22652020319613793935E+172, 1.32464181945182897449E+174, 1.44385958320249358220E+176,
	1.58824554152274294042E+178, 1.76295255109024466387E+180, 1.97450685722107402353E+182, 2.23119274865981364659E+184, 2.54355973347218755712E+186,
	2.92509369349301569068E+188, 3.39310868445189820119E+190, 3.96993716080872089540E+192, 4.68452584975429065657E+194, 5.57458576120760588132E+196,
	6.68950291344912705758E+198, 8.09429852527344373968E+200, 9.87504420083360136241E+202, 1.21463043670253296757E+205, 1.50614174151114087979E+207,
	1.88267717688892609974E+209, 2.37217324288004688567E+211, 3.01266001845765954481E+213, 3.85620482362580421735E+215, 4.97450422247728744039E+217,
	6.46685548922047367250E+219, 8.47158069087882051098E+221, 1.11824865119600430745E+224, 1.48727070609068572890E+226, 1.99294274616151887673E+228,
	2.69047270731805048359E+230, 3.65904288195254865768E+232, 5.01288874827499166103E+234, 6.91778647261948849222E+236, 9.61572319694108900419E+238,
	1.34620124757175246058E+241, 1.89814375907617096942E+243, 2.69536413788816277658E+245, 3.85437071718007277052E+247, 5.55029383273930478955E+249,
	8.04792605747199194484E+251, 1.17499720439091082394E+254, 1.72724589045463891120E+256, 2.55632391787286558858E+258, 3.80892263763056972698E+260,
	5.71338395644585459047E+262, 8.62720977423324043162E+264, 1.31133588568345254560E+267, 2.00634390509568239477E+269, 3.08976961384735088795E+271,
	4.78914290146339387633E+273, 7.47106292628289444708E+275, 1.17295687942641442819E+278, 1.85327186949373479654E+280, 2.94670227249503832650E+282,
	4.71472363599206132240E+284, 7.59070505394721872907E+286, 1.22969421873944943411E+289, 2.00440157654530257759E+291, 3.28721858553429622726E+293,
	5.42391066613158877498E+295, 9.00369170577843736647E+297, 1.50361651486499904020E+300, 2.52607574497319838753E+302, 4.26906800900470527493E+304,
	7.25741561530799896739E+306};

void getijk(double x1, double x2, double x3, double *xi, double *xj, double *xk, int c)
{
    switch (c) {
    case 1:
        *xi = x1;
        *xj = x2;
        *xk = x3;
        break;
    case 2:
        *xi = x3;
        *xj = x1;
        *xk = x2;
        break;
    case 3:
        *xi = x2;
        *xj = x3;
        *xk = x1;
        break;
    default:
        break;
    }
}

void get123(double *x1, double *x2, double *x3, double xi, double xj, double xk, int c)
{
    switch (c) {
    case 1:
        *x1 = xi;
        *x2 = xj;
        *x3 = xk;
        break;
    case 2:
        *x1 = xj;
        *x2 = xk;
        *x3 = xi;
        break;
    case 3:
        *x1 = xk;
        *x2 = xi;
        *x3 = xj;
        break;
    default:
        break;
    }
}

void sumckdk_scdk_vtype(sumckdk_scdk *scdk, vtype vsodt, int fg, basis_list qnlist, matrix_t **mlsj, int nf, int nfp, int ni, int nip, pthread_mutex_t *mutex, int lock, int c)
{
    int f, i, mf, mg;
    double msa1, msa2, msa3, msb1, msb2, msb3, msai, msaj, msak, msbi, msbj, msbk;
    int l1, l2, l3, l4, ml1, ml2, ml3, ml4;
    double coef, coei;

    if (1 == lock) {
        pthread_mutex_lock(mutex);
    }
    sumckdk_scdk_init(scdk);
    if (1 == lock) {
        pthread_mutex_unlock(mutex);
    }

    l1 = qnlist.qnum[nf][nfp].lrho;
    l2 = qnlist.qnum[nf][nfp].llam;
    l3 = qnlist.qnum[ni][nip].lrho;
    l4 = qnlist.qnum[ni][nip].llam;

    for (f = 0; f < mlsj[nf][nfp].n; f++) {
        for (i = 0; i < mlsj[ni][nip].n; i++) {
            coef = mlsj[nf][nfp].p[f][0];
            msa1 = mlsj[nf][nfp].p[f][1];
            msa2 = mlsj[nf][nfp].p[f][2];
            msa3 = mlsj[nf][nfp].p[f][3];
            ml1 = (int)mlsj[nf][nfp].p[f][4];
            ml2 = (int)mlsj[nf][nfp].p[f][5];

            coei = mlsj[ni][nip].p[i][0];
            msb1 = mlsj[ni][nip].p[i][1];
            msb2 = mlsj[ni][nip].p[i][2];
            msb3 = mlsj[ni][nip].p[i][3];
            ml3 = (int)mlsj[ni][nip].p[i][4];
            ml4 = (int)mlsj[ni][nip].p[i][5];

            getijk(msa1, msa2, msa3, &msai, &msaj, &msak, c);
            getijk(msb1, msb2, msb3, &msbi, &msbj, &msbk, c);

            vsodt(scdk, coef, coei, msai, msaj, msak, msbi, msbj, msbk, l1, l2, l3, l4, ml1, ml2, ml3, ml4, qnlist.qnum[nf][nfp], qnlist.qnum[ni][nip], mutex, lock, c);
        }
    }

    for (i = 0; i < scdk->len; i++) {
        mf = scdk->nfgh[i][0];
        mg = scdk->nfgh[i][5];
        switch (fg) {
        case 1:
            scdk->coe[i] *= 1 / factorial[2 * mf + 1] / factorial[mg];
            scdk->nfgh[i][18] = 2 * mf + 2;
            break;

        case 2:
            scdk->coe[i] *= factorial[mf + 2] / (factorial[2 * mf + 5] * factorial[mf]) / factorial[mg];
            scdk->nfgh[i][18] = 2 * mf + 4;
            break;

        case 3:
            scdk->coe[i] *= factorial[mf + 1] / (factorial[2 * mf + 3] * factorial[mf]) / factorial[mg];
            scdk->nfgh[i][18] = 2 * mf + 4;
            break;

        default:
            break;
        }
    }
}

void vcent(sumckdk_scdk *scdk, double coef, double coei, double msai, double msaj, double msak, double msbi, double msbj, double msbk, int l1, int l2, int l3, int l4, int ml1, int ml2, int ml3, int ml4, basis_qnum qa, basis_qnum qb, pthread_mutex_t *mutex, int lock, int c)
{
    double pi = 3.1415926535897932385;
    if (msai == msbi && msaj == msbj && msak == msbk) {
        sumckdk_scdk_calc(scdk, 4 * pi * coef * coei, l1, l2, l3, l4, 0, 0, ml1, ml2, ml3, ml4, 0, 0, mutex, lock);
    }
}

void vcont(sumckdk_scdk *scdk, double coef, double coei, double msai, double msaj, double msak, double msbi, double msbj, double msbk, int l1, int l2, int l3, int l4, int ml1, int ml2, int ml3, int ml4, basis_qnum qa, basis_qnum qb, pthread_mutex_t *mutex, int lock, int c)
{
    int M;
    double pi = 3.1415926535897932385;
    double coe = 0;
    double si = 0.5, sj = 0.5;
    if (msak == msbk) {
        for (M = -1; M <= 1; M++) {
            coe = coe - sqrt(3) * sqrt(si * (si + 1)) * sqrt(sj * (sj + 1)) * clebsch_gordan(1, M, 1, -M, 0, 0) * clebsch_gordan(si, msbi, 1, M, si, msai) * clebsch_gordan(sj, msbj, 1, -M, sj, msaj);
        }
    }
    if (0 != coe) {
        sumckdk_scdk_calc(scdk, 4 * pi * coe * coef * coei, l1, l2, l3, l4, 0, 0, ml1, ml2, ml3, ml4, 0, 0, mutex, lock);
    }
}

void vtens(sumckdk_scdk *scdk, double coef, double coei, double msai, double msaj, double msak, double msbi, double msbj, double msbk, int l1, int l2, int l3, int l4, int ml1, int ml2, int ml3, int ml4, basis_qnum qa, basis_qnum qb, pthread_mutex_t *mutex, int lock, int c)
{
    int M, mu;
    double pi = 3.1415926535897932385;
    double coe;
    double si = 0.5, sj = 0.5;
    if (msak == msbk) {
        for (M = -2; M <= 2; M++) {
            for (mu = -1; mu <= 1; mu++) {
                coe = clebsch_gordan(2, -M, 2, M, 0, 0) * clebsch_gordan(1, mu, 1, -M - mu, 2, -M) * sqrt(si * (si + 1)) * sqrt(sj * (sj + 1)) * clebsch_gordan(si, msbi, 1, mu, si, msai) * clebsch_gordan(sj, msbj, 1, -M - mu, sj, msaj);
                if (0 != coe) {
                    sumckdk_scdk_calc(scdk, sqrt(4 * pi) * coe * coef * coei, l1, l2, l3, l4, 2, 0, ml1, ml2, ml3, ml4, M, 0, mutex, lock);
                }
            }
        }
    }
}

void vsorp(sumckdk_scdk *scdk, double coef, double coei, double msai, double msaj, double msak, double msbi, double msbj, double msbk, int l1, int l2, int l3, int l4, int ml1, int ml2, int ml3, int ml4, basis_qnum qa, basis_qnum qb, pthread_mutex_t *mutex, int lock, int c)
{
    int M, mu;
    double coe;
    double si = 0.5, sj = 0.5;
    double mi, mj, mk;
    getijk(qa.m1, qa.m2, qa.m3, &mi, &mj, &mk, c);
    if (msak == msbk) {
        for (M = -1; M <= 1; M++) {
            for (mu = -1; mu <= 1; mu++) {
                coe = clebsch_gordan(1, -M, 1, M, 0, 0) * clebsch_gordan(1, mu, 1, M, 1, mu + M) * (sqrt(si * (si + 1)) * clebsch_gordan(si, msbi, 1, -M, si, msai) + sqrt(sj * (sj + 1)) * clebsch_gordan(sj, msbj, 1, -M, sj, msaj));
                if (0 != coe) {
                    sumckdk_scdk_calc(scdk, pow(-1, mu) * coe * coef * coei, l1, l2, l3, l4, 1, 1, ml1, ml2, ml3, ml4, -mu, mu + M, mutex, lock);
                }
            }
        }
    }
}

void vsorn(sumckdk_scdk *scdk, double coef, double coei, double msai, double msaj, double msak, double msbi, double msbj, double msbk, int l1, int l2, int l3, int l4, int ml1, int ml2, int ml3, int ml4, basis_qnum qa, basis_qnum qb, pthread_mutex_t *mutex, int lock, int c)
{
    int M, mu;
    double coe;
    double si = 0.5, sj = 0.5;
    double mi, mj, mk;
    getijk(qa.m1, qa.m2, qa.m3, &mi, &mj, &mk, c);
    if (msak == msbk) {
        for (M = -1; M <= 1; M++) {
            for (mu = -1; mu <= 1; mu++) {
                coe = clebsch_gordan(1, -M, 1, M, 0, 0) * clebsch_gordan(1, mu, 1, M, 1, mu + M) * (sqrt(si * (si + 1)) * clebsch_gordan(si, msbi, 1, -M, si, msai) - sqrt(sj * (sj + 1)) * clebsch_gordan(sj, msbj, 1, -M, sj, msaj));
                if (0 != coe) {
                    sumckdk_scdk_calc(scdk, pow(-1, mu) * coe * coef * coei, l1, l2, l3, l4, 1, 1, ml1, ml2, ml3, ml4, -mu, mu + M, mutex, lock);
                }
            }
        }
    }
}

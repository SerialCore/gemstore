/*
 * Copyright (C) 2025, Wen-Xuan Zhang <serialcore@outlook.com>
 * Copyright (C) 2025, Si-Qiang Luo <luosq15@lzu.edu.cn>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <sumckdk.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

static inline double almj(int l, int m, int j);
static inline double almj(int l, int m, int j)
{
	return sqrt((2 * l + 1) * factorial[l - m] / factorial[l + m]) * factorial[l + m] * pow(-1, j + m) * pow(2, -1 - 2 * j - m) / (sqrt(M_PI) * factorial[j] * factorial[l - 2 * j - m] * factorial[j + m]);
}

static inline double binomial(int n, int m);
static inline double binomial(int n, int m)
{
	return factorial[n] / (factorial[n - m] * factorial[m]);
}

void sumckdk_cdlm_init(sumckdk_cdlm *cdlm)
{
	cdlm->len = 0;
	cdlm->c = (double *)malloc(sizeof(double) * 0);
	cdlm->d = (double **)malloc(sizeof(double *) * 0);
}

void sumckdk_cdlm_logs(sumckdk_cdlm cdlm)
{
	int i;
	for (i = 0; i < cdlm.len; i++) {
		printf("%4d: %15.10f | %15.10f %15.10f %15.10f\n", i, cdlm.c[i], cdlm.d[i][0], cdlm.d[i][1], cdlm.d[i][2]);
	}
}

void sumckdk_cdlm_push(sumckdk_cdlm *cdlm, double c, double dx, double dy, double dz)
{
	cdlm->c = (double *)realloc(cdlm->c, sizeof(double) * (cdlm->len + 1));
	cdlm->c[cdlm->len] = c;
	cdlm->d = (double **)realloc(cdlm->d, sizeof(double *) * (cdlm->len + 1));
	cdlm->d[cdlm->len] = (double *)malloc(sizeof(double) * 3);
	cdlm->d[cdlm->len][0] = dx;
	cdlm->d[cdlm->len][1] = dy;
	cdlm->d[cdlm->len][2] = dz;
	cdlm->len++;
}

void sumckdk_cdlm_cali(sumckdk_cdlm *cdlm, int l, int m)
{
	double pi = 3.1415926535897932385;
	int j, p, q, s, t, u;
	double c;
	sumckdk_cdlm_init(cdlm);
	if (0 == l && 0 == m) {
		sumckdk_cdlm_push(cdlm, 1 / sqrt(4 * pi), 0, 0, 0);
		return;
	}
	for (j = 0; j <= (l - m) / 2; j++) {
		p = l - m - 2 * j;
		q = m + j;
		for (s = 0; s <= p; s++) {
			for (t = 0; t <= q; t++) {
				for (u = 0; u <= j; u++) {
					c = pow(l * 0.25, l) * almj(l, m, j) * binomial(p, s) * binomial(q, t) * binomial(j, u) * pow(-1, s + t + u);
					sumckdk_cdlm_push(cdlm, c, (j + q - 2 * t - 2 * u) * 1.0 / l, +(-j + q - 2 * t + 2 * u) * 1.0 / l, (p - 2 * s) * 1.0 / l);
				}
			}
		}
	}
}

void sumckdk_cdlm_calf(sumckdk_cdlm *cdlm, int l, int m)
{
	double pi = 3.1415926535897932385;
	int j, p, q, s, t, u;
	double c;
	sumckdk_cdlm_init(cdlm);
	if (0 == l && 0 == m) {
		sumckdk_cdlm_push(cdlm, 1 / sqrt(4 * pi), 0, 0, 0);
		return;
	}
	for (j = 0; j <= (l - m) / 2; j++) {
		p = l - m - 2 * j;
		q = m + j;
		for (s = 0; s <= p; s++) {
			for (t = 0; t <= q; t++) {
				for (u = 0; u <= j; u++) {
					c = pow(l * 0.25, l) * almj(l, m, j) * binomial(p, s) * binomial(q, t) * binomial(j, u) * pow(-1, s + t + u);
					sumckdk_cdlm_push(cdlm, c, (j + q - 2 * t - 2 * u) * 1.0 / l, -(-j + q - 2 * t + 2 * u) * 1.0 / l, (p - 2 * s) * 1.0 / l);
				}
			}
		}
	}
}

void sumckdk_cdlm_free(sumckdk_cdlm *cdlm)
{
	int i;
	for (i = 0; i < cdlm->len; i++) {
		free(cdlm->d[i]);
	}
	free(cdlm->d);
	free(cdlm->c);
	cdlm->len = 0;
}

void sumckdk_scdk_init(sumckdk_scdk *scdk)
{
	scdk->len = 0;
	scdk->coe = (double *)malloc(sizeof(double) * 0);
	scdk->nfgh = (int **)malloc(sizeof(int *) * 0);
}

void sumckdk_scdk_free(sumckdk_scdk *scdk)
{
	int i;
	for (i = 0; i < scdk->len; i++) {
		free(scdk->nfgh[i]);
	}
	free(scdk->nfgh);
	free(scdk->coe);
	scdk->len = 0;
}

void sumckdk_scdk_logs(sumckdk_scdk scdk)
{
	int i;
	for (i = 0; i < scdk.len; i++) {
		printf("%4d: %18.10E |nf:%2d f:%2d %2d %2d %2d |ng:%2d g:%2d %2d %2d %2d|h1= %2d %2d %2d %2d |h2= %2d %2d %2d %2d|int=%2d\n", i, scdk.coe[i], scdk.nfgh[i][0], scdk.nfgh[i][1], scdk.nfgh[i][2], scdk.nfgh[i][3], scdk.nfgh[i][4], scdk.nfgh[i][5], scdk.nfgh[i][6], scdk.nfgh[i][7], scdk.nfgh[i][8], scdk.nfgh[i][9], scdk.nfgh[i][10], scdk.nfgh[i][11], scdk.nfgh[i][12], scdk.nfgh[i][13], scdk.nfgh[i][14], scdk.nfgh[i][15], scdk.nfgh[i][16], scdk.nfgh[i][17], scdk.nfgh[i][18]);
	}
}

void sumckdk_scdk_push(sumckdk_scdk *scdk, double coe, int nf, int nf1, int nf2, int nf3, int nf4, int ng, int ng1, int ng2, int ng3, int ng4, int nh11, int nh12, int nh13, int nh14, int nh21, int nh22, int nh23, int nh24)
{
	int i, push;
	push = 1;
	for (i = 0; i < scdk->len; i++)
	{
		if (
			scdk->nfgh[i][0] == nf &&
			scdk->nfgh[i][1] == nf1 &&
			scdk->nfgh[i][2] == nf2 &&
			scdk->nfgh[i][3] == nf3 &&
			scdk->nfgh[i][4] == nf4 &&
			scdk->nfgh[i][5] == ng &&
			scdk->nfgh[i][6] == ng1 &&
			scdk->nfgh[i][7] == ng2 &&
			scdk->nfgh[i][8] == ng3 &&
			scdk->nfgh[i][9] == ng4 &&
			scdk->nfgh[i][10] == nh11 &&
			scdk->nfgh[i][11] == nh12 &&
			scdk->nfgh[i][12] == nh13 &&
			scdk->nfgh[i][13] == nh14 &&
			scdk->nfgh[i][14] == nh21 &&
			scdk->nfgh[i][15] == nh22 &&
			scdk->nfgh[i][16] == nh23 &&
			scdk->nfgh[i][17] == nh24 &&
			1 == 1)
		{
			scdk->coe[i] += coe;
			push = 0;
			break;
		}
	}
	if (1 == push)
	{
		scdk->coe = (double *)realloc(scdk->coe, sizeof(double) * (scdk->len + 1));
		scdk->coe[scdk->len] = coe;
		scdk->nfgh = (int **)realloc(scdk->nfgh, sizeof(int *) * (scdk->len + 1));
		scdk->nfgh[scdk->len] = (int *)malloc(sizeof(int) * 19);
		scdk->nfgh[scdk->len][0] = nf;
		scdk->nfgh[scdk->len][1] = nf1;
		scdk->nfgh[scdk->len][2] = nf2;
		scdk->nfgh[scdk->len][3] = nf3;
		scdk->nfgh[scdk->len][4] = nf4;
		scdk->nfgh[scdk->len][5] = ng;
		scdk->nfgh[scdk->len][6] = ng1;
		scdk->nfgh[scdk->len][7] = ng2;
		scdk->nfgh[scdk->len][8] = ng3;
		scdk->nfgh[scdk->len][9] = ng4;
		scdk->nfgh[scdk->len][10] = nh11;
		scdk->nfgh[scdk->len][11] = nh12;
		scdk->nfgh[scdk->len][12] = nh13;
		scdk->nfgh[scdk->len][13] = nh14;
		scdk->nfgh[scdk->len][14] = nh21;
		scdk->nfgh[scdk->len][15] = nh22;
		scdk->nfgh[scdk->len][16] = nh23;
		scdk->nfgh[scdk->len][17] = nh24;
		scdk->nfgh[scdk->len][18] = -1;
		scdk->len++;
	}
}

void sumckdk_scdk_spfy(sumckdk_scdk *scdk)
{
	int i, j, k;
	double eps = 1.1102230246251565404E-16;
	j = -1;
	for (i = 0; i < scdk->len; i++)
	{
		if (fabs(scdk->coe[i]) > eps)
		{
			j++;
			if (i != j)
			{
				scdk->coe[j] = scdk->coe[i];
				for (k = 0; k < 18; k++)
				{
					scdk->nfgh[j][k] = scdk->nfgh[i][k];
				}
			}
		}
	}
	for (i = j + 1; i < scdk->len; i++)
	{
		free(scdk->nfgh[i]);
	}

	scdk->len = j + 1;
	scdk->coe = (double *)realloc(scdk->coe, sizeof(double) * (scdk->len));
	scdk->nfgh = (int **)realloc(scdk->nfgh, sizeof(int *) * (scdk->len));
}

void sumckdk_scdk_calc(sumckdk_scdk *scdk, double cit, int l1, int l2, int l3, int l4, int lh1, int lh2, int m1, int m2, int m3, int m4, int mh1, int mh2, pthread_mutex_t *mutex, int lock)
{
	int n;
	sumckdk_cdlm cdlm1, cdlm2, cdlm3, cdlm4;
	int nf, ng;
	int f12, f13, f14, f23, f24, f34;
	int g12, g13, g14, g23, g24, g34;
	int k1, k2, k3, k4, k;
	int nf1, nf2, nf3, nf4;
	int ng1, ng2, ng3, ng4;
	int sign;
	double d12, d13, d14, d23, d24, d34;
	double coe, sc;
	int j1, j2;
	int ha11, ha12, ha13, ha14;
	int hb11, hb12, hb13, hb14;
	int hc11, hc12, hc13, hc14;
	int ha21, ha22, ha23, ha24;
	int hb21, hb22, hb23, hb24;
	int hc21, hc22, hc23, hc24;
	int nh11, nh12, nh13, nh14;
	int nh21, nh22, nh23, nh24;

	if (1 == lock)
	{
		pthread_mutex_lock(mutex);
	}
	sumckdk_cdlm_calf(&cdlm1, l1, m1);
	sumckdk_cdlm_calf(&cdlm2, l2, m2);
	sumckdk_cdlm_cali(&cdlm3, l3, m3);
	sumckdk_cdlm_cali(&cdlm4, l4, m4);
	if (1 == lock)
	{
		pthread_mutex_unlock(mutex);
	}

	if ((l1 + l2 + l3 + l4 - lh1 - lh2) % 2 != 0)
	{
		printf("error_sumckdk_scdk_calc\n");
		return;
	}
	n = (l1 + l2 + l3 + l4 - lh1 - lh2) / 2;
	for (nf = 0; nf <= n; nf++)
	{
		ng = n - nf;
		for (f12 = 0; f12 <= nf; f12++)
		{
			for (f13 = 0; f13 <= nf - f12; f13++)
			{
				for (f14 = 0; f14 <= nf - f12 - f13; f14++)
				{
					for (f23 = 0; f23 <= nf - f12 - f13 - f14; f23++)
					{
						for (f24 = 0; f24 <= nf - f12 - f13 - f14 - f23; f24++)
						{
							f34 = nf - f12 - f13 - f14 - f23 - f24;

							for (g12 = 0; g12 <= ng; g12++)
							{
								for (g13 = 0; g13 <= ng - g12; g13++)
								{
									for (g14 = 0; g14 <= ng - g12 - g13; g14++)
									{
										for (g23 = 0; g23 <= ng - g12 - g13 - g14; g23++)
										{
											for (g24 = 0; g24 <= ng - g12 - g13 - g14 - g23; g24++)
											{
												g34 = ng - g12 - g13 - g14 - g23 - g24;

												for (j1 = 0; j1 <= (lh1 - mh1) / 2; j1++)
												{

													for (ha11 = 0; ha11 <= lh1 - mh1 - 2 * j1; ha11++)
													{
														for (ha12 = 0; ha12 <= lh1 - mh1 - 2 * j1 - ha11; ha12++)
														{
															for (ha13 = 0; ha13 <= lh1 - mh1 - 2 * j1 - ha11 - ha12; ha13++)
															{
																ha14 = lh1 - mh1 - 2 * j1 - ha11 - ha12 - ha13;

																for (hb11 = 0; hb11 <= mh1 + j1; hb11++)
																{
																	for (hb12 = 0; hb12 <= mh1 + j1 - hb11; hb12++)
																	{
																		for (hb13 = 0; hb13 <= mh1 + j1 - hb11 - hb12; hb13++)
																		{
																			hb14 = mh1 + j1 - hb11 - hb12 - hb13;

																			for (hc11 = 0; hc11 <= j1; hc11++)
																			{
																				for (hc12 = 0; hc12 <= j1 - hc11; hc12++)
																				{
																					for (hc13 = 0; hc13 <= j1 - hc11 - hc12; hc13++)
																					{
																						hc14 = j1 - hc11 - hc12 - hc13;

																						for (j2 = 0; j2 <= (lh2 - mh2) / 2; j2++)
																						{
																							for (ha21 = 0; ha21 <= lh2 - mh2 - 2 * j2; ha21++)
																							{
																								for (ha22 = 0; ha22 <= lh2 - mh2 - 2 * j2 - ha21; ha22++)
																								{
																									for (ha23 = 0; ha23 <= lh2 - mh2 - 2 * j2 - ha21 - ha22; ha23++)
																									{
																										ha24 = lh2 - mh2 - 2 * j2 - ha21 - ha22 - ha23;

																										for (hb21 = 0; hb21 <= mh2 + j2; hb21++)
																										{
																											for (hb22 = 0; hb22 <= mh2 + j2 - hb21; hb22++)
																											{
																												for (hb23 = 0; hb23 <= mh2 + j2 - hb21 - hb22; hb23++)
																												{
																													hb24 = mh2 + j2 - hb21 - hb22 - hb23;

																													for (hc21 = 0; hc21 <= j2; hc21++)
																													{
																														for (hc22 = 0; hc22 <= j2 - hc21; hc22++)
																														{
																															for (hc23 = 0; hc23 <= j2 - hc21 - hc22; hc23++)
																															{
																																hc24 = j2 - hc21 - hc22 - hc23;

																																nf1 = f12 + f13 + f14;
																																ng1 = g12 + g13 + g14;
																																nh11 = ha11 + hb11 + hc11;
																																nh21 = ha21 + hb21 + hc21;
																																nf2 = f12 + f23 + f24;
																																ng2 = g12 + g23 + g24;
																																nh12 = ha12 + hb12 + hc12;
																																nh22 = ha22 + hb22 + hc22;
																																nf3 = f13 + f23 + f34;
																																ng3 = g13 + g23 + g34;
																																nh13 = ha13 + hb13 + hc13;
																																nh23 = ha23 + hb23 + hc23;
																																nf4 = f14 + f24 + f34;
																																ng4 = g14 + g24 + g34;
																																nh14 = ha14 + hb14 + hc14;
																																nh24 = ha24 + hb24 + hc24;
																																if (nf1 + ng1 + nh11 + nh21 == l1 && nf2 + ng2 + nh12 + nh22 == l2 && nf3 + ng3 + nh13 + nh23 == l3 && nf4 + ng4 + nh14 + nh24 == l4)
																																{
																																	coe = binomial(nf, f12) * binomial(nf - f12, f13) * binomial(nf - f12 - f13, f14) * binomial(nf - f12 - f13 - f14, f23) * binomial(nf - f12 - f13 - f14 - f23, f24) *
																																		  binomial(ng, g12) * binomial(ng - g12, g13) * binomial(ng - g12 - g13, g14) * binomial(ng - g12 - g13 - g14, g23) * binomial(ng - g12 - g13 - g14 - g23, g24) *
																																		  binomial(lh1 - mh1 - 2 * j1, ha11) * binomial(lh1 - mh1 - 2 * j1 - ha11, ha12) * binomial(lh1 - mh1 - 2 * j1 - ha11 - ha12, ha13) *
																																		  binomial(lh2 - mh2 - 2 * j2, ha21) * binomial(lh2 - mh2 - 2 * j2 - ha21, ha22) * binomial(lh2 - mh2 - 2 * j2 - ha21 - ha22, ha23) *
																																		  binomial(mh1 + j1, hb11) * binomial(mh1 + j1 - hb11, hb12) * binomial(mh1 + j1 - hb11 - hb12, hb13) *
																																		  binomial(mh2 + j2, hb21) * binomial(mh2 + j2 - hb21, hb22) * binomial(mh2 + j2 - hb21 - hb22, hb23) *
																																		  binomial(j1, hc11) * binomial(j1 - hc11, hc12) * binomial(j1 - hc11 - hc12, hc13) *
																																		  binomial(j2, hc21) * binomial(j2 - hc21, hc22) * binomial(j2 - hc21 - hc22, hc23);

																																	sc = 0;
																																	for (k1 = 0; k1 < cdlm1.len; k1++)
																																	{
																																		for (k2 = 0; k2 < cdlm2.len; k2++)
																																		{
																																			for (k3 = 0; k3 < cdlm3.len; k3++)
																																			{
																																				for (k4 = 0; k4 < cdlm4.len; k4++)
																																				{
																																					d12 = 0;
																																					d13 = 0;
																																					d14 = 0;
																																					d23 = 0;
																																					d24 = 0;
																																					d34 = 0;
																																					for (k = 0; k < 3; k++)
																																					{
																																						if (1 == k)
																																						{
																																							sign = -1;
																																						}
																																						else
																																						{
																																							sign = +1;
																																						}
																																						d12 += sign * cdlm1.d[k1][k] * cdlm2.d[k2][k];
																																						d13 += sign * cdlm1.d[k1][k] * cdlm3.d[k3][k];
																																						d14 += sign * cdlm1.d[k1][k] * cdlm4.d[k4][k];
																																						d23 += sign * cdlm2.d[k2][k] * cdlm3.d[k3][k];
																																						d24 += sign * cdlm2.d[k2][k] * cdlm4.d[k4][k];
																																						d34 += sign * cdlm3.d[k3][k] * cdlm4.d[k4][k];
																																					}
																																					sc += cdlm1.c[k1] * cdlm2.c[k2] * cdlm3.c[k3] * cdlm4.c[k4] *

																																						  pow(cdlm1.d[k1][2], ha11) *
																																						  pow(cdlm2.d[k2][2], ha12) *
																																						  pow(cdlm3.d[k3][2], ha13) *
																																						  pow(cdlm4.d[k4][2], ha14) *
																																						  pow(cdlm1.d[k1][2], ha21) *
																																						  pow(cdlm2.d[k2][2], ha22) *
																																						  pow(cdlm3.d[k3][2], ha23) *
																																						  pow(cdlm4.d[k4][2], ha24) *

																																						  pow(cdlm1.d[k1][0] - cdlm1.d[k1][1], hb11) *
																																						  pow(cdlm2.d[k2][0] - cdlm2.d[k2][1], hb12) *
																																						  pow(cdlm3.d[k3][0] - cdlm3.d[k3][1], hb13) *
																																						  pow(cdlm4.d[k4][0] - cdlm4.d[k4][1], hb14) *
																																						  pow(cdlm1.d[k1][0] - cdlm1.d[k1][1], hb21) *
																																						  pow(cdlm2.d[k2][0] - cdlm2.d[k2][1], hb22) *
																																						  pow(cdlm3.d[k3][0] - cdlm3.d[k3][1], hb23) *
																																						  pow(cdlm4.d[k4][0] - cdlm4.d[k4][1], hb24) *

																																						  pow(cdlm1.d[k1][0] + cdlm1.d[k1][1], hc11) *
																																						  pow(cdlm2.d[k2][0] + cdlm2.d[k2][1], hc12) *
																																						  pow(cdlm3.d[k3][0] + cdlm3.d[k3][1], hc13) *
																																						  pow(cdlm4.d[k4][0] + cdlm4.d[k4][1], hc14) *
																																						  pow(cdlm1.d[k1][0] + cdlm1.d[k1][1], hc21) *
																																						  pow(cdlm2.d[k2][0] + cdlm2.d[k2][1], hc22) *
																																						  pow(cdlm3.d[k3][0] + cdlm3.d[k3][1], hc23) *
																																						  pow(cdlm4.d[k4][0] + cdlm4.d[k4][1], hc24) *

																																						  pow(d12, f12 + g12) *
																																						  pow(d13, f13 + g13) *
																																						  pow(d14, f14 + g14) *
																																						  pow(d23, f23 + g23) *
																																						  pow(d24, f24 + g24) *
																																						  pow(d34, f34 + g34);
																																				}
																																			}
																																		}
																																	}

																																	if (1 == lock)
																																	{
																																		pthread_mutex_lock(mutex);
																																	}
																																	sumckdk_scdk_push(scdk, cit * coe * almj(lh1, mh1, j1) * almj(lh2, mh2, j2) * sc * pow(2, nf + ng), nf, nf1, nf2, nf3, nf4, ng, ng1, ng2, ng3, ng4, nh11, nh12, nh13, nh14, nh21, nh22, nh23, nh24);
																																	if (1 == lock)
																																	{
																																		pthread_mutex_unlock(mutex);
																																	}
																																}
																															}
																														}
																													}
																												}
																											}
																										}
																									}
																								}
																							}
																						}
																					}
																				}
																			}
																		}
																	}
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	pthread_mutex_lock(mutex);
	sumckdk_cdlm_free(&cdlm1);
	sumckdk_cdlm_free(&cdlm2);
	sumckdk_cdlm_free(&cdlm3);
	sumckdk_cdlm_free(&cdlm4);
	sumckdk_scdk_spfy(scdk);
	pthread_mutex_unlock(mutex);
}

double sumckdk_scdk_vals(sumckdk_scdk scdk, double coef, double f1, double f2, double f3, double f4, double coeg, double g1, double g2, double g3, double g4, double h11, double h12, double h13, double h14, double h21, double h22, double h23, double h24)
{
	int i;
	double sum;
	sum = 0;
	for (i = 0; i < scdk.len; i++)
	{
		sum += scdk.coe[i] *
			   pow(coef, scdk.nfgh[i][0]) *
			   sin(scdk.nfgh[i][0]) *
			   pow(f1, scdk.nfgh[i][1]) *
			   pow(f2, scdk.nfgh[i][2]) *
			   pow(f3, scdk.nfgh[i][3]) *
			   pow(f4, scdk.nfgh[i][4]) *
			   pow(coeg, scdk.nfgh[i][5]) *
			   cos(scdk.nfgh[i][5]) *
			   pow(g1, scdk.nfgh[i][6]) *
			   pow(g2, scdk.nfgh[i][7]) *
			   pow(g3, scdk.nfgh[i][8]) *
			   pow(g4, scdk.nfgh[i][9]) *
			   pow(h11, scdk.nfgh[i][10]) *
			   pow(h12, scdk.nfgh[i][11]) *
			   pow(h13, scdk.nfgh[i][12]) *
			   pow(h14, scdk.nfgh[i][13]) *
			   pow(h21, scdk.nfgh[i][14]) *
			   pow(h22, scdk.nfgh[i][15]) *
			   pow(h23, scdk.nfgh[i][16]) *
			   pow(h24, scdk.nfgh[i][17]) *
			   1;
	}
	return sum;
}

/*
 * Copyright (C) 2026, Wen-Xuan Zhang <serialcore@outlook.com>
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Radial integrals of the baryon GI operators (recycle inteCenV).
 * deltaij/betaij use 0.5+epsilon. They do not add the epsilon^2 term in GIVdelta*.
 */

#include <gemstore/model/giradial.h>
#include <gemstore/basis/channel.h>
#include <gemstore/model/gimodel.h>

#include <math.h>
#include <stdio.h>

void scdk_vargs_from_model(scdk_vargs_t *varg, const argsGIModel_t *model, model_type_t mtype)
{
    varg->mn = model->mn;
    varg->ms = model->ms;
    varg->mc = model->mc;
    varg->mb = model->mb;
    varg->alpha[0] = GI_ALPHA_K[0];
    varg->alpha[1] = GI_ALPHA_K[1];
    varg->alpha[2] = GI_ALPHA_K[2];
    varg->gamma[0] = GI_GAMMA_K[0];
    varg->gamma[1] = GI_GAMMA_K[1];
    varg->gamma[2] = GI_GAMMA_K[2];
    varg->b = model->b;
    varg->c = model->c;
    varg->sigma0 = model->sigma_0;
    varg->s = model->s;
    varg->f = 0.5;
    varg->mu = model->mu;
    varg->econt = model->epsilon_cont;
    varg->etens = model->epsilon_tens;
    varg->esov = model->epsilon_sov;
    varg->esos = model->epsilon_sos;
    varg->eCoul = 0.0;
    varg->model = mtype;
}


typedef double (*vCenPart)(double,basis_3b,basis_3b,scdk_vargs_t,int);

static double getsigmaij(double sigma0,double s,double mi,double mj)
{
	return sqrt(sigma0 * sigma0 * (0.5 + 0.5 * pow(4.0 * mi * mj / ((mi + mj) * (mi + mj)), 4.0))+ s * s * pow(2.0 * mi * mj / (mi + mj), 2.0));
}

static void getsigmakij(double *sigmakij, double *gammak,double sigma0,double s,double mi,double mj)
{
	int k;
	double sigmaij2=pow(getsigmaij(sigma0,s,mi,mj),2);
	for(k=0;k<3;k++)
	{
		sigmakij[k]=1/sqrt(1/pow(gammak[k],2)+1/sigmaij2);
	}
}

/* √(p²+m²). c = 1, 2, 3 reads qf.m1, m2, m3. qi is unused. */
static double Ti(double p,basis_3b qf,basis_3b qi,scdk_vargs_t varg,int c)
{
	double mc;
	switch(c)
	{
		case 1:
			mc=qf.m1;
			break;
		case 2:
			mc=qf.m2;
			break;
		case 3:
			mc=qf.m3;
			break;
		default:
			printf("error_Ti\n");
			mc=0;
			break;
	}
	return  sqrt(p*p+mc*mc);
}


/* 50-point recycle quadrature for ∫ x^n exp(-b11 x^2) K(x) dx.
 * The nodes are divided by √b11. This table is not integral_exp_rn. */
static double InteCenV(double b11,int n,basis_3b qf,basis_3b qi,vCenPart v,scdk_vargs_t varg,int c)
{
	#define InteLength 50
	
	const double InteWeight[InteLength]={
		 9.49074370469973154861E-003, 2.20248042432282550986E-002, 3.43744122779471232838E-002, 4.62954423315795426704E-002, 5.74270673524102611643E-002,
		 6.72616743347232737151E-002, 7.51655105427606561651E-002, 8.04487210189303589897E-002, 8.24866653780293862703E-002, 8.08772841479732575406E-002,
		 7.55968131214373958980E-002, 6.70988969096135470989E-002, 5.63037514101928791197E-002, 4.44533647270852744933E-002, 3.28602655529441193691E-002,
		 2.26281047022441036531E-002, 1.44421205638529698147E-002, 8.49988917441698033564E-003, 4.58981582947691516137E-003, 2.26248032155746286264E-003,
		 1.01296345857313010366E-003, 4.09855051025779821045E-004, 1.49103396554347145977E-004, 4.85209725035065972790E-005, 1.40499187954310526217E-005,
		 3.60057445584128798598E-006, 8.12054885660543383779E-007, 1.60236226323028537512E-007, 2.74914498210976715853E-008, 4.07393870688591939297E-009,
		 5.17733909405469218567E-010, 5.59883137228935952164E-011, 5.10833440280922385331E-012, 3.89531180859741328527E-013, 2.45631708318249480239E-014,
		 1.26561330926208323536E-015, 5.25582653443700296954E-017, 1.73149373479240142712E-018, 4.44198243358375823341E-020, 8.68013859892778184465E-022,
		 1.25805074779480909285E-023, 1.30871663315222277155E-025, 9.37690669168311366587E-028, 4.38351942427205493247E-030, 1.25991102592426549802E-032,
		 1.45161317869704793270E-035, 1.44829917427803858346E-037,-2.07564613115271761842E-039, 1.59488572462612746697E-041,-3.26209228139812174140E-044
	};
	const double InteNode[InteLength]={
		 3.69941668941189387078E-003, 1.94655785593001262377E-002, 4.77232575110775569884E-002, 8.82994741423326136810E-002, 1.40944249831512837831E-001,
		 2.05343824732975154536E-001, 2.81130944731155938791E-001, 3.67895701867100816092E-001, 4.65196632382567480466E-001, 5.72571575839739193266E-001,
		 6.89547864293400636314E-001, 8.15651524129035456870E-001, 9.50415293767051842487E-001, 1.09338537085220496150E+000, 1.24412689394047714663E+000,
		 1.40222823252836353434E+000, 1.56730420579590844586E+000, 1.73899837732924138424E+000, 1.91698458425367595630E+000, 2.10096785886658098596E+000,
		 2.29068489295538163123E+000, 2.48590418287233886279E+000, 2.68642597976929227008E+000, 2.89208215616667258410E+000, 3.10273608869143607692E+000,
		 3.31828264842635672287E+000, 3.53864838570112628139E+000, 3.76379199610129929525E+000, 3.99370515987453091261E+000, 4.22841385900133938902E+000,
		 4.46798029678176976527E+000, 4.71250557663335036490E+000, 4.96213334417320913961E+000, 5.21705466625021294844E+000, 5.47751452299329293275E+000,
		 5.74382044124152601111E+000, 6.01635402812180604868E+000, 6.29558651982983683470E+000, 6.58210002638651010927E+000, 6.87661707965657488199E+000,
		 7.18004266516692291795E+000, 7.49352570514240963406E+000, 7.81855215039574067889E+000, 8.15709210450168878089E+000, 8.51184526474371849140E+000,
		 8.88668005924412897894E+000, 9.28749674141648604654E+000, 9.72416586588463146083E+000, 1.02158862585784281522E+001, 1.08129860729453608573E+001
	};
	double r,res;
	int i;
	res=0;
	for(i=0;i<InteLength;i++)
	{
		r=InteNode[i]/sqrt(b11);
		res=res+InteWeight[i]*v(r,qf,qi,varg,c)*pow(r,n);
	}
	res=res/sqrt(b11);
	return res;
}

/* ∫_0^∞ r^n exp(-a r^2) dr = (1/2) a^{-(n+1)/2} Γ((n+1)/2).
 * GammaSemiInt[n+1] stores that gamma value. */
static double InteExprn(double a,int n)
{
	double GammaSemiInt[100]={
 0.00000000000000000000e+000, 1.77245385090551602729e+000, 1.00000000000000000000e+000, 8.86226925452758013649e-001, 1.00000000000000000000e+000, 1.32934038817913702047e+000, 2.00000000000000000000e+000, 3.32335097044784255118e+000, 6.00000000000000000000e+000, 1.16317283965674489291e+001, 2.40000000000000000000e+001, 5.23427777845535201811e+001, 1.20000000000000000000e+002, 2.87885277815044360996e+002, 7.20000000000000000000e+002, 1.87125430579778834647e+003, 5.04000000000000000000e+003, 1.40344072934834125985e+004, 4.03200000000000000000e+004, 1.19292461994609007087e+005,
 3.62880000000000000000e+005, 1.13327838894878556733e+006, 3.62880000000000000000e+006, 1.18994230839622484570e+007, 3.99168000000000000000e+007, 1.36843365465565857255e+008, 4.79001600000000000000e+008, 1.71054206831957321569e+009, 6.22702080000000000000e+009, 2.30923179223142384118e+010, 8.71782912000000000000e+010, 3.34838609873556456972e+011, 1.30767436800000000000e+012, 5.18999845304012508307e+012, 2.09227898880000000000e+013, 8.56349744751620638707e+013, 3.55687428096000000000e+014, 1.49861205331533611773e+015, 6.40237370572800000000e+015, 2.77243229863337181781e+016,
 1.21645100408832000000e+017, 5.40624298233507504473e+017, 2.43290200817664000000e+018, 1.10827981137869038417e+019, 5.10909421717094400000e+019, 2.38280159446418432596e+020, 1.12400072777760768000e+021, 5.36130358754441473342e+021, 2.58520167388849766400e+022, 1.25990634307293746235e+023, 6.20448401733239439360e+023, 3.08677054052869678277e+024, 1.55112100433309859840e+025, 7.87126487834817679606e+025, 4.03291461126605635584e+026, 2.08588519276226685095e+027, 1.08888694504183521607e+028, 5.73618428009623384013e+028, 3.04888344611713860501e+029, 1.63481251982742664443e+030,
 8.84176199373970195454e+030, 4.82269693349090860109e+031, 2.65252859812191058636e+032, 1.47092256471472712333e+033, 8.22283865417792281772e+033, 4.63340607885139043849e+034, 2.63130836933693530167e+035, 1.50585697562670189251e+036, 8.68331761881188649551e+036, 5.04462086834945133991e+037, 2.95232799039604140847e+038, 1.74039419958056071227e+039, 1.03331479663861449296e+040, 6.17839940851099052856e+040, 3.71993326789901217468e+041, 2.25511578410651154292e+042, 1.37637530912263450463e+043, 8.45668419039941828596e+043, 5.23022617466601111760e+044, 3.25582341330377604009e+045,
 2.03978820811974433586e+046, 1.28605024825499153583e+047, 8.15915283247897734345e+047, 5.20850350543271572014e+048, 3.34525266131638071081e+049, 2.16152895475457702386e+050, 1.40500611775287989854e+051, 9.18649805770695235140e+051, 6.04152630633738356373e+052, 3.99612665510252427286e+053, 2.65827157478844876804e+054, 1.77827636152062330142e+055, 1.19622220865480194561e+056, 8.09115744491883602147e+056, 5.50262215981208894985e+057, 3.76238821188725874998e+058, 2.58623241511168180642e+059, 1.78713440064644790624e+060, 1.24139155925360726708e+061, 8.66760184313527234528e+061};
	return 0.5e0*pow(a,-0.5*(1+n))*GammaSemiInt[n+1];
}

/* Closed form of ∫ r^n exp(-a r^2) erf(b r) dr for n = 1, 3, 5, 7, 9.
 * Any other n prints and returns 0. */
static double InteExpErfrn(double a, double b, int n)
{
    double result = 0.0;
    
    switch(n)
    {
        case 1:
            result = b / (2.0 * a * sqrt(a + b * b));
            break;
            
        case 3:
            result = (3.0 * a * b + 2.0 * pow(b, 3)) / (4.0 * pow(a, 2) * pow(a + b * b, 1.5));
            break;
            
        case 5:
            result = (15.0 * pow(a, 2) * b + 20.0 * a * pow(b, 3) + 8.0 * pow(b, 5)) / 
                     (8.0 * pow(a, 3) * pow(a + b * b, 2.5));
            break;
            
        case 7:
            result = (3.0 * b * (35.0 * pow(a, 3) + 70.0 * pow(a, 2) * pow(b, 2) + 
                     56.0 * a * pow(b, 4) + 16.0 * pow(b, 6))) / 
                     (16.0 * pow(a, 4) * pow(a + b * b, 3.5));
            break;
            
        case 9:
            result = (3.0 * b * (315.0 * pow(a, 4) + 840.0 * pow(a, 3) * pow(b, 2) + 
                     1008.0 * pow(a, 2) * pow(b, 4) + 576.0 * a * pow(b, 6) + 128.0 * pow(b, 8))) / 
                     (32.0 * pow(a, 5) * pow(a + b * b, 4.5));
            break;
            
        default:
            printf("only n=1,3,5,7,9\n");
            result = 0.0;
            break;
    }
    
    return result;
}

/* Callbacks below. Names and the GISTRING / GISCREEN split are in giradial.h. */
double inteNfi(double b11,int n,basis_3b qf,basis_3b qi,void *ctx,int c)
{
	(void)ctx;

	return InteExprn(b11,n);
}

double inteVogeG(double b11,int n,basis_3b qf,basis_3b qi,void *ctx,int c)
{
	scdk_vargs_t varg = *(const scdk_vargs_t *)ctx;

	double sigmakij[3];
	int k;
	double mi,mj,mk;
	double res=0;
	channel_3b_ijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	getsigmakij(sigmakij,varg.gamma,varg.sigma0,varg.s,mi,mj);
	for(k=0;k<3;k++)
	{
		res=res-2*varg.alpha[k]*InteExpErfrn(b11,sigmakij[k],n-1)/3.0;
	}
	return res;
}

double inteVcont(double b11,int n,basis_3b qf,basis_3b qi,void *ctx,int c)
{
	scdk_vargs_t varg = *(const scdk_vargs_t *)ctx;

	double pi=3.1415926535897932385;
	double sigmakij[3];
	int k;
	double mi,mj,mk;
	double res=0;
	channel_3b_ijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	getsigmakij(sigmakij,varg.gamma,varg.sigma0,varg.s,mi,mj);
	for(k=0;k<3;k++)
	{
		res=res+16*varg.alpha[k]*pow(sigmakij[k],3)*InteExprn(pow(sigmakij[k],2)+b11,n)/(9*mi*mj*sqrt(pi));
	}
	return res;
}

double inteVtens(double b11,int n,basis_3b qf,basis_3b qi,void *ctx,int c)
{
	scdk_vargs_t varg = *(const scdk_vargs_t *)ctx;

	double pi=3.1415926535897932385;
	double sigmakij[3];
	int k;
	double mi,mj,mk;
	double res=0;
	channel_3b_ijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	getsigmakij(sigmakij,varg.gamma,varg.sigma0,varg.s,mi,mj);
	for(k=0;k<3;k++)
	{
		res=res-4*varg.alpha[k]*sigmakij[k]*InteExprn(pow(sigmakij[k],2)+b11,n-2)/(3*mi*mj*sqrt(pi))-8*varg.alpha[k]*pow(sigmakij[k],3)*InteExprn(pow(sigmakij[k],2)+b11,n)/(9*mi*mj*sqrt(pi))+2*varg.alpha[k]*InteExpErfrn(b11,sigmakij[k],n-3)/(3*mi*mj);

	}
	return res;
}

double inteVsovii(double b11,int n,basis_3b qf,basis_3b qi,void *ctx,int c)
{
	scdk_vargs_t varg = *(const scdk_vargs_t *)ctx;

	double pi=3.1415926535897932385;
	double sigmakij[3];
	int k;
	double mi,mj,mk;
	double res=0;
	channel_3b_ijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	getsigmakij(sigmakij,varg.gamma,varg.sigma0,varg.s,mi,mj);
	for(k=0;k<3;k++)
	{
		res=res-2*varg.alpha[k]*sigmakij[k]*InteExprn(pow(sigmakij[k],2)+b11,n-2)/(3*sqrt(pi)*mi*mi)+varg.alpha[k]*InteExpErfrn(b11,sigmakij[k],n-3)/(3*mi*mi);
	}
	return res;
}

double inteVsovjj(double b11,int n,basis_3b qf,basis_3b qi,void *ctx,int c)
{
	scdk_vargs_t varg = *(const scdk_vargs_t *)ctx;

	double pi=3.1415926535897932385;
	double sigmakij[3];
	int k;
	double mi,mj,mk;
	double res=0;
	channel_3b_ijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	getsigmakij(sigmakij,varg.gamma,varg.sigma0,varg.s,mi,mj);
	for(k=0;k<3;k++)
	{
		res=res+2*varg.alpha[k]*sigmakij[k]*InteExprn(pow(sigmakij[k],2)+b11,n-2)/(3*sqrt(pi)*mj*mj)-varg.alpha[k]*InteExpErfrn(b11,sigmakij[k],n-3)/(3*mj*mj);
	}
	return res;
}

double inteVsovij(double b11,int n,basis_3b qf,basis_3b qi,void *ctx,int c)
{
	scdk_vargs_t varg = *(const scdk_vargs_t *)ctx;

	double pi=3.1415926535897932385;
	double sigmakij[3];
	int k;
	double mi,mj,mk;
	double res=0;
	channel_3b_ijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	getsigmakij(sigmakij,varg.gamma,varg.sigma0,varg.s,mi,mj);
	for(k=0;k<3;k++)
	{
		res=res-4*varg.alpha[k]*sigmakij[k]*InteExprn(pow(sigmakij[k],2)+b11,n-2)/(3*sqrt(pi)*mi*mj)+2*varg.alpha[k]*InteExpErfrn(b11,sigmakij[k],n-3)/(3*mi*mj);
	}
	return res;
}

double inteVsovji(double b11,int n,basis_3b qf,basis_3b qi,void *ctx,int c)
{
	scdk_vargs_t varg = *(const scdk_vargs_t *)ctx;

	double pi=3.1415926535897932385;
	double sigmakij[3];
	int k;
	double mi,mj,mk;
	double res=0;
	channel_3b_ijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	getsigmakij(sigmakij,varg.gamma,varg.sigma0,varg.s,mi,mj);
	for(k=0;k<3;k++)
	{
		res=res+4*varg.alpha[k]*sigmakij[k]*InteExprn(pow(sigmakij[k],2)+b11,n-2)/(3*sqrt(pi)*mi*mj)-2*varg.alpha[k]*InteExpErfrn(b11,sigmakij[k],n-3)/(3*mi*mj);
	}
	return res;
}

/* Baryon pair colour Cij = −2/3, OCent = 1. Matches recycle c/2 on the constant. */
static double baryon_vconf_r(double r, basis_3b qf, basis_3b qi, scdk_vargs_t varg, int c)
{
	double mi, mj, mk;
	channel_3b_ijk(qf.m1, qf.m2, qf.m3, &mi, &mj, &mk, c);
	double sigmaij = getsigmaij(varg.sigma0, varg.s, mi, mj);
	double Cij = -2.0 / 3.0;
	(void)qi;
	if (r == 0.0) {
		return -0.75 * Cij * varg.c;
	}
	if (varg.model == MODEL_GISCREEN && varg.mu > 0.0) {
		double sig2 = sigmaij * sigmaij;
		double mu = varg.mu;
		double mu2_4sig2 = mu * mu / (4.0 * sig2);
		double mu_m = mu - 2.0 * r * sig2;
		double mu_p = mu + 2.0 * r * sig2;
		double pref = -3.0 * Cij * varg.b * exp(-r * mu) / (16.0 * r * mu * sig2);
		double inner1 = 4.0 * r * sig2 * exp(r * mu);
		double inner2 = mu_m * exp(mu2_4sig2) * erfc(mu_m / (2.0 * sigmaij));
		double inner3 = mu_p * exp(mu2_4sig2 + 2.0 * r * mu) * erfc(mu_p / (2.0 * sigmaij));
		return pref * (inner1 + inner2 - inner3) - 0.75 * Cij * varg.c;
	}
	double rsig = r * sigmaij;
	double pref = -3.0 * Cij * varg.b / (8.0 * r * sigmaij * sigmaij);
	double inner1 = (2.0 / sqrt(3.1415926535897932385)) * r * sigmaij * exp(-rsig * rsig);
	double inner2 = (1.0 + 2.0 * rsig * rsig) * erf(rsig);
	return pref * (inner1 + inner2) - 0.75 * Cij * varg.c;
}

static double baryon_vsosii_r(double r, basis_3b qf, basis_3b qi, scdk_vargs_t varg, int c)
{
	double mi, mj, mk;
	channel_3b_ijk(qf.m1, qf.m2, qf.m3, &mi, &mj, &mk, c);
	if (r == 0.0) {
		return 0.0;
	}
	double h = 1e-6 * (r > 1e-8 ? r : 1.0);
	double dV = (baryon_vconf_r(r + h, qf, qi, varg, c) - baryon_vconf_r(r - h, qf, qi, varg, c)) / (2.0 * h);
	return -dV / (2.0 * r * mi * mi);
}

static double baryon_vsosjj_r(double r, basis_3b qf, basis_3b qi, scdk_vargs_t varg, int c)
{
	double mi, mj, mk;
	channel_3b_ijk(qf.m1, qf.m2, qf.m3, &mi, &mj, &mk, c);
	if (r == 0.0) {
		return 0.0;
	}
	double h = 1e-6 * (r > 1e-8 ? r : 1.0);
	double dV = (baryon_vconf_r(r + h, qf, qi, varg, c) - baryon_vconf_r(r - h, qf, qi, varg, c)) / (2.0 * h);
	return dV / (2.0 * r * mj * mj);
}

double inteVstring(double b11,int n,basis_3b qf,basis_3b qi,void *ctx,int c)
{
	scdk_vargs_t varg = *(const scdk_vargs_t *)ctx;

	double pi=3.1415926535897932385;
	double sigmaij;
	double mi,mj,mk;
	channel_3b_ijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	if (varg.model == MODEL_GISCREEN && varg.mu > 0.0) {
		return InteCenV(b11, n, qf, qi, baryon_vconf_r, varg, c);
	}
	sigmaij=getsigmaij(varg.sigma0,varg.s,mi,mj);
	return varg.b*varg.f*InteExprn(pow(sigmaij,2)+b11,n)/(sqrt(pi)*sigmaij)+varg.b*varg.f*InteExpErfrn(b11,sigmaij,n+1)+varg.b*varg.f*InteExpErfrn(b11,sigmaij,n-1)/(2*pow(sigmaij,2))+varg.c*InteExprn(b11,n)/2;
}

double inteVsosii(double b11,int n,basis_3b qf,basis_3b qi,void *ctx,int c)
{
	scdk_vargs_t varg = *(const scdk_vargs_t *)ctx;

	double pi=3.1415926535897932385;
	double sigmaij;
	double mi,mj,mk;
	channel_3b_ijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	if (varg.model == MODEL_GISCREEN && varg.mu > 0.0) {
		return InteCenV(b11, n, qf, qi, baryon_vsosii_r, varg, c);
	}
	sigmaij=getsigmaij(varg.sigma0,varg.s,mi,mj);
	return -varg.b*varg.f*InteExprn(pow(sigmaij,2)+b11,n-2)/(2*sqrt(pi)*sigmaij*pow(mi,2))-varg.b*varg.f*InteExpErfrn(b11,sigmaij,n-1)/(2*pow(mi,2))+varg.b*varg.f*InteExpErfrn(b11,sigmaij,n-3)/(4*pow(sigmaij,2)*pow(mi,2));
}

double inteVsosjj(double b11,int n,basis_3b qf,basis_3b qi,void *ctx,int c)
{
	scdk_vargs_t varg = *(const scdk_vargs_t *)ctx;

	double pi=3.1415926535897932385;
	double sigmaij;
	double mi,mj,mk;
	channel_3b_ijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	if (varg.model == MODEL_GISCREEN && varg.mu > 0.0) {
		return InteCenV(b11, n, qf, qi, baryon_vsosjj_r, varg, c);
	}
	sigmaij=getsigmaij(varg.sigma0,varg.s,mi,mj);
	return  varg.b*varg.f*InteExprn(pow(sigmaij,2)+b11,n-2)/(2*sqrt(pi)*sigmaij*pow(mj,2))+varg.b*varg.f*InteExpErfrn(b11,sigmaij,n-1)/(2*pow(mj,2))-varg.b*varg.f*InteExpErfrn(b11,sigmaij,n-3)/(4*pow(sigmaij,2)*pow(mj,2));
}

static double betaij(double r,double mi,double mj,double epsilon)
{
	return pow(1+r*r/(sqrt(r*r+mi*mi)*sqrt(r*r+mj*mj)),0.5+epsilon);
}

static double deltaij(double r,double mi,double mj,double epsilon)
{
	return pow(mi*mj/(sqrt(r*r+mi*mi)*sqrt(r*r+mj*mj)),0.5+epsilon);
}

static double pogeG(double r,basis_3b qf,basis_3b qi,scdk_vargs_t varg,int c)
{
	double mi,mj,mk;
	channel_3b_ijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	return betaij(r,mi,mj,varg.eCoul);
}

static double pcont(double r,basis_3b qf,basis_3b qi,scdk_vargs_t varg,int c)
{
	double mi,mj,mk;
	channel_3b_ijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	return deltaij(r,mi,mj,varg.econt);
}

static double ptens(double r,basis_3b qf,basis_3b qi,scdk_vargs_t varg,int c)
{
	double mi,mj,mk;
	channel_3b_ijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	return deltaij(r,mi,mj,varg.etens);
}

static double psovii(double r,basis_3b qf,basis_3b qi,scdk_vargs_t varg,int c)
{
	double mi,mj,mk;
	channel_3b_ijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	return deltaij(r,mi,mi,varg.esov);
}

static double psovjj(double r,basis_3b qf,basis_3b qi,scdk_vargs_t varg,int c)
{
	double mi,mj,mk;
	channel_3b_ijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	return deltaij(r,mj,mj,varg.esov);
}

static double psovji(double r,basis_3b qf,basis_3b qi,scdk_vargs_t varg,int c)
{
	double mi,mj,mk;
	channel_3b_ijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	return deltaij(r,mj,mi,varg.esov);
}

static double psovij(double r,basis_3b qf,basis_3b qi,scdk_vargs_t varg,int c)
{
	double mi,mj,mk;
	channel_3b_ijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	return deltaij(r,mi,mj,varg.esov);
}

static double psosii(double r,basis_3b qf,basis_3b qi,scdk_vargs_t varg,int c)
{
	double mi,mj,mk;
	channel_3b_ijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	return deltaij(r,mi,mi,varg.esos);
}

static double psosjj(double r,basis_3b qf,basis_3b qi,scdk_vargs_t varg,int c)
{
	double mi,mj,mk;
	channel_3b_ijk(qf.m1,qf.m2,qf.m3,&mi,&mj,&mk,c);
	return deltaij(r,mj,mj,varg.esos);
}

double intepogeG(double b11,int n,basis_3b qf,basis_3b qi,void *ctx,int c)
{
	scdk_vargs_t varg = *(const scdk_vargs_t *)ctx;

	return InteCenV(b11,n,qf,qi,pogeG,varg,c);
}

double intepcont(double b11,int n,basis_3b qf,basis_3b qi,void *ctx,int c)
{
	scdk_vargs_t varg = *(const scdk_vargs_t *)ctx;

	return InteCenV(b11,n,qf,qi,pcont,varg,c);
}

double inteptens(double b11,int n,basis_3b qf,basis_3b qi,void *ctx,int c)
{
	scdk_vargs_t varg = *(const scdk_vargs_t *)ctx;

	return InteCenV(b11,n,qf,qi,ptens,varg,c);
}

double intepsovii(double b11,int n,basis_3b qf,basis_3b qi,void *ctx,int c)
{
	scdk_vargs_t varg = *(const scdk_vargs_t *)ctx;

	return InteCenV(b11,n,qf,qi,psovii,varg,c);
}

double intepsovjj(double b11,int n,basis_3b qf,basis_3b qi,void *ctx,int c)
{
	scdk_vargs_t varg = *(const scdk_vargs_t *)ctx;

	return InteCenV(b11,n,qf,qi,psovjj,varg,c);
}

double intepsovji(double b11,int n,basis_3b qf,basis_3b qi,void *ctx,int c)
{
	scdk_vargs_t varg = *(const scdk_vargs_t *)ctx;

	return InteCenV(b11,n,qf,qi,psovji,varg,c);
}

double intepsovij(double b11,int n,basis_3b qf,basis_3b qi,void *ctx,int c)
{
	scdk_vargs_t varg = *(const scdk_vargs_t *)ctx;

	return InteCenV(b11,n,qf,qi,psovij,varg,c);
}

double intepsosii(double b11,int n,basis_3b qf,basis_3b qi,void *ctx,int c)
{
	scdk_vargs_t varg = *(const scdk_vargs_t *)ctx;

	return InteCenV(b11,n,qf,qi,psosii,varg,c);
}

double intepsosjj(double b11,int n,basis_3b qf,basis_3b qi,void *ctx,int c)
{
	scdk_vargs_t varg = *(const scdk_vargs_t *)ctx;

	return InteCenV(b11,n,qf,qi,psosjj,varg,c);
}

double inteTi(double b11,int n,basis_3b qf,basis_3b qi,void *ctx,int c)
{
	scdk_vargs_t varg = *(const scdk_vargs_t *)ctx;

	return InteCenV(b11,n,qf,qi,Ti,varg,c);
}


double inteRMS(double b11,int n,basis_3b qf,basis_3b qi,void *ctx,int c)
{
	(void)ctx;

	return InteExprn(b11,n+2);
}



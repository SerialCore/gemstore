void getq(double mn,double ms,double mc,double mb,int q,double *mi)
{
	switch(q)
	{
		case 1:
			*mi=mn;
			break;

		case 2:
			*mi=ms;
			break;

		case 3:
			*mi=mc;
			break;

		case 4:
			*mi=mb;
			break;
			
		default:
			printf("error_getq\n");
			*mi=0;
			break;
	}
}

void debug01(int q1,int q2,int q3,int f12,double J,int P,int Lmax,double jl)
{
	double GeV=1;
	double fm=5.0676896/GeV;

	double mn=0.220*GeV;
	double ms=0.419*GeV;
	double mc=1.628*GeV;
	double mb=4.977*GeV;

	double alpha1= 0.25;
	double alpha2= 0.15;
	double alpha3= 0.20;
	double gamma1=0.5;
	double gamma2=sqrt(2.5);
	double gamma3=5*sqrt(10.0);

        double      b= 0.141*GeV*GeV;
        double      c=-0.204*GeV;
	double sigma0= 1.889*GeV;
	double      s= 1.422;
	double      f= 0.5;

	double econt=-0.156;
	double etens= 0.379;
	double esov = 0.006;
	double esos = 0.449;
	double eCoul= 0;

	double rmin=0.2*fm;
	double rmax=2.0*fm;
	int nmax=6;

	int i,j;

	double m1,m2,m3;
	char cP;

	margs marg;
	vargs varg;
	
	varg.mn=mn;
	varg.ms=ms;
	varg.mc=mc;
	varg.mb=mb;

	varg.alpha[0]=alpha1;
	varg.alpha[1]=alpha2;
	varg.alpha[2]=alpha3;
	varg.gamma[0]=gamma1;
	varg.gamma[1]=gamma2;
	varg.gamma[2]=gamma3;

	varg.b=b;
	varg.c=c;
	varg.sigma0=sigma0;
	varg.s=s;
	varg.f=f;

	varg.econt=econt;
	varg.etens=etens;
	varg.esov =esov;
	varg.esos =esos;
	varg.eCoul=eCoul;

	marg.varg=varg;

	getq(mn,ms,mc,mb,q1,&m1);
	getq(mn,ms,mc,mb,q2,&m2);
	getq(mn,ms,mc,mb,q3,&m3);

	if(P>0)
	{
		cP='+';
	}
	else if(P<0)
	{
		cP='-';
	}
	else
	{
		cP='?';
	}
	printf("%d/2^%c:",(int)(2*J),cP);

	getlsj_jl(&marg, m1,m2,m3, rmin,rmax,nmax, f12, J,P, Lmax,jl);

//	basis_list_logs(marg.qnlist_spfy);
//	basis_list_logs(marg.qnlist_full);

	malloc_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_cent_1));
	malloc_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_cent_2));
	malloc_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_cent_3));
	malloc_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_cont_1));
	malloc_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_cont_2));
	malloc_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_cont_3));
	malloc_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_tens_1));
	malloc_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_tens_2));
	malloc_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_tens_3));
	malloc_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_soii_1));
	malloc_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_soii_2));
	malloc_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_soii_3));
	malloc_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_sojj_1));
	malloc_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_sojj_2));
	malloc_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_sojj_3));
	malloc_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_soji_1));
	malloc_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_soji_2));
	malloc_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_soji_3));
	malloc_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_soij_1));
	malloc_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_soij_2));
	malloc_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_soij_3));

	mt_load(21*marg.qnlist_spfy.len_list,calc_scdk_mt,&marg,numberProcessors());

	initmatrix(&(marg.Nfi),     marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.VogeG1),  marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.VogeG2),  marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.VogeG3),  marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.Vcont1),  marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.Vcont2),  marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.Vcont3),  marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.Vtens1),  marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.Vtens2),  marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.Vtens3),  marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.Vsovii1), marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.Vsovii2), marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.Vsovii3), marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.Vsovjj1), marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.Vsovjj2), marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.Vsovjj3), marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.Vsovji1), marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.Vsovji2), marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.Vsovji3), marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.Vsovij1), marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.Vsovij2), marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.Vsovij3), marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.Vstring1),marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.Vstring2),marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.Vstring3),marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.Vsosii1), marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.Vsosii2), marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.Vsosii3), marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.Vsosjj1), marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.Vsosjj2), marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.Vsosjj3), marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.pogeG1),  marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.pogeG2),  marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.pogeG3),  marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.pcont1),  marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.pcont2),  marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.pcont3),  marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.ptens1),  marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.ptens2),  marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.ptens3),  marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.psovii1), marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.psovii2), marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.psovii3), marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.psovjj1), marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.psovjj2), marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.psovjj3), marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.psovji1), marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.psovji2), marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.psovji3), marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.psovij1), marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.psovij2), marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.psovij3), marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.psosii1), marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.psosii2), marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.psosii3), marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.psosjj1), marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.psosjj2), marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.psosjj3), marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.T1),      marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.T2),      marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.T3),      marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.rmsr12),  marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.rmsr13),  marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.rmsr23),  marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.rmsl12),  marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.rmsl13),  marg.qnlist_full.len_list,marg.qnlist_full.len_list);
        initmatrix(&(marg.rmsl23),  marg.qnlist_full.len_list,marg.qnlist_full.len_list);
	initmatrix(&(marg.Hfi),     marg.qnlist_full.len_list,marg.qnlist_full.len_list);
	initmatrix(&(marg.u),       marg.qnlist_full.len_list,marg.qnlist_full.len_list);
	initmatrix(&(marg.v),       marg.qnlist_full.len_list,marg.qnlist_full.len_list);

	marg.e1=(double*)malloc(sizeof(double)*marg.qnlist_full.len_list);
	marg.e2=(double*)malloc(sizeof(double)*marg.qnlist_full.len_list);

	mt_load(marg.qnlist_full.len_list,getmfi,&marg,numberProcessors());

	eigsys(&marg);

	free_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_cent_1));
	free_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_cent_2));
	free_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_cent_3));
	free_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_cont_1));
	free_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_cont_2));
	free_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_cont_3));
	free_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_tens_1));
	free_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_tens_2));
	free_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_tens_3));
	free_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_soii_1));
	free_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_soii_2));
	free_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_soii_3));
	free_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_sojj_1));
	free_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_sojj_2));
	free_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_sojj_3));
	free_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_soji_1));
	free_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_soji_2));
	free_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_soji_3));
	free_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_soij_1));
	free_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_soij_2));
	free_scdk(marg.qnlist_spfy.len_part,marg.qnlist_spfy.len_list,&(marg.scdk_soij_3));

	for(i=0;i<marg.qnlist_spfy.len_list;i++)
	{
		for(j=0;j<marg.qnlist_spfy.len_part[i];j++)
		{
			freematrix(&marg.mlsj[i][j]);
		}
		free(marg.mlsj[i]);
	}
	free(marg.mlsj);

	for(i=0;i<marg.qnlist_spfy.len_list;i++)
	{
		free(marg.qnlist_spfy.qnum[i]);
	}
	free(marg.qnlist_spfy.qnum);
	free(marg.qnlist_spfy.len_part);
	
	for(i=0;i<marg.qnlist_full.len_list;i++)
	{
		free(marg.qnlist_full.qnum[i]);
	}
	free(marg.qnlist_full.qnum);
	free(marg.qnlist_full.len_part);
	
	freematrix(&(marg.Nfi)     );
        freematrix(&(marg.VogeG1)  );
        freematrix(&(marg.VogeG2)  );
        freematrix(&(marg.VogeG3)  );
        freematrix(&(marg.Vcont1)  );
        freematrix(&(marg.Vcont2)  );
        freematrix(&(marg.Vcont3)  );
        freematrix(&(marg.Vtens1)  );
        freematrix(&(marg.Vtens2)  );
        freematrix(&(marg.Vtens3)  );
        freematrix(&(marg.Vsovii1) );
        freematrix(&(marg.Vsovii2) );
        freematrix(&(marg.Vsovii3) );
        freematrix(&(marg.Vsovjj1) );
        freematrix(&(marg.Vsovjj2) );
        freematrix(&(marg.Vsovjj3) );
        freematrix(&(marg.Vsovji1) );
        freematrix(&(marg.Vsovji2) );
        freematrix(&(marg.Vsovji3) );
        freematrix(&(marg.Vsovij1) );
        freematrix(&(marg.Vsovij2) );
        freematrix(&(marg.Vsovij3) );
        freematrix(&(marg.Vstring1));
        freematrix(&(marg.Vstring2));
        freematrix(&(marg.Vstring3));
        freematrix(&(marg.Vsosii1) );
        freematrix(&(marg.Vsosii2) );
        freematrix(&(marg.Vsosii3) );
        freematrix(&(marg.Vsosjj1) );
        freematrix(&(marg.Vsosjj2) );
        freematrix(&(marg.Vsosjj3) );
        freematrix(&(marg.pogeG1)  );
        freematrix(&(marg.pogeG2)  );
        freematrix(&(marg.pogeG3)  );
        freematrix(&(marg.pcont1)  );
        freematrix(&(marg.pcont2)  );
        freematrix(&(marg.pcont3)  );
        freematrix(&(marg.ptens1)  );
        freematrix(&(marg.ptens2)  );
        freematrix(&(marg.ptens3)  );
        freematrix(&(marg.psovii1) );
        freematrix(&(marg.psovii2) );
        freematrix(&(marg.psovii3) );
        freematrix(&(marg.psovjj1) );
        freematrix(&(marg.psovjj2) );
        freematrix(&(marg.psovjj3) );
        freematrix(&(marg.psovji1) );
        freematrix(&(marg.psovji2) );
        freematrix(&(marg.psovji3) );
        freematrix(&(marg.psovij1) );
        freematrix(&(marg.psovij2) );
        freematrix(&(marg.psovij3) );
        freematrix(&(marg.psosii1) );
        freematrix(&(marg.psosii2) );
        freematrix(&(marg.psosii3) );
        freematrix(&(marg.psosjj1) );
        freematrix(&(marg.psosjj2) );
        freematrix(&(marg.psosjj3) );
        freematrix(&(marg.T1)      );
        freematrix(&(marg.T2)      );
        freematrix(&(marg.T3)      );
        freematrix(&(marg.rmsr12)  );
        freematrix(&(marg.rmsr13)  );
        freematrix(&(marg.rmsr23)  );
        freematrix(&(marg.rmsl12)  );
        freematrix(&(marg.rmsl13)  );
        freematrix(&(marg.rmsl23)  );
	
	freematrix(&(marg.Hfi)     );
	
	freematrix(&(marg.u)       );
	
	freematrix(&(marg.v)       );

	free(marg.e1);
	free(marg.e2);

}

void debug02(int q1,int q2,int q3)
{ 

 //       debug01(q1,q2,q3,+1,1.5,-1,1,2);
 
	debug01(q1,q2,q3,-1,0.5,+1,0,0);
        debug01(q1,q2,q3,-1,0.5,-1,1,1);
        debug01(q1,q2,q3,-1,1.5,-1,1,1);
        debug01(q1,q2,q3,-1,1.5,+1,2,2);
        debug01(q1,q2,q3,-1,2.5,+1,2,2);
        debug01(q1,q2,q3,-1,2.5,-1,3,3);
        debug01(q1,q2,q3,-1,3.5,-1,3,3);

        debug01(q1,q2,q3,+1,0.5,+1,0,1);
        debug01(q1,q2,q3,+1,1.5,+1,0,1);
	debug01(q1,q2,q3,+1,0.5,-1,1,0);
        debug01(q1,q2,q3,+1,0.5,-1,1,1);
        debug01(q1,q2,q3,+1,1.5,-1,1,1);
        debug01(q1,q2,q3,+1,1.5,-1,1,2);
        debug01(q1,q2,q3,+1,2.5,-1,1,2);
        debug01(q1,q2,q3,+1,0.5,+1,2,1);
        debug01(q1,q2,q3,+1,1.5,+1,2,1);
        debug01(q1,q2,q3,+1,1.5,+1,2,2);
        debug01(q1,q2,q3,+1,2.5,+1,2,2);
        debug01(q1,q2,q3,+1,2.5,+1,2,3);
        debug01(q1,q2,q3,+1,3.5,+1,2,3);
        debug01(q1,q2,q3,+1,1.5,-1,3,2);
        debug01(q1,q2,q3,+1,2.5,-1,3,2);
        debug01(q1,q2,q3,+1,2.5,-1,3,3);
        debug01(q1,q2,q3,+1,3.5,-1,3,3);
        debug01(q1,q2,q3,+1,3.5,-1,3,4);
        debug01(q1,q2,q3,+1,4.5,-1,3,4);

}

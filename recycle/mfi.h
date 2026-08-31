typedef struct{
        basis_list qnlist_spfy;
        basis_list qnlist_full;
        matrix **mlsj;
        sumckdk_scdk ****scdk_cent_1;
        sumckdk_scdk ****scdk_cent_2;
        sumckdk_scdk ****scdk_cent_3;
        sumckdk_scdk ****scdk_cont_1;
        sumckdk_scdk ****scdk_cont_2;
        sumckdk_scdk ****scdk_cont_3;
        sumckdk_scdk ****scdk_tens_1;
        sumckdk_scdk ****scdk_tens_2;
        sumckdk_scdk ****scdk_tens_3;
        sumckdk_scdk ****scdk_soii_1;
        sumckdk_scdk ****scdk_soii_2;
        sumckdk_scdk ****scdk_soii_3;
        sumckdk_scdk ****scdk_soij_1;
        sumckdk_scdk ****scdk_soij_2;
        sumckdk_scdk ****scdk_soij_3;
        sumckdk_scdk ****scdk_soji_1;
        sumckdk_scdk ****scdk_soji_2;
        sumckdk_scdk ****scdk_soji_3;
        sumckdk_scdk ****scdk_sojj_1;
        sumckdk_scdk ****scdk_sojj_2;
        sumckdk_scdk ****scdk_sojj_3;
	vargs varg;
        
	matrix Nfi;
	matrix VogeG1;
	matrix VogeG2;
	matrix VogeG3;
	matrix Vcont1;
	matrix Vcont2;
	matrix Vcont3;
	matrix Vtens1;
	matrix Vtens2;
	matrix Vtens3;
	matrix Vsovii1;
	matrix Vsovii2;
	matrix Vsovii3;
	matrix Vsovjj1;
	matrix Vsovjj2;
	matrix Vsovjj3;
	matrix Vsovji1;
	matrix Vsovji2;
	matrix Vsovji3;
	matrix Vsovij1;
	matrix Vsovij2;
	matrix Vsovij3;
	matrix Vstring1;
	matrix Vstring2;
	matrix Vstring3;
	matrix Vsosii1;
	matrix Vsosii2;
	matrix Vsosii3;
	matrix Vsosjj1;
	matrix Vsosjj2;
	matrix Vsosjj3;
	matrix pogeG1;
	matrix pogeG2;
	matrix pogeG3;
	matrix pcont1;
	matrix pcont2;
	matrix pcont3;
	matrix ptens1;
	matrix ptens2;
	matrix ptens3;
	matrix psovii1;
	matrix psovii2;
	matrix psovii3;
	matrix psovjj1;
	matrix psovjj2;
	matrix psovjj3;
	matrix psovji1;
	matrix psovji2;
	matrix psovji3;
	matrix psovij1;
	matrix psovij2;
	matrix psovij3;
	matrix psosii1;
	matrix psosii2;
	matrix psosii3;
	matrix psosjj1;
	matrix psosjj2;
	matrix psosjj3;
	matrix T1;
	matrix T2;
	matrix T3;
	matrix rmsr12;
	matrix rmsr13;
	matrix rmsr23;
	matrix rmsl12;
	matrix rmsl13;
	matrix rmsl23;

	matrix Hfi;
	
	matrix u;
		
	matrix v;
	double *e1;
	double *e2;
}margs;




void malloc_scdk(int *len_part,int len_list,sumckdk_scdk *****scdk)
{
	int i,j,k;
	*scdk=(sumckdk_scdk****)malloc(sizeof(sumckdk_scdk***)*len_list);
	for(i=0;i<len_list;i++)
	{
		(*scdk)[i]=(sumckdk_scdk***)malloc(sizeof(sumckdk_scdk**)*len_part[i]);
		for(j=0;j<len_part[i];j++)
		{
			(*scdk)[i][j]=(sumckdk_scdk**)malloc(sizeof(sumckdk_scdk*)*len_list);
			for(k=0;k<len_list;k++)
			{
				(*scdk)[i][j][k]=(sumckdk_scdk*)malloc(sizeof(sumckdk_scdk)*len_part[k]);
			}
		}
	}
}

void free_scdk(int *len_part,int len_list,sumckdk_scdk *****scdk)
{
	int i,j,k,l;
	for(i=0;i<len_list;i++)
	{
		for(j=0;j<len_part[i];j++)
		{
			for(k=0;k<len_list;k++)
			{
				for(l=0;l<len_part[k];l++)
				{
					sumckdk_scdk_free(&((*scdk)[i][j][k][l]));
				}
				free((*scdk)[i][j][k]);
			}
			free((*scdk)[i][j]);
		}
		free((*scdk)[i]);
	}
	free(*scdk);
}




void* calc_scdk(void *args)
{
	mt_args *mtargs = (mt_args*)args;
	int ith=mtargs->ith;
	pthread_mutex_t *mutex=mtargs->mutex;
	int lock=mtargs->lock;
	margs *arg = (margs*)mtargs->p;
	int  len_list=arg->qnlist_spfy.len_list;
	int *len_part=arg->qnlist_spfy.len_part;
	int nf=ith;
	int nfp,ni,nip;

	for(nfp=0;nfp<len_part[nf];nfp++)
	{
		for(ni=0;ni<len_list;ni++)
		{
			for(nip=0;nip<len_part[ni];nip++)
			{
				sumckdk_scdk_vtype(&(arg->scdk_cent_1[nf][nfp][ni][nip]),vcent,1,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,1);
				sumckdk_scdk_vtype(&(arg->scdk_cent_2[nf][nfp][ni][nip]),vcent,1,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,2);
				sumckdk_scdk_vtype(&(arg->scdk_cent_3[nf][nfp][ni][nip]),vcent,1,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,3);
				
				sumckdk_scdk_vtype(&(arg->scdk_cont_1[nf][nfp][ni][nip]),vcont,1,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,1);
				sumckdk_scdk_vtype(&(arg->scdk_cont_2[nf][nfp][ni][nip]),vcont,1,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,2);
				sumckdk_scdk_vtype(&(arg->scdk_cont_3[nf][nfp][ni][nip]),vcont,1,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,3);
				
				sumckdk_scdk_vtype(&(arg->scdk_tens_1[nf][nfp][ni][nip]),vtens,2,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,1);
				sumckdk_scdk_vtype(&(arg->scdk_tens_2[nf][nfp][ni][nip]),vtens,2,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,2);
				sumckdk_scdk_vtype(&(arg->scdk_tens_3[nf][nfp][ni][nip]),vtens,2,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,3);
			
				sumckdk_scdk_vtype(&(arg->scdk_soii_1[nf][nfp][ni][nip]),vsoii,3,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,1);
				sumckdk_scdk_vtype(&(arg->scdk_soii_2[nf][nfp][ni][nip]),vsoii,3,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,2);
				sumckdk_scdk_vtype(&(arg->scdk_soii_3[nf][nfp][ni][nip]),vsoii,3,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,3);
			
				sumckdk_scdk_vtype(&(arg->scdk_sojj_1[nf][nfp][ni][nip]),vsojj,4,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,1);
				sumckdk_scdk_vtype(&(arg->scdk_sojj_2[nf][nfp][ni][nip]),vsojj,4,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,2);
				sumckdk_scdk_vtype(&(arg->scdk_sojj_3[nf][nfp][ni][nip]),vsojj,4,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,3);
			
				sumckdk_scdk_vtype(&(arg->scdk_soji_1[nf][nfp][ni][nip]),vsoji,5,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,1);
				sumckdk_scdk_vtype(&(arg->scdk_soji_2[nf][nfp][ni][nip]),vsoji,5,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,2);
				sumckdk_scdk_vtype(&(arg->scdk_soji_3[nf][nfp][ni][nip]),vsoji,5,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,3);
			
				sumckdk_scdk_vtype(&(arg->scdk_soij_1[nf][nfp][ni][nip]),vsoij,6,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,1);
				sumckdk_scdk_vtype(&(arg->scdk_soij_2[nf][nfp][ni][nip]),vsoij,6,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,2);
				sumckdk_scdk_vtype(&(arg->scdk_soij_3[nf][nfp][ni][nip]),vsoij,6,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,3);
			}
		}
	}
	printf("ith=%d\n",ith);
	return NULL;
}


void* calc_scdk_mt(void *args)
{
	mt_args *mtargs = (mt_args*)args;
	int ith=mtargs->ith;
	pthread_mutex_t *mutex=mtargs->mutex;
	int lock=mtargs->lock;
	margs *arg = (margs*)mtargs->p;
	int  len_list=arg->qnlist_spfy.len_list;
	int *len_part=arg->qnlist_spfy.len_part;
	int nf=ith/21;
	int nfp,ni,nip;

	for(nfp=0;nfp<len_part[nf];nfp++)
	{
		for(ni=0;ni<len_list;ni++)
		{
			for(nip=0;nip<len_part[ni];nip++)
			{
				switch(ith%21)
				{
					case  0: sumckdk_scdk_vtype(&(arg->scdk_cent_1[nf][nfp][ni][nip]),vcent,1,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,1);break;
					case  1: sumckdk_scdk_vtype(&(arg->scdk_cent_2[nf][nfp][ni][nip]),vcent,1,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,2);break;
					case  2: sumckdk_scdk_vtype(&(arg->scdk_cent_3[nf][nfp][ni][nip]),vcent,1,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,3);break;
				
					case  3: sumckdk_scdk_vtype(&(arg->scdk_cont_1[nf][nfp][ni][nip]),vcont,1,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,1);break;
					case  4: sumckdk_scdk_vtype(&(arg->scdk_cont_2[nf][nfp][ni][nip]),vcont,1,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,2);break;
					case  5: sumckdk_scdk_vtype(&(arg->scdk_cont_3[nf][nfp][ni][nip]),vcont,1,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,3);break;
				
					case  6: sumckdk_scdk_vtype(&(arg->scdk_tens_1[nf][nfp][ni][nip]),vtens,2,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,1);break;
					case  7: sumckdk_scdk_vtype(&(arg->scdk_tens_2[nf][nfp][ni][nip]),vtens,2,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,2);break;
					case  8: sumckdk_scdk_vtype(&(arg->scdk_tens_3[nf][nfp][ni][nip]),vtens,2,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,3);break;
				
					case  9: sumckdk_scdk_vtype(&(arg->scdk_soii_1[nf][nfp][ni][nip]),vsoii,3,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,1);break;
					case 10: sumckdk_scdk_vtype(&(arg->scdk_soii_2[nf][nfp][ni][nip]),vsoii,3,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,2);break;
					case 11: sumckdk_scdk_vtype(&(arg->scdk_soii_3[nf][nfp][ni][nip]),vsoii,3,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,3);break;
				
					case 12: sumckdk_scdk_vtype(&(arg->scdk_sojj_1[nf][nfp][ni][nip]),vsojj,4,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,1);break;
					case 13: sumckdk_scdk_vtype(&(arg->scdk_sojj_2[nf][nfp][ni][nip]),vsojj,4,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,2);break;
					case 14: sumckdk_scdk_vtype(&(arg->scdk_sojj_3[nf][nfp][ni][nip]),vsojj,4,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,3);break;
				
					case 15: sumckdk_scdk_vtype(&(arg->scdk_soji_1[nf][nfp][ni][nip]),vsoji,5,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,1);break;
					case 16: sumckdk_scdk_vtype(&(arg->scdk_soji_2[nf][nfp][ni][nip]),vsoji,5,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,2);break;
					case 17: sumckdk_scdk_vtype(&(arg->scdk_soji_3[nf][nfp][ni][nip]),vsoji,5,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,3);break;
				
					case 18: sumckdk_scdk_vtype(&(arg->scdk_soij_1[nf][nfp][ni][nip]),vsoij,6,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,1);break;
					case 19: sumckdk_scdk_vtype(&(arg->scdk_soij_2[nf][nfp][ni][nip]),vsoij,6,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,2);break;
					case 20: sumckdk_scdk_vtype(&(arg->scdk_soij_3[nf][nfp][ni][nip]),vsoij,6,arg->qnlist_spfy,arg->mlsj,nf,nfp,ni,nip,mutex,lock,3);break;
					
					default:printf("error_calc_scdk_mt\n");break;
				}
			}
		}
	}
	
	return NULL;
}


void basis_mlsj_jl(margs *marg)
{
	int i,j;
	int lrho,llam,L,c;
	double sij,jl,J,MJ;
	double s1,s2,s3,si,sj,sk,ms1,ms2,ms3,msi,msj,msk,mrho,mlam;
	double coe,cgf;
	basis_list *qnlist=&(marg->qnlist_spfy);
	for(i=0;i<qnlist->len_list;i++)
	{
		for(j=0;j<qnlist->len_part[i];j++)
		{
			initmatrix(&(marg->mlsj[i][j]),0,6);
			 coe=qnlist->qnum[i][j].coe;
			  s1=qnlist->qnum[i][j].s1;
			  s2=qnlist->qnum[i][j].s2;
			  s3=qnlist->qnum[i][j].s3;
			   c=qnlist->qnum[i][j].c;
			lrho=qnlist->qnum[i][j].lrho;
			llam=qnlist->qnum[i][j].llam;
		 	   L=qnlist->qnum[i][j].L;
			 sij=qnlist->qnum[i][j].sij;
			  jl=qnlist->qnum[i][j].jl;
			   J=qnlist->qnum[i][j].J;
			  MJ=J;

			getijk(s1,s2,s3,&si,&sj,&sk,c);
			for(msi=-si;msi<=si;msi++)
			{
				for(msj=-sj;msj<=sj;msj++)
				{
					for(msk=-sk;msk<=sk;msk++)
					{
						for(mrho=-lrho;mrho<=lrho;mrho++)
						{
							for(mlam=-llam;mlam<=llam;mlam++)
							{
								cgf=coe*clebschGordan(si,msi,sj,msj,sij,msi+msj)
								       *clebschGordan(lrho,mrho,llam,mlam,L,mrho+mlam)
								       *clebschGordan(sij,msi+msj,L,mrho+mlam,jl,msi+msj+mrho+mlam)
								       *clebschGordan(jl,msi+msj+mrho+mlam,sk,msk,J,MJ);
								if(0!=cgf)
								{
									get123(&ms1,&ms2,&ms3,msi,msj,msk,c);
									pushmatrix(&(marg->mlsj[i][j]));
									marg->mlsj[i][j].p[marg->mlsj[i][j].n-1][0]=cgf;
									marg->mlsj[i][j].p[marg->mlsj[i][j].n-1][1]=ms1;
									marg->mlsj[i][j].p[marg->mlsj[i][j].n-1][2]=ms2;
									marg->mlsj[i][j].p[marg->mlsj[i][j].n-1][3]=ms3;
									marg->mlsj[i][j].p[marg->mlsj[i][j].n-1][4]=mrho;	
								       	marg->mlsj[i][j].p[marg->mlsj[i][j].n-1][5]=mlam;
                                                                }
                                                        }
                                                }
                                        }
                                }
                        }
                }
        }
}



void getlsj_jl(margs *marg,double m1,double m2,double m3,double rmin,double rmax,int nmax,int f12,double J,int P,int Lmax,double jl)
{
        int i,lrho,llam,L;
        double si,sj,sk,sij,t1,t2,t3,tij,T;
	int c;

	basis_list_init(&(marg->qnlist_spfy));
	basis_list_init(&(marg->qnlist_full));

        si=0.5;
        sj=0.5;
        sk=0.5;
        t1=0;
        t2=0;
        t3=0;
        tij=0;
        T=0;

        for(lrho=0;lrho<=0;lrho++)
        {
                for(llam=Lmax;llam<=Lmax;llam++)
                {
                        for(L=abs(lrho-llam);L<=lrho+llam;L++)
                        {
                                for(sij=0;sij<=1;sij++)
                                {
					if(lrho+llam<=Lmax && P==pow(-1,lrho+llam) && fabs(L-sij)<=jl && jl<=L+sij && fabs(jl-sk)<=J && J<=jl+sk)
					{
						if(1==f12*pow(-1,1+sij+lrho))
						{
							c=1;basis_list_push(&(marg->qnlist_spfy),1,-1,-1,1.0,m1,m2,m3,si,sj,sk,t1,t2,t3,tij,T,c,lrho,llam,L,sij,jl,J,0,0,0,0);
						}
						//	c=2;basis_list_push(&(marg->qnlist_spfy),1,-1,-1,1.0,                   m1,m2,m3,si,sj,sk,t1,t2,t3,tij,T,c,lrho,llam,L,sij,jl,J,0,0,0,0);
						//	c=3;basis_list_push(&(marg->qnlist_spfy),0,-1,-1,f12*pow(-1,1+sij+lrho),m1,m2,m3,si,sj,sk,t1,t2,t3,tij,T,c,lrho,llam,L,sij,jl,J,0,0,0,0);
					}
                                }
                        }
                }
        }

	basis_list_push_full(&(marg->qnlist_spfy),&(marg->qnlist_full),rmin,rmax,nmax);

	marg->mlsj=(matrix**)malloc(sizeof(matrix*)*(marg->qnlist_spfy.len_list));
	for(i=0;i<marg->qnlist_spfy.len_list;i++)
	{
		marg->mlsj[i]=(matrix*)malloc(sizeof(matrix)*(marg->qnlist_spfy.len_part[i]));
	}
	basis_mlsj_jl(marg);
}




void *getmfi(void *args)
{
	mt_args *mtarg =(mt_args*)args;
	margs *arg =(margs*)mtarg->p;
	int nf=mtarg->ith;
	int nfp,ni,nip;
	int mapf1,mapf2,mapi1,mapi2;
	
	for(ni=0;ni<arg->qnlist_full.len_list;ni++)
	{
		for(nfp=0;nfp<arg->qnlist_full.len_part[nf];nfp++)
		{
			for(nip=0;nip<arg->qnlist_full.len_part[ni];nip++)
			{
				mapf1=arg->qnlist_full.qnum[nf][nfp].map1;
				mapf2=arg->qnlist_full.qnum[nf][nfp].map2;
				mapi1=arg->qnlist_full.qnum[ni][nip].map1;
				mapi2=arg->qnlist_full.qnum[ni][nip].map2;
				
				arg->Nfi   .   p[nf][ni]=inteVcenPartA(tir_cent,arg->scdk_cent_1[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteNfi,   1);

				arg->VogeG1.   p[nf][ni]=inteVcenPartA(tir_cent,arg->scdk_cent_1[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteVogeG, 1);
				arg->VogeG2.   p[nf][ni]=inteVcenPartA(tir_cent,arg->scdk_cent_2[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteVogeG, 2);
				arg->VogeG3.   p[nf][ni]=inteVcenPartA(tir_cent,arg->scdk_cent_3[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteVogeG, 3);

				arg->Vcont1.   p[nf][ni]=inteVcenPartA(tir_cent,arg->scdk_cont_1[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteVcont, 1);
				arg->Vcont2.   p[nf][ni]=inteVcenPartA(tir_cent,arg->scdk_cont_2[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteVcont, 2);
				arg->Vcont3.   p[nf][ni]=inteVcenPartA(tir_cent,arg->scdk_cont_3[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteVcont, 3);
				
				arg->Vtens1.  p[nf][ni]=inteVcenPartA(tir_tens,arg->scdk_tens_1[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteVtens,  1);
				arg->Vtens2.  p[nf][ni]=inteVcenPartA(tir_tens,arg->scdk_tens_2[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteVtens,  2);
				arg->Vtens3.  p[nf][ni]=inteVcenPartA(tir_tens,arg->scdk_tens_3[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteVtens,  3);
				
				arg->Vsovii1. p[nf][ni]=inteVcenPartA(tir_soii,arg->scdk_soii_1[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteVsovii, 1);
				arg->Vsovii2. p[nf][ni]=inteVcenPartA(tir_soii,arg->scdk_soii_2[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteVsovii, 2);
				arg->Vsovii3. p[nf][ni]=inteVcenPartA(tir_soii,arg->scdk_soii_3[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteVsovii, 3);
				
				arg->Vsovjj1. p[nf][ni]=inteVcenPartA(tjr_sojj,arg->scdk_sojj_1[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteVsovjj, 1);
				arg->Vsovjj2. p[nf][ni]=inteVcenPartA(tjr_sojj,arg->scdk_sojj_2[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteVsovjj, 2);
				arg->Vsovjj3. p[nf][ni]=inteVcenPartA(tjr_sojj,arg->scdk_sojj_3[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteVsovjj, 3);
				
				arg->Vsovji1. p[nf][ni]=inteVcenPartA(t1r_soji,arg->scdk_soji_1[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteVsovji, 1);
				arg->Vsovji2. p[nf][ni]=inteVcenPartA(t1r_soji,arg->scdk_soji_2[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteVsovji, 2);
				arg->Vsovji3. p[nf][ni]=inteVcenPartA(t1r_soji,arg->scdk_soji_3[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteVsovji, 3);
				
				arg->Vsovij1. p[nf][ni]=inteVcenPartA(tir_soij,arg->scdk_soij_1[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteVsovij, 1);
				arg->Vsovij2. p[nf][ni]=inteVcenPartA(tir_soij,arg->scdk_soij_2[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteVsovij, 2);
				arg->Vsovij3. p[nf][ni]=inteVcenPartA(tir_soij,arg->scdk_soij_3[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteVsovij, 3);
				
				arg->Vstring1.p[nf][ni]=inteVcenPartA(tir_cent,arg->scdk_cent_1[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteVstring,1);
				arg->Vstring2.p[nf][ni]=inteVcenPartA(tir_cent,arg->scdk_cent_2[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteVstring,2);
				arg->Vstring3.p[nf][ni]=inteVcenPartA(tir_cent,arg->scdk_cent_3[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteVstring,3);
				
				arg->Vsosii1. p[nf][ni]=inteVcenPartA(tir_soii,arg->scdk_soii_1[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteVsosii, 1);
				arg->Vsosii2. p[nf][ni]=inteVcenPartA(tir_soii,arg->scdk_soii_2[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteVsosii, 2);
				arg->Vsosii3. p[nf][ni]=inteVcenPartA(tir_soii,arg->scdk_soii_3[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteVsosii, 3);
				
				arg->Vsosjj1. p[nf][ni]=inteVcenPartA(tjr_sojj,arg->scdk_sojj_1[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteVsosjj, 1);
				arg->Vsosjj2. p[nf][ni]=inteVcenPartA(tjr_sojj,arg->scdk_sojj_2[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteVsosjj, 2);
				arg->Vsosjj3. p[nf][ni]=inteVcenPartA(tjr_sojj,arg->scdk_sojj_3[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteVsosjj, 3);
                                
				arg->pogeG1.  p[nf][ni]=inteVcenPartA(t1p_cent,arg->scdk_cent_1[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,intepogeG,  1);
                                arg->pogeG2.  p[nf][ni]=inteVcenPartA(t1p_cent,arg->scdk_cent_2[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,intepogeG,  2);
                                arg->pogeG3.  p[nf][ni]=inteVcenPartA(t1p_cent,arg->scdk_cent_3[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,intepogeG,  3);

                                arg->pcont1.  p[nf][ni]=inteVcenPartA(t1p_cent,arg->scdk_cent_1[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,intepcont,  1);
                                arg->pcont2.  p[nf][ni]=inteVcenPartA(t1p_cent,arg->scdk_cent_2[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,intepcont,  2);
                                arg->pcont3.  p[nf][ni]=inteVcenPartA(t1p_cent,arg->scdk_cent_3[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,intepcont,  3);

                                arg->ptens1.  p[nf][ni]=inteVcenPartA(t1p_cent,arg->scdk_cent_1[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteptens,  1);
                                arg->ptens2.  p[nf][ni]=inteVcenPartA(t1p_cent,arg->scdk_cent_2[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteptens,  2);
                                arg->ptens3.  p[nf][ni]=inteVcenPartA(t1p_cent,arg->scdk_cent_3[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteptens,  3);

                                arg->psovii1. p[nf][ni]=inteVcenPartA(t1p_cent,arg->scdk_cent_1[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,intepsovii, 1);
                                arg->psovii2. p[nf][ni]=inteVcenPartA(t1p_cent,arg->scdk_cent_2[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,intepsovii, 2);
                                arg->psovii3. p[nf][ni]=inteVcenPartA(t1p_cent,arg->scdk_cent_3[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,intepsovii, 3);

                                arg->psovjj1. p[nf][ni]=inteVcenPartA(t1p_cent,arg->scdk_cent_1[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,intepsovjj, 1);
                                arg->psovjj2. p[nf][ni]=inteVcenPartA(t1p_cent,arg->scdk_cent_2[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,intepsovjj, 2);
                                arg->psovjj3. p[nf][ni]=inteVcenPartA(t1p_cent,arg->scdk_cent_3[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,intepsovjj, 3);

                                arg->psovji1. p[nf][ni]=inteVcenPartA(t1p_cent,arg->scdk_cent_1[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,intepsovji, 1);
                                arg->psovji2. p[nf][ni]=inteVcenPartA(t1p_cent,arg->scdk_cent_2[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,intepsovji, 2);
                                arg->psovji3. p[nf][ni]=inteVcenPartA(t1p_cent,arg->scdk_cent_3[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,intepsovji, 3);

                                arg->psovij1. p[nf][ni]=inteVcenPartA(t1p_cent,arg->scdk_cent_1[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,intepsovij, 1);
                                arg->psovij2. p[nf][ni]=inteVcenPartA(t1p_cent,arg->scdk_cent_2[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,intepsovij, 2);
                                arg->psovij3. p[nf][ni]=inteVcenPartA(t1p_cent,arg->scdk_cent_3[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,intepsovij, 3);

                                arg->psosii1. p[nf][ni]=inteVcenPartA(t1p_cent,arg->scdk_cent_1[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,intepsosii, 1);
                                arg->psosii2. p[nf][ni]=inteVcenPartA(t1p_cent,arg->scdk_cent_2[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,intepsosii, 2);
                                arg->psosii3. p[nf][ni]=inteVcenPartA(t1p_cent,arg->scdk_cent_3[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,intepsosii, 3);

                                arg->psosjj1. p[nf][ni]=inteVcenPartA(t1p_cent,arg->scdk_cent_1[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,intepsosjj, 1);
                                arg->psosjj2. p[nf][ni]=inteVcenPartA(t1p_cent,arg->scdk_cent_2[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,intepsosjj, 2);
                                arg->psosjj3. p[nf][ni]=inteVcenPartA(t1p_cent,arg->scdk_cent_3[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,intepsosjj, 3);

				arg->T1.      p[nf][ni]=inteVcenPartA(tpi_cent,arg->scdk_cent_1[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteTi,     1);
				arg->T2.      p[nf][ni]=inteVcenPartA(tpi_cent,arg->scdk_cent_2[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteTi,     2);
				arg->T3.      p[nf][ni]=inteVcenPartA(tpi_cent,arg->scdk_cent_3[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteTi,     3);
				
				arg->rmsr12.  p[nf][ni]=inteVcenPartA(tir_cent,arg->scdk_cent_1[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteRMS,    1);
				arg->rmsr13.  p[nf][ni]=inteVcenPartA(tir_cent,arg->scdk_cent_2[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteRMS,    2);
				arg->rmsr23.  p[nf][ni]=inteVcenPartA(tir_cent,arg->scdk_cent_3[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteRMS,    3);
				
				arg->rmsl12.  p[nf][ni]=inteVcenPartA(t2r_cent,arg->scdk_cent_1[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteRMS,    1);
				arg->rmsl13.  p[nf][ni]=inteVcenPartA(t2r_cent,arg->scdk_cent_2[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteRMS,    2);
				arg->rmsl23.  p[nf][ni]=inteVcenPartA(t2r_cent,arg->scdk_cent_3[mapf1][mapf2][mapi1][mapi2],arg->qnlist_full.qnum[nf][nfp],arg->qnlist_full.qnum[ni][nip],arg->varg,inteRMS,    3);


				}
		}
	}
	return NULL;
}


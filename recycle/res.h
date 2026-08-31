void matmuln(matrix c,matrix a,matrix b)
{
	int i,j,k;
	for(i=0;i<c.n;i++)
	{
		for(j=0;j<c.m;j++)
		{
			c.p[i][j]=0;
		}
	}
	for(i=0;i<a.n;i++)
	{
		for(k=0;k<a.m;k++)
		{
			for(j=0;j<b.m;j++)
			{
				c.p[i][j]+=a.p[i][k]*b.p[k][j];
			}
		}
	}
}

void matmult(matrix c,matrix a,matrix b)
{
	int i,j,k;
	for(i=0;i<c.n;i++)
	{
		for(j=0;j<c.m;j++)
		{
			c.p[i][j]=0;
		}
	}
	for(i=0;i<a.n;i++)
	{
		for(j=0;j<b.n;j++)
		{
			for(k=0;k<b.m;k++)
			{
				c.p[i][j]+=a.p[i][k]*b.p[j][k];
			}
		}
	}
}
void transU(matrix a,matrix u)
{
	int n=a.n;
	matrix temp;
	initmatrix(&temp,n,n);
	matmuln(temp,u,a);
	matmult(a,temp,u);
	freematrix(&temp);
}

void transP(matrix a,matrix p)
{
	int n=a.n;
	matrix temp;
	initmatrix(&temp,n,n);
	matmult(temp,p,a);
	matmult(a,temp,p);
	freematrix(&temp);
}

void matrixadd(matrix a,matrix b)
{
	int i,j;
	for(i=0;i<a.n;i++)
	{
		for(j=0;j<b.m;j++)
		{
			a.p[i][j]+=b.p[i][j];
		}
	}
}

void transUm(pthread_mutex_t *mutex,int lock,matrix a,matrix u)
{
	int n=a.n;
	matrix temp;
	if(1==lock)
	{
		pthread_mutex_lock(mutex);
	}
	initmatrix(&temp,n,n);
	if(1==lock)
	{
		pthread_mutex_unlock(mutex);
	}
	matmuln(temp,u,a);
	matmult(a,temp,u);
	if(1==lock)
	{
		pthread_mutex_lock(mutex);
	}
	freematrix(&temp);
	if(1==lock)
	{
		pthread_mutex_unlock(mutex);
	}
}

void transPm(pthread_mutex_t *mutex,int lock,matrix a,matrix p)
{
	int n=a.n;
	matrix temp;
	if(1==lock)
	{
		pthread_mutex_lock(mutex);
	}
	initmatrix(&temp,n,n);
	if(1==lock)
	{
		pthread_mutex_unlock(mutex);
	}
	matmult(temp,p,a);
	matmult(a,temp,p);
	if(1==lock)
	{
		pthread_mutex_lock(mutex);
	}
	freematrix(&temp);
	if(1==lock)
	{
		pthread_mutex_unlock(mutex);
	}
}

void *transUP(void *args)
{
	
	mt_args *mtarg =(mt_args*)args;
	int t=mtarg->ith;
	pthread_mutex_t *mutex=mtarg->mutex;
	int lock=mtarg->lock;
	margs *arg=(margs*)mtarg->p;
	
	switch (t)
	{
		case 0:
		transUm(mutex,lock,arg->VogeG1,arg->u);
		transUm(mutex,lock,arg->pogeG1,arg->u);
		transPm(mutex,lock,arg->VogeG1,arg->pogeG1);
		break;

		case 1:
		transUm(mutex,lock,arg->VogeG2,arg->u);
		transUm(mutex,lock,arg->pogeG2,arg->u);
		transPm(mutex,lock,arg->VogeG2,arg->pogeG2);
		break;

		case 2:
		transUm(mutex,lock,arg->VogeG3,arg->u);
		transUm(mutex,lock,arg->pogeG3,arg->u);
		transPm(mutex,lock,arg->VogeG3,arg->pogeG3);
		break;

		case 3:
		transUm(mutex,lock,arg->Vcont1,arg->u);
		transUm(mutex,lock,arg->pcont1,arg->u);
		transPm(mutex,lock,arg->Vcont1,arg->pcont1);
		break;

		case 4:
		transUm(mutex,lock,arg->Vcont2,arg->u);
		transUm(mutex,lock,arg->pcont2,arg->u);
		transPm(mutex,lock,arg->Vcont2,arg->pcont2);
		break;

		case 5:
		transUm(mutex,lock,arg->Vcont3,arg->u);
		transUm(mutex,lock,arg->pcont3,arg->u);
		transPm(mutex,lock,arg->Vcont3,arg->pcont3);
		break;

		case 6:
		transUm(mutex,lock,arg->Vtens1,arg->u);
		transUm(mutex,lock,arg->ptens1,arg->u);
		transPm(mutex,lock,arg->Vtens1,arg->ptens1);
		break;

		case 7:
		transUm(mutex,lock,arg->Vtens2,arg->u);
		transUm(mutex,lock,arg->ptens2,arg->u);
		transPm(mutex,lock,arg->Vtens2,arg->ptens2);
		break;

		case 8:
		transUm(mutex,lock,arg->Vtens3,arg->u);
		transUm(mutex,lock,arg->ptens3,arg->u);
		transPm(mutex,lock,arg->Vtens3,arg->ptens3);
		break;

		case 9:
		transUm(mutex,lock,arg->Vsovii1,arg->u);
		transUm(mutex,lock,arg->psovii1,arg->u);
		transPm(mutex,lock,arg->Vsovii1,arg->psovii1);
		break;

		case 10:
		transUm(mutex,lock,arg->Vsovii2,arg->u);
		transUm(mutex,lock,arg->psovii2,arg->u);
		transPm(mutex,lock,arg->Vsovii2,arg->psovii2);
		break;

		case 11:
		transUm(mutex,lock,arg->Vsovii3,arg->u);
		transUm(mutex,lock,arg->psovii3,arg->u);
		transPm(mutex,lock,arg->Vsovii3,arg->psovii3);
		break;

		case 12:
		transUm(mutex,lock,arg->Vsovjj1,arg->u);
		transUm(mutex,lock,arg->psovjj1,arg->u);
		transPm(mutex,lock,arg->Vsovjj1,arg->psovjj1);
		break;

		case 13:
		transUm(mutex,lock,arg->Vsovjj2,arg->u);
		transUm(mutex,lock,arg->psovjj2,arg->u);
		transPm(mutex,lock,arg->Vsovjj2,arg->psovjj2);
		break;

		case 14:
		transUm(mutex,lock,arg->Vsovjj3,arg->u);
		transUm(mutex,lock,arg->psovjj3,arg->u);
		transPm(mutex,lock,arg->Vsovjj3,arg->psovjj3);
		break;

		case 15:
		transUm(mutex,lock,arg->Vsovji1,arg->u);
		transUm(mutex,lock,arg->psovji1,arg->u);
		transPm(mutex,lock,arg->Vsovji1,arg->psovji1);
		break;

		case 16:
		transUm(mutex,lock,arg->Vsovji2,arg->u);
		transUm(mutex,lock,arg->psovji2,arg->u);
		transPm(mutex,lock,arg->Vsovji2,arg->psovji2);
		break;

		case 17:
		transUm(mutex,lock,arg->Vsovji3,arg->u);
		transUm(mutex,lock,arg->psovji3,arg->u);
		transPm(mutex,lock,arg->Vsovji3,arg->psovji3);
		break;

		case 18:
		transUm(mutex,lock,arg->Vsovij1,arg->u);
		transUm(mutex,lock,arg->psovij1,arg->u);
		transPm(mutex,lock,arg->Vsovij1,arg->psovij1);
		break;

		case 19:
		transUm(mutex,lock,arg->Vsovij2,arg->u);
		transUm(mutex,lock,arg->psovij2,arg->u);
		transPm(mutex,lock,arg->Vsovij2,arg->psovij2);
		break;

		case 20:
		transUm(mutex,lock,arg->Vsovij3,arg->u);
		transUm(mutex,lock,arg->psovij3,arg->u);
		transPm(mutex,lock,arg->Vsovij3,arg->psovij3);
		break;

		case 21:
		transUm(mutex,lock,arg->Vstring1,arg->u);
		transUm(mutex,lock,arg->Vstring2,arg->u);
		transUm(mutex,lock,arg->Vstring3,arg->u);
		break;

		case 22:
		transUm(mutex,lock,arg->Vsosii1,arg->u);
		transUm(mutex,lock,arg->psosii1,arg->u);
		transPm(mutex,lock,arg->Vsosii1,arg->psosii1);
		break;

		case 23:
		transUm(mutex,lock,arg->Vsosii2,arg->u);
		transUm(mutex,lock,arg->psosii2,arg->u);
		transPm(mutex,lock,arg->Vsosii2,arg->psosii2);
		break;

		case 24:
		transUm(mutex,lock,arg->Vsosii3,arg->u);
		transUm(mutex,lock,arg->psosii3,arg->u);
		transPm(mutex,lock,arg->Vsosii3,arg->psosii3);
		break;

		case 25:
		transUm(mutex,lock,arg->Vsosjj1,arg->u);
		transUm(mutex,lock,arg->psosjj1,arg->u);
		transPm(mutex,lock,arg->Vsosjj1,arg->psosjj1);
		break;

		case 26:
		transUm(mutex,lock,arg->Vsosjj2,arg->u);
		transUm(mutex,lock,arg->psosjj2,arg->u);
		transPm(mutex,lock,arg->Vsosjj2,arg->psosjj2);
		break;

		case 27:
		transUm(mutex,lock,arg->Vsosjj3,arg->u);
		transUm(mutex,lock,arg->psosjj3,arg->u);
		transPm(mutex,lock,arg->Vsosjj3,arg->psosjj3);
		break;

		case 28:
		transUm(mutex,lock,arg->T1,arg->u);
		transUm(mutex,lock,arg->T2,arg->u);
		transUm(mutex,lock,arg->T3,arg->u);
		break;

		case 29:
		transUm(mutex,lock,arg->rmsr12,arg->u);
		transUm(mutex,lock,arg->rmsr13,arg->u);
		transUm(mutex,lock,arg->rmsr23,arg->u);
		break;

		case 30:
		transUm(mutex,lock,arg->rmsl12,arg->u);
		transUm(mutex,lock,arg->rmsl13,arg->u);
		transUm(mutex,lock,arg->rmsl23,arg->u);
		break;

		default:
		break;

	}
	return NULL;
}

void eigsys(margs *arg)
{
	int n=arg->qnlist_full.len_list;
	int info;
	int i,j,k;
	double rmsrho,rmslam;
	matrix temp;
	double *et1,*et2;
 
//        printmatrix(arg->Vsovjj1) ;
//        printmatrix(arg->Vsovjj2) ;
//        printmatrix(arg->Vsovjj3) ;
/*
        printmatrix(arg->Nfi)     ;
        printmatrix(arg->VogeG1)  ;
        printmatrix(arg->VogeG2)  ;
        printmatrix(arg->VogeG3)  ;
        printmatrix(arg->Vcont1)  ;
        printmatrix(arg->Vcont2)  ;
        printmatrix(arg->Vcont3)  ;
        printmatrix(arg->Vtens1)  ;
        printmatrix(arg->Vtens2)  ;
        printmatrix(arg->Vtens3)  ;
        printmatrix(arg->Vsovii1) ;
        printmatrix(arg->Vsovii2) ;
        printmatrix(arg->Vsovii3) ;
        printmatrix(arg->Vsovjj1) ;
        printmatrix(arg->Vsovjj2) ;
        printmatrix(arg->Vsovjj3) ;
        printmatrix(arg->Vsovji1) ;
        printmatrix(arg->Vsovji2) ;
        printmatrix(arg->Vsovji3) ;
        printmatrix(arg->Vsovij1) ;
        printmatrix(arg->Vsovij2) ;
        printmatrix(arg->Vsovij3) ;
        printmatrix(arg->Vstring1);
        printmatrix(arg->Vstring2);
        printmatrix(arg->Vstring3);
        printmatrix(arg->Vsosii1) ;
        printmatrix(arg->Vsosii2) ;
        printmatrix(arg->Vsosii3) ;
        printmatrix(arg->Vsosjj1) ;
        printmatrix(arg->Vsosjj2) ;
        printmatrix(arg->Vsosjj3) ;
        printmatrix(arg->pogeG1)  ;
        printmatrix(arg->pogeG2)  ;
        printmatrix(arg->pogeG3)  ;
        printmatrix(arg->pcont1)  ;
        printmatrix(arg->pcont2)  ;
        printmatrix(arg->pcont3)  ;
        printmatrix(arg->ptens1)  ;
        printmatrix(arg->ptens2)  ;
        printmatrix(arg->ptens3)  ;
        printmatrix(arg->psovii1) ;
        printmatrix(arg->psovii2) ;
        printmatrix(arg->psovii3) ;
        printmatrix(arg->psovjj1) ;
        printmatrix(arg->psovjj2) ;
        printmatrix(arg->psovjj3) ;
        printmatrix(arg->psovji1) ;
        printmatrix(arg->psovji2) ;
        printmatrix(arg->psovji3) ;
        printmatrix(arg->psovij1) ;
        printmatrix(arg->psovij2) ;
        printmatrix(arg->psovij3) ;
        printmatrix(arg->psosii1) ;
        printmatrix(arg->psosii2) ;
        printmatrix(arg->psosii3) ;
        printmatrix(arg->psosjj1) ;
        printmatrix(arg->psosjj2) ;
        printmatrix(arg->psosjj3) ;
        printmatrix(arg->T1)      ;
        printmatrix(arg->T2)      ;
        printmatrix(arg->T3)      ;
        printmatrix(arg->rmsr12)  ;
        printmatrix(arg->rmsr13)  ;
        printmatrix(arg->rmsr23)  ;
        printmatrix(arg->rmsl12)  ;
        printmatrix(arg->rmsl13)  ;
        printmatrix(arg->rmsl23)  ;
*/
	
	initmatrix(&temp,n,n);
	et1=(double*)malloc(sizeof(double)*n);
	et2=(double*)malloc(sizeof(double)*n);
	
	for(i=0;i<n;i++)
	{
		for(j=i;j<n;j++)
		{
			temp.p[i][j]=sin(n+4321*i+1234*j+0.0);
			temp.p[j][i]=temp.p[i][j];
		}
	}

	eigv2Mul(temp.p,arg->Nfi.p,n,et1,et2,arg->u.p,n,&info);


	mt_load(31,transUP,arg,numberProcessors());
//	mt_load(31,transUP,arg,2);

/*
        transU(arg->VogeG1,arg->u);
        transU(arg->pogeG1,arg->u);
        transP(arg->VogeG1,arg->pogeG1);

        transU(arg->VogeG2,arg->u);
        transU(arg->pogeG2,arg->u);
        transP(arg->VogeG2,arg->pogeG2);

        transU(arg->VogeG3,arg->u);
        transU(arg->pogeG3,arg->u);
        transP(arg->VogeG3,arg->pogeG3);

        transU(arg->Vcont1,arg->u);
        transU(arg->pcont1,arg->u);
        transP(arg->Vcont1,arg->pcont1);

        transU(arg->Vcont2,arg->u);
        transU(arg->pcont2,arg->u);
        transP(arg->Vcont2,arg->pcont2);

        transU(arg->Vcont3,arg->u);
        transU(arg->pcont3,arg->u);
        transP(arg->Vcont3,arg->pcont3);

        transU(arg->Vtens1,arg->u);
        transU(arg->ptens1,arg->u);
        transP(arg->Vtens1,arg->ptens1);

        transU(arg->Vtens2,arg->u);
        transU(arg->ptens2,arg->u);
        transP(arg->Vtens2,arg->ptens2);

        transU(arg->Vtens3,arg->u);
        transU(arg->ptens3,arg->u);
        transP(arg->Vtens3,arg->ptens3);

        transU(arg->Vsovii1,arg->u);
        transU(arg->psovii1,arg->u);
        transP(arg->Vsovii1,arg->psovii1);

        transU(arg->Vsovii2,arg->u);
        transU(arg->psovii2,arg->u);
        transP(arg->Vsovii2,arg->psovii2);

        transU(arg->Vsovii3,arg->u);
        transU(arg->psovii3,arg->u);
        transP(arg->Vsovii3,arg->psovii3);

        transU(arg->Vsovjj1,arg->u);
        transU(arg->psovjj1,arg->u);
        transP(arg->Vsovjj1,arg->psovjj1);

        transU(arg->Vsovjj2,arg->u);
        transU(arg->psovjj2,arg->u);
        transP(arg->Vsovjj2,arg->psovjj2);

        transU(arg->Vsovjj3,arg->u);
        transU(arg->psovjj3,arg->u);
        transP(arg->Vsovjj3,arg->psovjj3);

        transU(arg->Vsovji1,arg->u);
        transU(arg->psovji1,arg->u);
        transP(arg->Vsovji1,arg->psovji1);

        transU(arg->Vsovji2,arg->u);
        transU(arg->psovji2,arg->u);
        transP(arg->Vsovji2,arg->psovji2);

        transU(arg->Vsovji3,arg->u);
        transU(arg->psovji3,arg->u);
        transP(arg->Vsovji3,arg->psovji3);

        transU(arg->Vsovij1,arg->u);
        transU(arg->psovij1,arg->u);
        transP(arg->Vsovij1,arg->psovij1);

        transU(arg->Vsovij2,arg->u);
        transU(arg->psovij2,arg->u);
        transP(arg->Vsovij2,arg->psovij2);

        transU(arg->Vsovij3,arg->u);
        transU(arg->psovij3,arg->u);
        transP(arg->Vsovij3,arg->psovij3);

	transU(arg->Vstring1,arg->u);
	transU(arg->Vstring2,arg->u);
	transU(arg->Vstring3,arg->u);
        
	transU(arg->Vsosii1,arg->u);
        transU(arg->psosii1,arg->u);
        transP(arg->Vsosii1,arg->psosii1);

        transU(arg->Vsosii2,arg->u);
        transU(arg->psosii2,arg->u);
        transP(arg->Vsosii2,arg->psosii2);

        transU(arg->Vsosii3,arg->u);
        transU(arg->psosii3,arg->u);
        transP(arg->Vsosii3,arg->psosii3);

        transU(arg->Vsosjj1,arg->u);
        transU(arg->psosjj1,arg->u);
        transP(arg->Vsosjj1,arg->psosjj1);

        transU(arg->Vsosjj2,arg->u);
        transU(arg->psosjj2,arg->u);
        transP(arg->Vsosjj2,arg->psosjj2);

        transU(arg->Vsosjj3,arg->u);
        transU(arg->psosjj3,arg->u);
        transP(arg->Vsosjj3,arg->psosjj3);


	transU(arg->T1,arg->u);
	transU(arg->T2,arg->u);
	transU(arg->T3,arg->u);

	transU(arg->rmsr12,arg->u);
	transU(arg->rmsr13,arg->u);
	transU(arg->rmsr23,arg->u);

	transU(arg->rmsl12,arg->u);
	transU(arg->rmsl13,arg->u);
	transU(arg->rmsl23,arg->u);
*/
	matrixadd(arg->Hfi,arg->VogeG1);
	matrixadd(arg->Hfi,arg->VogeG2);
	matrixadd(arg->Hfi,arg->VogeG3);
	matrixadd(arg->Hfi,arg->Vcont1);
	matrixadd(arg->Hfi,arg->Vcont2);
	matrixadd(arg->Hfi,arg->Vcont3);
	matrixadd(arg->Hfi,arg->Vtens1);
	matrixadd(arg->Hfi,arg->Vtens2);
	matrixadd(arg->Hfi,arg->Vtens3);
	matrixadd(arg->Hfi,arg->Vsovii1);
	matrixadd(arg->Hfi,arg->Vsovii2);
	matrixadd(arg->Hfi,arg->Vsovii3);
	matrixadd(arg->Hfi,arg->Vsovjj1);
	matrixadd(arg->Hfi,arg->Vsovjj2);
	matrixadd(arg->Hfi,arg->Vsovjj3);
	matrixadd(arg->Hfi,arg->Vsovji1);
	matrixadd(arg->Hfi,arg->Vsovji2);
	matrixadd(arg->Hfi,arg->Vsovji3);
	matrixadd(arg->Hfi,arg->Vsovij1);
	matrixadd(arg->Hfi,arg->Vsovij2);
	matrixadd(arg->Hfi,arg->Vsovij3);
	matrixadd(arg->Hfi,arg->Vstring1);
	matrixadd(arg->Hfi,arg->Vstring2);
	matrixadd(arg->Hfi,arg->Vstring3);
	matrixadd(arg->Hfi,arg->Vsosii1);
	matrixadd(arg->Hfi,arg->Vsosii2);
	matrixadd(arg->Hfi,arg->Vsosii3);
	matrixadd(arg->Hfi,arg->Vsosjj1);
	matrixadd(arg->Hfi,arg->Vsosjj2);
	matrixadd(arg->Hfi,arg->Vsosjj3);
	matrixadd(arg->Hfi,arg->T1);
	matrixadd(arg->Hfi,arg->T2);
	matrixadd(arg->Hfi,arg->T3);

	eigv1Mul(arg->Hfi.p,n,arg->e1,arg->e2,arg->v.p,n);
	printarrayd1(arg->e1,3);

        for(i=0;i<3;i++)
        {
                rmsrho=0;
                rmslam=0;
                for(j=0;j<n;j++)
                {
                        for(k=0;k<n;k++)
                        {
                                rmsrho+=(arg->v.p[i][j])*(arg->rmsr12.p[j][k])*(arg->v.p[i][k]);
                                rmslam+=(arg->v.p[i][j])*(arg->rmsl12.p[j][k])*(arg->v.p[i][k]);
                        }
                }
//                printf("%d: %10.6f %10.6f %10.6f\n",i+1,arg->e1[i],1/sqrt(rmsrho),1/sqrt(rmslam));
        }
//	printf("\n");

	freematrix(&temp);
	free(et1);
	free(et2);
}

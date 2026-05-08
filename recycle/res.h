void eigsys(margs *arg)
{
	int n=arg->qnlist_full.len_list;
	double **Nfi=arg->Nfi.p;
	double **Hfi=arg->Hfi.p;
	double **v=arg->v.p;
	double *e1=arg->e1;
	double *e2=arg->e2;
	int info;
	/*
	printmatrix(arg->Nfi);
	printmatrix(arg->Hfi);*/
	eigv2Mul(Hfi,Nfi,n,e1,e2,v,n,&info);
	/*printf("info=%d\n",info);*/
	printarrayd1(e1,3);
}

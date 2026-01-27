#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
#include<ctime>
#include<cassert>
#include<pthread.h>


double chi2(const std::vector<double>& paras)
{
	const int len=10;
	double x[len]={-4.51978,-2.09609,1.2625,-4.20036,-4.73878,4.67549,-2.15036,-4.80579,-2.03644,4.19255};
	double res=0;
	int i;
	for(i=0;i<len;i++)
	{
		res=res+pow(paras[i]-x[i],2);
		printf("para_%02d=%10.6f x_%02d=%10.6f\n",i,paras[i],i,x[i]);
	}
	res=sqrt(res);
	printf("chi_2=%15.10f\n\n",res);
	getchar();
	return res;
}

#include "Minuit2/FCNBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnPrint.h"

class chi_2 : public ROOT::Minuit2::FCNBase {

public:
	chi_2();
	~chi_2() {}
	virtual double operator()(const std::vector<double>& ) const;
	virtual double Up() const { return ferrordef; }
	void SetErrorDef(double def) { ferrordef = def; }
private:
	double ferrordef;
	int Nparas;
	int Nerrs;
	int DOF;
};

chi_2::chi_2() : ferrordef(1.)
{
	Nparas=10;
	Nerrs=25;
	DOF=1;
}

bool openfile=true;

double chi_2::operator()(const std::vector<double>& paras) const {
	long double chi_square;
	assert(paras.size() == Nparas);
	chi_square=chi2(paras);
	return chi_square;
}

int main(void)
{
	srand(time(0)); 

	chi_2 MinuitFit;
	ROOT::Minuit2::MnUserParameters upar;

	double step = 5.0;
	double err=0.1;

	double xi = 0;

	upar.Add("x0",      xi,         err/step, -100, 100);
	upar.Add("x1",      xi,         err/step, -100, 100);
	upar.Add("x2",      xi,         err/step, -100, 100);
	upar.Add("x3",      xi,         err/step, -100, 100);
	upar.Add("x4",      xi,         err/step, -100, 100);
	upar.Add("x5",      xi,         err/step, -100, 100);
	upar.Add("x6",      xi,         err/step, -100, 100);
	upar.Add("x7",      xi,         err/step, -100, 100);
	upar.Add("x8",      xi,         err/step, -100, 100);
	upar.Add("x9",      xi,         err/step, -100, 100);

	ROOT::Minuit2::MnMigrad migrad(MinuitFit, upar, 0);

	ROOT::Minuit2::FunctionMinimum min2 = migrad();
	std::cout << min2.UserParameters() << std::endl;
	
	return 0;
}

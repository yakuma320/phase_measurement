#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <complex>
#include <ctime>
#include <iomanip>
#include <stdio.h> 
#include <TMath.h>

using namespace std;

const int Nvac=500;
double Evac[Nvac];
double CSvac[Nvac];

const double me=0.511E-3;
const double alpha=1.0/137.0;
const double alpha2=pow(1.0/137.0,2);
const double pi=TMath::Pi();
const double pi2=pow(TMath::Pi(),2);
const double Emin=1.65;
const	double epsilon = 1.0e-8;
const int Npoints=100;
const double hc=0.389379660E9;

double vacc(double W)
{
	  int bin=int((W-2.5)/0.0034);
		double vp=1.0;
		if(bin>=0&&bin<500)  vp=CSvac[bin]+(W-Evac[bin])*(CSvac[bin+1]-CSvac[bin])/(Evac[bin+1]-Evac[bin]);
		return vp;
}
double phi(double scos,double y)
{
    double theta=acos(scos);
    double sinp=sin(theta*(1-y));
    double ssin=sin(theta);
    double sinpy=sin(pi*y);
    return (pi*y*sinp)/(ssin*sinpy);
}
double rescs(double *ecms,double *par)
{
    double s=ecms[0]*ecms[0];
    double MJ=par[0];
    double GT=par[1];
    double GE=par[2];
    double res=par[3];
    double ecm=par[4];
    double phase=par[5]*TMath::Pi()/180.0;
    double ATE=par[6];
    double FF=par[7];
    double GF=par[2];
    double MJ2=par[0]*par[0];
    double GT2=par[1]*par[1];
    
    double beta = 2.0*alpha/pi*(2.0*log(ecms[0]/me)-1.0);
    double delta= 0.75*beta+alpha/pi*(pi2/3.0-0.5)+beta*beta*(9./32.-pi2/12.);
    double a2=pow((1.0-MJ2/s),2)+MJ2*GT2/s/s;
    double a=sqrt(a2);
    double scos=1./a*(MJ2/s-1.);
    double ssin=sqrt(1-scos*scos);
    double Xf=1.-4.*Emin*Emin/s;

    double C1=(8.*pi*alpha*sqrt(GE*GF)/MJ*( (s-MJ2)*(cos(phase)*ATE+1)+MJ*GT*sin(phase)*ATE)*sqrt(vacc(ecms[0])) +12*pi*(GE*GF/MJ2)*s*(1+ATE*ATE+2*ATE*cos(phase)))/s/s;
    double C2=(8.*pi*alpha*sqrt(GE*GF)/MJ*(cos(phase)*ATE+1)*sqrt(vacc(ecms[0]))+12*pi*(GE*GF/MJ2)*(1+ATE*ATE+2*ATE*cos(phase)))/s;
    
    double Xf0 = pow(Xf,beta)/(beta);
    double Xf1 = pow(Xf,beta-1.)/(beta-1.);
    double Xf2 = pow(Xf,beta-2.)/(beta-2.);
    double Xf3 = pow(Xf,beta-3.)/(beta-3.);
    double Xf4 = pow(Xf,beta-4.)/(beta-4.);
    double R2=2.*(s-MJ2)/s;
    double R3=a2*(4.*scos*scos-1.);

    double temp1=C1*(1+delta)*( pow(a,beta-2)*phi(scos,beta)+beta*(Xf2+Xf3*R2+Xf4*R3) );
    double temp2=( -1.0*beta*(1+delta)*C2-(beta+beta*beta/4.)*C1 )*( pow(a,beta-1)/(1+beta)*phi(scos,beta+1)+Xf1+Xf2*R2+Xf3*R3 );
    double temp3=( 0.5*log((Xf*Xf+2.*a*Xf*scos+a2)/a2)-scos/ssin*(atan((Xf+a*scos)/(a*ssin))-pi/2.+(asin(ssin)+acos(scos))/2.0) ) * ( (beta+beta*beta/4.)*C2+(beta/2.-beta*beta*3./8)*C1 ); 
    double inter=temp1+temp2+temp3;

    double dXS=(1./sqrt(2*pi)/res*exp(-pow((ecm-ecms[0]),2)/(2.*res*res)))*(inter)*hc;
    return dXS;

}
double qedxs(double s, double FF)
{
    double COEF=4.0*pi*alpha*alpha*hc/3.0/s*pow(FF/pow(s,1),2);
    return COEF;
}

double funs(double *x, double *par)
{
    double s=par[0]*par[0];
    double beta = 2.0*alpha/pi*(log(s/me/me)-1.0);
    double dXS = beta*TMath::Power(x[0],beta-1)*qedxs(s*(1-x[0]), par[1]);
    return dXS;
}

double funh(double *x, double *par)
{
    double s=par[0]*par[0];
    double beta = 2.0*alpha/pi*(log(s/me/me)-1.0);
    double Xf=1.-4.*Emin*Emin/s;

    double hard1=1.0-x[0]/2.0;
    double hard2=4.0*(2.0-x[0])*log(1.0/x[0])+1.0/x[0]*(1.0+3.0*(1.0-x[0])*(1.0-x[0]))*log(1.0/(1.0-x[0]))-6.0+x[0];
    double dXS=(-beta*hard1+beta*beta/8.0*hard2)*qedxs(s*(1-x[0]), par[1]);
    return dXS;
}

double intXS(double ecm, double par0, double par1, double par2, double par3, double par4, double par5, double par6)
{
    double res=par3;

		// the interference and resonance parts
    TF1 *fXS = new TF1("fXS", rescs, 3.580, 3.710, 7);
    fXS->SetParameter(0, par0); //mass of psip
    fXS->SetParameter(1, par1); //total width of psip
    fXS->SetParameter(2, par2); //Gamma_ee
    fXS->SetParameter(3, par3); //energy spread
    fXS->SetParameter(4, ecm); //initial energy
    fXS->SetParameter(5, par4); //phase
    fXS->SetParameter(6, par5); //ATE
    fXS->SetParameter(7, par6); //FF
    double my_intXS=fXS->Integral(TMath::Min(ecm,par0)-10*res, TMath::Max(ecm,par0)+10*res)*pow(par6/pow(ecm,2),2);

		// the QED (continuume) part are seperated into soft and hard radiation parts
    double beta = 2.0*alpha/pi*(2.0*log(ecm/me)-1.0);
    double delta= 0.75*beta+alpha/pi*(pi2/3.0-0.5)+beta*beta*(9./32.-pi2/12.); 
    double UPL = 1.0-4.0*Emin*Emin/(ecm*ecm);
    TF1 *fXS_funs = new TF1("fXS_funs", funs, 0.0, 1.0, 2);
    fXS_funs->SetParameter(0, ecm); //initial energy
    fXS_funs->SetParameter(1, par6); //FF
    double my_intXS_qed1 = fXS_funs->Integral(epsilon, UPL);
    double my_intXS_qed2 = TMath::Power(epsilon,beta)*qedxs(ecm*ecm, par6);
    TF1 *fXS_funh = new TF1("fXS_funh", funh, 0.0, 1.0, 2);
    fXS_funh->SetParameter(0, ecm); //initial energy
    fXS_funh->SetParameter(1, par6); //FF
    double my_intXS_qed3 = fXS_funh->Integral(0.0, UPL);
    double my_intXS_qed = ( (my_intXS_qed1 + my_intXS_qed2)*(1.0+delta) + my_intXS_qed3 );
    return my_intXS+my_intXS_qed;
}

int main()
{
    ifstream vacfile("vacc_nopsip_gen.dat");//,ios::in);
    for(int i=0;i<Nvac;i++)
    {
        vacfile>>Evac[i]>>CSvac[i];
    }

		double MJ=3.686;
		double GT=294E-6;
		double GE=2.33e-6;
		double phi=0;
		double res=1.38282e-03;
		char name[40];
		for(int num1=0; num1<10; num1++){
				for(int num2=0; num2<10; num2++){
						double ATE=(11.0-1.0)/10*num1+1.0;
						double FF=(2.1-0.1)/10*num2+0.1;
						sprintf(name, "cal_phase0_%d_%d_log", num1, num2);
						ofstream f_out;
						f_out.open(name, ios::out);

						double x[Npoints]={0}, R[Npoints]={0};
						double xs_ana[Npoints]={0}, exs_ana[Npoints]={0};
						for(int i=0; i<Npoints; i++){
								x[i]=3.67+(3.70-3.67)/Npoints*i;
								xs_ana[i]=intXS(x[i], MJ, GT, GE, res, phi, ATE, FF);
								f_out<<x[i]<<"   "<<xs_ana[i]<<endl;
						}
						f_out.close();
				}
		}
		return 1;
}

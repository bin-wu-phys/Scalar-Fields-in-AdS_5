using namespace std;
#include "../Chebyshev.h"
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#define VAC
#define M2 0 
#define L 2.5
#define EPS 0.5
#define NSAVE 1.0
#define B 10 
#define DU 0.001

double phi0(double t){
	return EPS*exp(-B*t*t)/B;
	
}

double phi0dot(double t){
		return -2.0*t*EPS*exp(-B*t*t);

}

int main(){
	ofstream out;
	ostringstream ostr;

	int n=256;double Dm=2.0-sqrt(4.0+M2),coef=1.0;
	DoubleVector x=L*(1.0+ChebPoints(n))+DU;
	DoubleMatrix D=ChebDiff(n)/L,d=D.sub(0,n-1,0,n-1);
	DoubleVector xx=x.sub(0,n-1);
	d.LUdcmp();

	// Initial conditions
	DoubleVector Phi=0*x,II=0*x,phi=0*x;
	
	DoubleVector delta,A,ddelta;DoubleMatrix Ap;
	delta=0*x;
	
	A=( delta + 1.0 );
	
	double dt=1.0e-5/NSAVE, tmax=3.0,t=-3.0, tsave=.1;
	unsigned int nSave=0,numsave=(int)(tsave/dt)+1;dt=tsave/numsave;
	
#ifdef VAC
        ostr << "VAC";
#endif
	ostr << "N" << n;
	ostr << "L" << setfill('0') <<  setw(4) << int(100.0*L);
	if(M2>=0)
		ostr << "m" << setw(4)<< int(M2);
	else
		ostr << "mm" << setw(4)<< int(-M2);
	
	ostr << "e"  << setfill('0') << setw(4) << int(10000.0*EPS);
	ostr << "a"  << setfill('0') << setw(3) << int(B);

	out.open(string(ostr.str()+".dat").c_str(),ios::app);
	out << "#n = " << n << ", dt = " << dt << ", B = " << B << ", L = " << L << ", EPS = " << EPS << "\n";
	out << "#t = " << t << "\n" << ((DoubleMatrix)x&Phi&II&delta&A&phi);//a-aold);//
	out.close();

	double running_time=clock()/CLOCKS_PER_SEC;
		
	// Time marching
	DoubleVector dII,dPhi,dA,dAnm1,dIInm1,dPhinm1,dAnm2,dIInm2,dPhinm2;//Time-marching
	DoubleVector Ad,AdII,Adphi,DAdphi;//temp vectors to speed up;
	
	dAnm2=Phi;dPhinm2=Phi;dIInm2=Phi;
	dAnm1=Phi;dPhinm1=Phi;dIInm1=Phi;
	
	double d4o3 = (4.0/3.0),ddto12=(dt/12.0),d2o3=(2.0/3.0);
	DoubleVector emd;
	
	do{
		// Adams-Bashforth 3rd
		emd=exp(-delta);
		Ad=A*emd;AdII=Ad*II;
		dA=( d4o3*x*A*Phi*AdII );
		dPhi=( D*(AdII) );
		Adphi = ( Ad*Phi ); DAdphi=(D*Adphi);
		dII = (( (DAdphi.sub(0,n-1)) - 3.0*Adphi.sub(0,n-1)/xx )&0)-M2*emd*(pow(xx,-2.0)&0)*phi;

#ifndef VAC             						
		A= ( A+ (ddto12*(23.0*dA -16.0*dAnm1+5.0*dAnm2)) );
#endif

		Phi= ( Phi+ ( ddto12*(23.0*dPhi -16.0*dPhinm1+5.0*dPhinm2)) );
		II= ( II+ (ddto12*(23.0*dII -16.0*dIInm1+5.0*dIInm2)) );
		
		//II[n]=phi0dot(t);//v3
		//II[n]=phi0dot(t);II[0]=0;//v1
		II[n]=coef*phi0dot(t);Phi[0]=0;//v2

#ifndef VAC             		
		ddelta=(d2o3*xx*( ((Phi*Phi)+(II*II)).sub(0,n-1) ) );
		delta=((d|ddelta)&0);
#endif

		phi=((d|Phi.sub(0,n-1))&0)+coef*phi0(t);
		
		for(int i=0;i<A.size();i++){
			if(A[i]<0.01){
				out.open(string(ostr.str()+".dat").c_str(),ios::app);
				out << "\n#t = " << t << "\n" << ((DoubleMatrix)x&Phi&II&delta&A&phi);
				out.close();
				t=tmax;break;nSave=0;				
			}else if(isnan(A[i])){
				t=tmax;break;nSave=0;
			}

		}
		
		//if(t<=0){	
		if(nSave==numsave){
				out.open(string(ostr.str()+".dat").c_str(),ios::app);
				out << "\n#t = " << t << "\n" << ((DoubleMatrix)x&Phi&II&delta&A&phi);
				out.close();
				nSave=0;
		}
		/*}else{
			if(nSave>=100*NSAVE){
				out << "\n#t = " << t << "\n" << ((DoubleMatrix)x&Phi&II&delta&A);
				nSave=0;
			}
		}*/
		t+=dt;nSave++;dPhinm2=dPhinm1;dPhinm1=dPhi;dIInm2=dIInm1;dIInm1=dII;dAnm2=dAnm1;dAnm1=dA;
		
	}while(t<=tmax);

	running_time=clock()/CLOCKS_PER_SEC - running_time;
	out.open(string(ostr.str()+".dat").c_str(),ios::app);
	if(running_time>=60.0){
		if(running_time>=3600){
			out << "#Calculation is done and the running time is " << running_time/3600.0 << " hours.\n\n\n";
		}
		else {
			out << "#Calculation is done and the running time is " << running_time/60.0 << " minutes.\n\n\n";
		}
	}
	else {
		out << "#Calculation is done and the running time is " << running_time << " seconds.\n\n\n";
	}
	out.close();
		
	return 0;
}

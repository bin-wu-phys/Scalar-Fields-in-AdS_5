using namespace std;
#include "../DoubleVector.h"
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

//#define VAC
#define LAMBDA 0.5	
#define M2 5 
#define DU 0.001
#define NSAVE 1.0

double MM=sqrt(4+M2);//mm used in Mathematica file
double B=10/(LAMBDA*LAMBDA); 
double EPS=0.005*pow(DU,2-MM)*(pow(LAMBDA,MM-4));
double L=3.0*LAMBDA;

double phi0(double t){
	return EPS*exp(-B*t*t)/B;
	
}

double phi0dot(double t){
	return -2.0*t*EPS*exp(-B*t*t);
	
}

int main(){
	ofstream out;
	ostringstream ostr;
	
	int n=(int)(NSAVE*10000),i;double dx=L/n;
	DoubleVector x(n+1);for(i=0;i<=n;i++) x[i]=dx*i+DU;

	// Initial conditions
	DoubleVector Phi=0*x,II=0*x,phi=Phi;
	
	DoubleVector delta,A,ddelta;DoubleMatrix Ap;
	delta=0*x;
	
	A=( delta + 1.0 );
	
	double dt=0.1*dx, tmax=LAMBDA*3.0,
t=-LAMBDA*3.0,tsave=LAMBDA*1.0;
	unsigned int nSave=0,numsave=(int)(tsave/dt)+1;dt=tsave/numsave;
#ifdef VAC
	ostr << "VAC";
#endif	
	ostr << "Nd" << n;
	ostr << "L" << setfill('0') <<  setw(4) << int(100.0*L);
	ostr << "U" << setfill('0') <<  setw(4) << int(10000.0*DU);
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
	DoubleVector dII,dPhi(n+1),dA,dAnm1,dIInm1,dPhinm1,dAnm2,dIInm2,dPhinm2;
	dPhi[n]=0;
	DoubleVector Ad,AdII,Adphi,DAdphi(n+1);//temp vectors to speed up;
	
	dAnm2=Phi;dPhinm2=Phi;dIInm2=Phi;dAnm1=Phi;dPhinm1=Phi;dIInm1=Phi;

	double d4o3 = (4.0/3.0),ddto12=dt/12.0,d2o3=(2.0/3.0);

	double umin=L;DoubleVector emd;
	
	do{
		// Adams-Bashforth 3rd
		emd=exp(-delta);
		Ad=A*emd;AdII=Ad*II;
		dA=( d4o3*x*A*Phi*AdII );
		for(i=0;i<n;i++)dPhi[i]=( AdII[i+1]-AdII[i] )/dx;
		Adphi = ( Ad*Phi ); for(i=1;i<=n;i++) DAdphi[i]=(Adphi[i]-Adphi[i-1])/dx;
		dII = (0&( (DAdphi.sub(1,n)) - 3.0*Adphi.sub(1,n)/x.sub(1,n) ))-M2*emd*(0&pow(x.sub(1,n),-2.0))*phi;
#ifndef VAC		
		A= ( A+ (ddto12*(23.0*dA -16.0*dAnm1+5.0*dAnm2)) );
#endif
		Phi= ( Phi+ ( ddto12*(23.0*dPhi -16.0*dPhinm1+5.0*dPhinm2)) );
		II= ( II+ (ddto12*(23.0*dII -16.0*dIInm1+5.0*dIInm2)) );
		
		II[0]=phi0dot(t);Phi[n]=0;//v2
		
		ddelta=(d2o3*x.sub(1,n)*( ((Phi*Phi)+(II*II)).sub(1,n) ) );
		phi[0]=phi0(t);
		for(i=1;i<=n;i++){
			int im=i-1;
#ifndef VAC 
			delta[i]=dx*ddelta[im]+delta[im];
#endif
			phi[i]=dx*Phi[im]+phi[im];
		}
		
		t+=dt;nSave++;dPhinm2=dPhinm1;dPhinm1=dPhi;dIInm2=dIInm1;dIInm1=dII;dAnm2=dAnm1;dAnm1=dA;
		
		for(int i=0;i<A.size();i++){
			if(A[i]<0.01&&x[i]<umin-0.1){
				umin=x[i];
				out.open(string(ostr.str()+".dat").c_str(),ios::app);
				out << "\n#t = " << t << "\n" << ((DoubleMatrix)x&Phi&II&delta&A&phi);
				out.close();
				t=tmax;break;nSave=0;				
			}else if(isnan(A[i])){
				t=tmax;break;nSave=0;
			}

		}
		if(nSave==numsave){
			out.open(string(ostr.str()+".dat").c_str(),ios::app);
			out << "\n#t = " << t << "\n" << ((DoubleMatrix)x&Phi&II&delta&A&phi);
			out.close();
			nSave=0;
		}
	}while(t<tmax);
	
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

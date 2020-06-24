using namespace std;
#include "../Chebyshev.h"

#define SIG 1.0/16.0
#define EPS 10.0
//#define ZERO 1.0e-50
#define NSAVE 1.0

double phi0dot(double t){
		return -2.0*t*EPS*exp(-10*t*t);
}

int main(){
	int n=64;
	DoubleVector x=0.25*PI*(1.0+ChebPoints(n));
	DoubleMatrix D=ChebDiff(n)/(0.25*PI);
	
	// Initial conditions
	DoubleVector Phi=0*x,II=0*x;//((2.0*EPS*exp(-4.0*(cot(x.sub(1,n))^2)/(PI*PI*SIG*SIG))/PI)&0);// Note: () for ^
	
	double dt=1.0e-5/NSAVE, tmax=5000000000.0*NSAVE*dt, t=-10;
	unsigned int nSave=0;
	
	cout << "#n = " << n << ", dt = " << dt << "\n\n";
	
	cout << "#t = " << t << "\n" << ((DoubleMatrix)x&Phi&II);//a-aold);//

		
	// Time marching
	DoubleVector dII,dPhi,dIInm1,dPhinm1,dIInm2,dPhinm2;
	
	dPhinm2=Phi;dIInm2=Phi;dPhinm1=Phi;dIInm1=Phi;
	
	do{
		// Adams-Bashforth 3rd
		dPhi=( (cos(x)^2)*D*(II) );
		dII = ( (cos(x)^2)*D*(Phi) );
						
		Phi= ( Phi+ ( (dt/12.0)*(23.0*dPhi -16.0*dPhinm1+5.0*dPhinm2)) );
		II= ( II+ ((dt/12.0)*(23.0*dII -16.0*dIInm1+5.0*dIInm2)) );
		
		//II[n]=phi0dot(t);Phi[n]=-II[n];
		II[n]=phi0dot(t);//v2
				
		t+=dt;nSave++;dPhinm2=dPhinm1;dPhinm1=dPhi;dIInm2=dIInm1;dIInm1=dII;
		
		if(nSave==1000*NSAVE){
			cout << "\n#t = " << t << "\n" << ((DoubleMatrix)x&Phi&II);
			nSave=0;
		}
		//getchar();
	}while(t<tmax);
	
	return 0;
}
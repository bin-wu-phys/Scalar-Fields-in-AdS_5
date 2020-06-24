=using namespace std;
#include "Chebyshev.h"

#define EPS 0.01
//#define ZERO 1.0e-50
#define NSAVE 1.0

double phi0(double t){
	return EPS*exp(-100*(t-1.0)*(t-1.0));
}


int main(){
	int n=128;
	DoubleVector x=0.25*PI*(1.0+ChebPoints(n));
	DoubleMatrix D=4.0*ChebDiff(n)/PI;

	// Initial conditions
	DoubleVector Phi=0*x;
	
	double dt=1.0e-3/NSAVE, tmax=100000000.0*NSAVE*dt, t=0;
	unsigned int nSave=0;
	
	cout << "#n = " << n << ", dt = " << dt << "\n\n";
	
	cout << "#t = " << t << "\n" << ((DoubleMatrix)x&Phi);//a-aold);//

		
	// Time marching
	DoubleVector dPhi,dPhinm1,dPhinm2;
	
	dPhinm2=Phi;dPhinm1=Phi;
	
	do{
		
		// Adams-Bashforth 3rd
		dPhi=-((cos(x)^2)*(D*Phi));
						
		Phi= ( Phi+ ( (dt/12.0)*(23.0*dPhi -16.0*dPhinm1+5.0*dPhinm2)) );
		
		Phi[n]=phi0(t);
		
		
		t+=dt;nSave++;dPhinm2=dPhinm1;dPhinm1=dPhi;
				
		if(nSave==100*NSAVE){
			cout << "\n#t = " << t << "\n" << ((DoubleMatrix)x&Phi);
			nSave=0;
		}
		//getchar();
	}while(t<tmax);

	return 0;
}

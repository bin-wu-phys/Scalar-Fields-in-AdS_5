using namespace std;
#include "../Chebyshev.h"

#define EPS 0.01
//#define ZERO 1.0e-50
#define NSAVE 1.0


int main(){
	int n=64;
	DoubleVector x=0.5*PI*ChebPoints(n);
	DoubleMatrix D=2.0*ChebDiff(n)/PI;

	// Initial conditions
	DoubleVector Phi=exp(-10.0*(tan(x)^2));
	
	double dt=1.0e-5/NSAVE, tmax=100000000.0*NSAVE*dt, t=0;
	unsigned int nSave=0;
	
	cout << "#n = " << n << ", dt = " << dt << "\n\n";
	
	cout << "#t = " << t << "\n" << ((DoubleMatrix)x&Phi);//a-aold);//

		
	// Time marching
	DoubleVector dPhi,dPhinm1,dPhinm2;
	
	dPhinm2=-((cos(x)^2)*(D*Phi));
	Phi=Phi+dt*dPhinm2;
	
	t+=dt;nSave++;
		
	dPhinm1=-((cos(x)^2)*(D*Phi));
	Phi=Phi+dt*dPhinm1;
	
	t+=dt;nSave++;
		
	do{
		
		// Adams-Bashforth 3rd
		dPhi=-((cos(x)^2)*(D*Phi));
						
		Phi= ( Phi+ ( (dt/12.0)*(23.0*dPhi -16.0*dPhinm1+5.0*dPhinm2)) );
				
		
		t+=dt;nSave++;dPhinm2=dPhinm1;dPhinm1=dPhi;
				
		if(nSave==1000*NSAVE){
			cout << "\n#t = " << t << "\n" << ((DoubleMatrix)x&Phi);
			nSave=0;
		}
		//getchar();
	}while(t<tmax);

	return 0;
}

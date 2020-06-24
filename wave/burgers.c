using namespace std;
#include "../Chebyshev.h"

#define EPS 0.01
//#define ZERO 1.0e-50
#define NSAVE 1.0


int main(){
	int n=128;
	DoubleVector x=ChebPoints(n);
	DoubleMatrix D=ChebDiff(n);

	// Initial conditions
	DoubleVector Phi=exp(-100.0*(x^2));
	
	double dt=1.0e-5/NSAVE, tmax=100000000.0*NSAVE*dt, t=0;
	unsigned int nSave=0;
	
	cout << "#n = " << n << ", dt = " << dt << "\n\n";
	
	cout << "#t = " << t << "\n" << ((DoubleMatrix)x&Phi);//a-aold);//

		
	// Time marching
	DoubleVector dPhi,dPhinm1,dPhinm2;
	
	dPhinm2=-0.5*(D*(Phi*Phi));
	Phi=Phi+dt*dPhinm2;
	
	t+=dt;nSave++;
		
	dPhinm1=-0.5*(D*(Phi*Phi));
	Phi=Phi+dt*dPhinm1;
	
	t+=dt;nSave++;
		
	do{
		
		// Adams-Bashforth 3rd
		dPhi=-0.5*(D*(Phi*Phi));
						
		Phi= ( Phi+ ( (dt/12.0)*(23.0*dPhi -16.0*dPhinm1+5.0*dPhinm2)) );
		
		Phi[n]=0;
				
		
		t+=dt;nSave++;dPhinm2=dPhinm1;dPhinm1=dPhi;
				
		if(nSave==1000*NSAVE){
			cout << "\n#t = " << t << "\n" << ((DoubleMatrix)x&Phi);
			nSave=0;
		}
		//getchar();
	}while(t<tmax);

	return 0;
}

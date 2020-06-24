using namespace std;
#include "../Chebyshev.h"

#define M2 100
#define L 1.5 
#define EPS 0.5
#define NSAVE 1.0
#define B 10.0

double phi0(double t){
	return EPS*exp(-B*t*t)/B;
	
}

double phi0dot(double t){
		return -2.0*t*EPS*exp(-B*t*t);

}

int main(){
	int n=128;
	DoubleVector x=L*(1.0+ChebPoints(n));
	DoubleMatrix D=ChebDiff(n)/L,d=D.sub(0,n-1,0,n-1);
	DoubleVector xx=x.sub(0,n-1);
	d.LUdcmp();

	// Initial conditions
	DoubleVector Phi=0*x,II=0*x,phi=0*x;
	
	double dt=1.0e-5/NSAVE, tmax=100.0,t=-3.0, tsave=1.0;
	unsigned int nSave=0,numsave=(int)(tsave/dt)+1;dt=tsave/numsave;
	
	cout << "#n = " << n << ", dt = " << dt << ", B = " << B << ", L = " << L << ", EPS = " << EPS << "\n";
	cout << "#t = " << t << "\n" << ((DoubleMatrix)x&Phi&II&phi);//a-aold);//

		
	// Time marching
	DoubleVector dII,dPhi,dIInm1,dPhinm1,dIInm2,dPhinm2;//Time-marching
	
	dPhinm2=Phi;dIInm2=Phi;
	dPhinm1=Phi;dIInm1=Phi;
	
	double d4o3 = (4.0/3.0),ddto12=(dt/12.0),d2o3=(2.0/3.0);
	
	do{
		// Adams-Bashforth 3rd
		dPhi=D*II;
		dII = (( (D*Phi).sub(0,n-1) - 3.0*Phi.sub(0,n-1)/xx )&0)-M2*pow(x,3.0)*phi;
						
		Phi= ( Phi+ ( ddto12*(23.0*dPhi -16.0*dPhinm1+5.0*dPhinm2)) );
		II= ( II+ (ddto12*(23.0*dII -16.0*dIInm1+5.0*dIInm2)) );
		
		//II[n]=phi0dot(t);//v3
		//II[n]=phi0dot(t);II[0]=0;//v1
		II[n]=phi0dot(t);Phi[0]=0;//v2
		
		phi=((d|Phi.sub(0,n-1))&0)+phi0(t);
		
		//if(t<=0){	
		if(nSave==numsave){
				cout << "\n#t = " << t << "\n" << ((DoubleMatrix)x&Phi&II&phi);
				nSave=0;
		}
		/*}else{
			if(nSave>=100*NSAVE){
				cout << "\n#t = " << t << "\n" << ((DoubleMatrix)x&Phi&II&delta&A);
				nSave=0;
			}
		}*/
		t+=dt;nSave++;dPhinm2=dPhinm1;dPhinm1=dPhi;dIInm2=dIInm1;dIInm1=dII;
		
	}while(t<=tmax);
	
	return 0;
}

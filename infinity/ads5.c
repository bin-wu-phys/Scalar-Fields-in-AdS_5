using namespace std;
#include "Chebyshev.h"

#define SIG 1.0/16.0
#define EPS 10.0
//#define ZERO 1.0e-50
#define NSAVE 1.0

double phi0(double t){
	return EPS*exp(-10000*t*t);
}

double phi0dot(double t){
		return -2.0*t*EPS*exp(-1.0*t*t);

}

int main(){
	int n=128;
	DoubleVector x=0.25*PI*(1.0+ChebPoints(n));
	DoubleVector xx=x.sub(0,n-1),xxx=x.sub(1,n-1);
	DoubleMatrix D=4.0*ChebDiff(n)/PI,d=D.sub(0,n-1,0,n-1),dd=d;
	for(int i=0;i<dd.col();i++) dd[0][i]*=5.0;
	dd=dd-(0&(8.0/sin(2.0*xxx)));
	//DoubleMatrix ddd=( (D.sub(1,n-1,1,n-1))+( 1.0+ (2.0*(sin(xxx)^2)))/( sin(xxx)*cos(xxx) ));
	d.LUdcmp();

	// Initial conditions
	DoubleVector Phi=0*x,II=0*x;//((2.0*EPS*exp(-4.0*(cot(x.sub(1,n))^2)/(PI*PI*SIG*SIG))/PI)&0);// Note: () for ^
	// Solve delta and A
	
	//cout << "# x n-by-n integration (n-1)-by-(n-1) err12 err13 err23\n";
	
	DoubleVector delta,A,a,ddelta;DoubleMatrix Ap;
	delta=0*x;
	//This is a better way to solve A
	
	A=( 0*x + 1.0 );//cout << ((DoubleMatrix)x&A);
	
	double dt=1.0e-5/NSAVE, tmax=5000000000.0*NSAVE*dt, t=-10;
	unsigned int nSave=0;
	
	cout << "#n = " << n << ", dt = " << dt << "\n\n";
	
	cout << "#t = " << t << "\n" << ((DoubleMatrix)x&Phi&II&delta&A);//a-aold);//

		
	// Time marching
	DoubleVector dII,dPhi,dA,dAnm1,dIInm1,dPhinm1,dAnm2,dIInm2,dPhinm2,dAnp1,dIInp1,dPhinp1, Adphi, DAdphi;
	
	dAnm2=Phi;dPhinm2=Phi;dIInm2=Phi;
	dAnm1=Phi;dPhinm1=Phi;dIInm1=Phi;
	
	do{
		// Adams-Bashforth 3rd
		dA=( 4.0*sin(x)*(cos(x)^3)*A*A*exp(-delta)*Phi*II );
		dPhi=( D*(A*(cos(x)^2)*exp(-delta)*II) );
		Adphi = ( A*exp(-delta)*(cos(x)^2)*Phi ); DAdphi=(D*Adphi);
		dII = (4.0*DAdphi[0])&( (DAdphi.sub(1,n-1)) - 6.0*Adphi.sub(1,n-1)/sin(2.0*xxx) )&(-2.0*DAdphi[n]);
						
		A= ( A+ ((dt/12.0)*(23.0*dA -16.0*dAnm1+5.0*dAnm2)) );
		Phi= ( Phi+ ( (dt/12.0)*(23.0*dPhi -16.0*dPhinm1+5.0*dPhinm2)) );
		II= ( II+ ((dt/12.0)*(23.0*dII -16.0*dIInm1+5.0*dIInm2)) );
		
		II[n]=phi0dot(t);//v3
		//II[n]=phi0dot(t);II[0]=0;//v1
		//II[n]=phi0dot(t);Phi[0]=0;//v2
		
		ddelta=(sin(2.0*xx)*( ((Phi*Phi)+(II*II)).sub(0,n-1) )/3.0 );
		delta=((d|ddelta)&0);
		
		t+=dt;nSave++;dPhinm2=dPhinm1;dPhinm1=dPhi;dIInm2=dIInm1;dIInm1=dII;dAnm2=dAnm1;dAnm1=dA;
		
		for(int i=0;i<A.size();i++){
			if(A[i]<0.001){
				cout << "\n#t = " << t << "\n" << ((DoubleMatrix)x&Phi&II&delta&A);
				t=tmax;break;nSave=0;				
			}else if(isnan(A[i])){
				t=tmax;break;nSave=0;
			}

		}
		
		if(nSave==100*NSAVE){
			//cout << ((DoubleMatrix)x&A&((A-AA)/A));
			cout << "\n#t = " << t << "\n" << ((DoubleMatrix)x&Phi&II&delta&A);
			nSave=0;
		}
		//getchar();
	}while(t<tmax);
	
	return 0;
}

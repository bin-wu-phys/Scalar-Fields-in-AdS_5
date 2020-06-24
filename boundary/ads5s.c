using namespace std;
#include "../Chebyshev.h"

#define EPS 0.01
#define NSAVE 1.0

double phi0dot(double t){
	return -2.0*t*EPS*exp(-10.0*t*t);
	
}

int main(){
	int n=128;
	DoubleVector x=0.25*PI*(1.0+ChebPoints(n));
	DoubleMatrix D=4.0*ChebDiff(n)/PI,d=D.sub(0,n-1,0,n-1);
	DoubleVector xx=x.sub(0,n-1),xxx=x.sub(1,n-1);
	d.LUdcmp();
	
	// Initial conditions
	DoubleVector Phi=0*x,II=0*x;
	
	DoubleVector delta,A,ddelta;DoubleMatrix Ap;
	delta=0*x;
	
	A=( delta + 1.0 );
	
	double dt=1.0e-5/NSAVE, tmax=10000000.0, t=-3;
	unsigned int nSave=0;
	
	cout << "#n = " << n << ", dt = " << dt << "\n";
	cout << "#t = " << t << "\n" << ((DoubleMatrix)x&Phi&II&delta&A);//a-aold);//
	
	
	// Time marching
	DoubleVector dII,dPhi,dA,dAnm1,dIInm1,dPhinm1,dAnm2,dIInm2,dPhinm2,dAnp1,dIInp1,dPhinp1, Adphi, DAdphi;
	
	dAnm2=Phi;dPhinm2=Phi;dIInm2=Phi;
	dAnm1=Phi;dPhinm1=Phi;dIInm1=Phi;
	
	do{
		// Adams-Bashforth 3rd
		dA=( (2.0/3.0)*sin(2.0*x)*A*A*exp(-delta)*Phi*II );
		dPhi=( D*(A*exp(-delta)*II) );
		Adphi = ( A*exp(-delta)*Phi ); DAdphi=(D*Adphi);
		dII = (4.0*DAdphi[0])&( (DAdphi.sub(1,n-1)) - 6.0*Adphi.sub(1,n-1)/sin(2.0*xxx) )&(-2.0*DAdphi[n]);
		
		A= ( A+ ((dt/12.0)*(23.0*dA -16.0*dAnm1+5.0*dAnm2)) );
		Phi= ( Phi+ ( (dt/12.0)*(23.0*dPhi -16.0*dPhinm1+5.0*dPhinm2)) );
		II= ( II+ ((dt/12.0)*(23.0*dII -16.0*dIInm1+5.0*dIInm2)) );
	
		II[n]=phi0dot(t);Phi[0]=0;
		
		ddelta=(sin(2.0*xx)*( ((Phi*Phi)+(II*II)).sub(0,n-1) )/3.0 );
		delta=((d|ddelta)&0);
		
		t+=dt;nSave++;dPhinm2=dPhinm1;dPhinm1=dPhi;dIInm2=dIInm1;dIInm1=dII;dAnm2=dAnm1;dAnm1=dA;
		
		for(int i=0;i<A.size();i++){
			if(A[i]<0.01){
				cout << "\n#t = " << t << "\n" << ((DoubleMatrix)x&Phi&II&delta&A);
				t=tmax;break;nSave=0;				
			}else if(isnan(A[i])){
				t=tmax;break;nSave=0;
			}

		}
		
		if(nSave==1000*NSAVE){
			double a4=0;DoubleVector A3=(D*(D*(D*A)));for(int i=0;i<=n;i++) a4+=D[n][i]*A3[i];a4/=24.0;
			cout << "\n#t = " << t << "\n" << ((DoubleMatrix)x&Phi&II&delta&A&(0*x+a4));
			nSave=0;
		}
	}while(t<tmax);
	
	return 0;
}

using namespace std;
#include "../Chebyshev.h"

#define SIG 1.0/16.0
#define EPS 40.0
#define NSAVE 1.0

int main(){
	int n=256;
	DoubleVector x=0.25*PI*(1.0+ChebPoints(n));
	DoubleVector xx=x.sub(1,n),xxx=x.sub(1,n-1);
	DoubleMatrix D=4.0*ChebDiff(n)/PI,d=D.sub(1,n,1,n),dd=d;
	for(int i=0;i<dd.col();i++) dd[n-1][i]*=-2.0;
	dd=dd-((2.0*(2.0+cos(2*xxx))/sin(2.0*xxx))&0);
	d.LUdcmp();

	// Initial conditions
	DoubleVector II= ((2.0*EPS*exp(-4.0*(cot(x.sub(0,n-1))^2)/(PI*PI*SIG*SIG))/PI)&0),Phi=0*x;//(0&-2.0*EPS*exp( -((cot(xxx)-R0)^2)/(SIG*SIG) )*(R0-cot(xxx))/(SIG*SIG*sqrt(3.0)*(sin(xxx)^2))&0)//(0&(2.0*EPS*exp(-4.0/(PI*PI*SIG*SIG*(sin(2.0*x.sub(1,n-1))^2)))/PI)&0);//

	// Solve delta and A	
	DoubleVector delta,A,a,ddelta;DoubleMatrix Ap;
	ddelta=(sin(2.0*xx)*( ((Phi*Phi)+(II*II)).sub(1,n) )/2.0 );
	delta=(0&(d|ddelta));
	Ap=(dd-ddelta);
	Ap.LUdcmp();
	a=(0&(Ap|(ddelta)));
	A=( a + 1.0 );//cout << ((DoubleMatrix)x&A);
	
	double dt=1.0e-5/NSAVE, tmax=5000000000.0*NSAVE*dt, t=0;
	unsigned int nSave=0;
	
	cout << "#n = " << n << ", dt = " << dt << "\n";
	
	cout << "#t = " << t << "\n" << ((DoubleMatrix)x&Phi&II&delta&A);//a-aold);//

		
	// Time marching
	DoubleVector dII,dPhi,dA,dAnm1,dIInm1,dPhinm1,dAnm2,dIInm2,dPhinm2;
	
	dAnm2=( sin(2.0*x)*A*A*exp(-delta)*Phi*II );
	dPhinm2=( D*(A*exp(-delta)*II) );
	DoubleVector Adphi = ( A*exp(-delta)*Phi ), DAdphi=(D*Adphi);
	dIInm2 = (3.0*DAdphi[0])&( (DAdphi.sub(1,n-1)) - 4.0*Adphi.sub(1,n-1)/sin(2.0*xxx) )&(-1.0*DAdphi[n]);
	
	dIInm2[n]=0;dPhinm2[0]=0;
	
	A= ( A+( dt*dAnm2 ) );
	Phi= ( Phi+( dt*dPhinm2 ) );
	II= ( II+( dt*dIInm2 ) );
	
	ddelta=(sin(2.0*xx)*( ((Phi*Phi)+(II*II)).sub(1,n) )/2.0 );
	delta=(0&(d|ddelta));
	
	t+=dt;nSave++;//cout << ((DoubleMatrix)x&Phi&II&delta&A&(ddelta&0));//cout << "The first step forward is done\n";

	dAnm1=( sin(2.0*x)*A*A*exp(-delta)*Phi*II );
	dPhinm1=( D*(A*exp(-delta)*II) );
	Adphi = ( A*exp(-delta)*Phi ); DAdphi=(D*Adphi);
	dIInm1 = (3.0*DAdphi[0])&( (DAdphi.sub(1,n-1)) - 4.0*Adphi.sub(1,n-1)/sin(2.0*xxx) )&(-1.0*DAdphi[n]);
	
	dIInm1[n]=0;dPhinm1[0]=0;	

	A= ( A+( dt*dAnm1 ) );
	Phi= ( Phi+( dt*dPhinm1 ) );
	II= ( II+( dt*dIInm1 ) );
	
	ddelta=(sin(2.0*xx)*( ((Phi*Phi)+(II*II)).sub(1,n) )/2.0 );
	delta=(0&(d|ddelta));

	t+=dt;nSave++;//cout << ((DoubleMatrix)x&Phi&II&delta&A&(ddelta&0));//cout << "The second step forward is done\n";
	
	do{
		// Adams-Bashforth 3rd
		dA=( sin(2.0*x)*A*A*exp(-delta)*Phi*II );
		dPhi=( D*(A*exp(-delta)*II) );
		Adphi = ( A*exp(-delta)*Phi ); DAdphi=(D*Adphi);
		dII = (3.0*DAdphi[0])&( (DAdphi.sub(1,n-1)) - 4.0*Adphi.sub(1,n-1)/sin(2.0*xxx) )&(-1.0*DAdphi[n]);
		
		dII[n]=0;dPhi[0]=0;
				
		A= ( A+ ((dt/12.0)*(23.0*dA -16.0*dAnm1+5.0*dAnm2)) );
		Phi= ( Phi+ ( (dt/12.0)*(23.0*dPhi -16.0*dPhinm1+5.0*dPhinm2)) );
		II= ( II+ ((dt/12.0)*(23.0*dII -16.0*dIInm1+5.0*dIInm2)) );
	
		ddelta=(sin(2.0*xx)*( ((Phi*Phi)+(II*II)).sub(1,n) )/2.0 );
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
			cout << "\n#t = " << t << "\n" << ((DoubleMatrix)x&Phi&II&delta&A);
			nSave=0;
		}
	}while(t<tmax);
	
	return 0;
}

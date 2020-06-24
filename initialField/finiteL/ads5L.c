using namespace std;
#include "Chebyshev.h"

#define L 10
#define EPS 1.0
//#define ZERO 1.0e-50
#define NSAVE 1.0

int main(){
	int n=6;
	DoubleVector x=L*(1.0+ChebPoints(n)),xx=x.sub(0,n-1);
	DoubleMatrix D=ChebDiff(n)/L,d=D.sub(0,n-1,0,n-1),dd=d;
	dd=dd-(4.0/xx);
	d.LUdcmp();

	// Initial conditions
	DoubleVector Phi=0*x,II=EPS*exp(-10*((x-4)^2));// Note: () for ^
	// Solve delta and A
	
	DoubleVector delta,A,a,ddelta;DoubleMatrix Ap;

	ddelta=((2.0/3.0)*xx*( ((Phi*Phi)+(II*II)).sub(0,n-1) ) );
	delta=((d|ddelta)&0);

	Ap=(dd-ddelta);
	Ap.LUdcmp();
	a=((Ap|(ddelta))&0);
	A=( a + 1.0 );
	cout << ((DoubleMatrix)x&A);
	cout << d<<dd<<Ap;
	cout << Ap.det();
	/*(
	double dt=1.0e-5/NSAVE, tmax=10.0, t=0;
	unsigned int nSave=0;
	
	cout << "#n = " << n << ", dt = " << dt << "\n";
	
	cout << "#t = " << t << "\n" << ((DoubleMatrix)x&Phi&II&delta&A);//a-aold);//
		
	// Time marching
	DoubleVector dII,dPhi,dA,dAnm1,dIInm1,dPhinm1,dAnm2,dIInm2,dPhinm2;
	
	dAnm2=( (4.0/3.0)*x*A*A*exp(-delta)*Phi*II );
	dPhinm2=( D*(A*exp(-delta)*II) );
	DoubleVector Adphi = ( A*exp(-delta)*Phi ), DAdphi=(D*Adphi);
	dIInm2 = (( (DAdphi.sub(0,n-1)) - 3.0*Adphi.sub(0,n-1)/xx )&(-2.0*DAdphi[n]));
		
	dIInm2[n]=0;dPhinm2[0]=0;
	
	A= ( A+( dt*dAnm2 ) );
	Phi= ( Phi+( dt*dPhinm2 ) );
	II= ( II+( dt*dIInm2 ) );
	
	ddelta=((2.0/3.0)*xx*( ((Phi*Phi)+(II*II)).sub(0,n-1) ) );
	delta=((d|ddelta)&0);
	
	t+=dt;nSave++;//cout << ((DoubleMatrix)x&Phi&II&delta&A&(ddelta&0));//cout << "The first step forward is done\n";

	dAnm1=( (4.0/3.0)*x*A*A*exp(-delta)*Phi*II );
	dPhinm1=( D*(A*exp(-delta)*II) );
	Adphi = ( A*exp(-delta)*Phi ); DAdphi=(D*Adphi);
	dIInm1 = (( (DAdphi.sub(0,n-1)) - 3.0*Adphi.sub(0,n-1)/xx )&(-2.0*DAdphi[n]));
	
	dIInm1[n]=0;dPhinm1[0]=0;	

	A= ( A+( dt*dAnm1 ) );
	Phi= ( Phi+( dt*dPhinm1 ) );
	II= ( II+( dt*dIInm1 ) );
	
	ddelta=((2.0/3.0)*xx*( ((Phi*Phi)+(II*II)).sub(0,n-1) ) );
	delta=((d|ddelta)&0);

	t+=dt;nSave++;//cout << ((DoubleMatrix)x&Phi&II&delta&A&(ddelta&0));//cout << "The second step forward is done\n";
	
	do{
		// Adams-Bashforth 3rd
		dA=( (4.0/3.0)*x*A*A*exp(-delta)*Phi*II );
		dPhi=( D*(A*exp(-delta)*II) );
		Adphi = ( A*exp(-delta)*Phi ); DAdphi=(D*Adphi);
		dII = (( (DAdphi.sub(0,n-1)) - 3.0*Adphi.sub(0,n-1)/xx )&(-2.0*DAdphi[n]));
						
		A= ( A+ ((dt/12.0)*(23.0*dA -16.0*dAnm1+5.0*dAnm2)) );
		Phi= ( Phi+ ( (dt/12.0)*(23.0*dPhi -16.0*dPhinm1+5.0*dPhinm2)) );
		II= ( II+ ((dt/12.0)*(23.0*dII -16.0*dIInm1+5.0*dIInm2)) );

		II[n]=0;Phi[0]=0;//v2
	
		ddelta=((2.0/3.0)*xx*( ((Phi*Phi)+(II*II)).sub(0,n-1) ) );
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
		
		if(nSave==100*NSAVE){
			cout << "\n#t = " << t << "\n" << ((DoubleMatrix)x&Phi&II&delta&A);
			nSave=0;
		}
	}while(t<tmax);
	*/
	return 0;
}

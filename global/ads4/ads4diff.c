using namespace std;
#include "../DoubleVector.h"

#define SIG 1.0/16.0
#define EPS 45.0
//#define ZERO 1.0e-50
#define NSAVE 1.0

int main(){
	int n=(int)(NSAVE*10000),i;double dx=L/n;
	DoubleVector x(n+1);for(i=0;i<=n;i++) x[i]=dx*i+0.005;

	// Initial conditions
	DoubleVector Phi=0*x,II=(0&(2.0*EPS*exp(-4.0*(tan(x.sub(1,n))^2)/(PI*PI*SIG*SIG))/PI));// Note: () for ^
	// Solve delta and A
	
	//cout << "# x n-by-n integration (n-1)-by-(n-1) err12 err13 err23\n";
	
	DoubleVector delta,A,a,ddelta,ddel;DoubleMatrix Ap;
	ddel=( -(sin(x)*cos(x)*( (Phi*Phi)+(II*II) )) );
	ddelta=( ddel.sub(0,n-1) );
	delta=((d|ddelta)&0);
	//This is a better way to solve A
	Ap=(dd-ddelta);
	Ap.LUdcmp();
	a=((Ap|(ddelta))&0);
	A=( a + 1.0 );//cout << ((DoubleMatrix)x&A);
	
	double dt=1.0e-6/NSAVE, tmax=5000000000.0*NSAVE*dt, t=0;
	unsigned int nSave=0;
	
	cout << "#n = " << n << ", dt = " << dt << "\n\n";
	
	cout << "#t = " << 0 << "\n" << ((DoubleMatrix)x&Phi&II&delta&A);//a-aold);//

		
	// Time marching
	DoubleVector dII,dPhi,dA,dAnm1,dIInm1,dPhinm1,dAnm2,dIInm2,dPhinm2,dAnp1,dIInp1,dPhinp1;

	dAnm2=( (-2.0)*sin(x)*cos(x)*A*A*exp(-delta)*Phi*II );
	dPhinm2=( D*(A*exp(-delta)*II) );
	DoubleVector Adphi = ( A*exp(-delta)*Phi ), DAdphi=(D*Adphi);
	dIInm2 = (-DAdphi[0])&( ( 2.0*(Adphi.sub(1,n-1))/(sin(xxx)*cos(xxx)) ) + (DAdphi.sub(1,n-1)) )&(3.0*DAdphi[n]);
	
	
	dIInm2[0]=0;dAnm2[0]=0;dAnm2[n]=0;dPhinm2[n]=0;//dPhinm2[0]=0;
	
	A= ( A+( dt*dAnm2 ) );
	Phi= ( Phi+( dt*dPhinm2 ) );
	II= ( II+( dt*dIInm2 ) );
	
	ddelta=( -(sin(xx)*cos(xx)*( ((Phi*Phi)+(II*II)).sub(0,n-1) ) ) );
	delta=((d|ddelta)&0);
	//Ap=(dd-ddelta);Ap.LUdcmp();a=((Ap|(ddelta))&0);A=( a + 1.0 );//cout << ((DoubleMatrix)x&A);
	
	t+=dt;nSave++;//cout << ((DoubleMatrix)x&Phi&II&delta&A&(ddelta&0));//scout << "The first step forward is done\n";

	dAnm1=( (-2.0)*sin(x)*cos(x)*A*A*exp(-delta)*Phi*II );
	dPhinm1=( D*(A*exp(-delta)*II) );
	Adphi = ( A*exp(-delta)*Phi ); DAdphi=(D*Adphi);
	dIInm1 = (-DAdphi[0])&( ( 2.0*(Adphi.sub(1,n-1))/(sin(xxx)*cos(xxx)) ) + (DAdphi.sub(1,n-1)) )&(3.0*DAdphi[n]);
	
	dIInm1[0]=0;dAnm1[0]=0;dAnm1[n]=0;dPhinm1[n]=0;//dPhinm1[0]=0;

	A= ( A+( dt*dAnm1 ) );
	Phi= ( Phi+( dt*dPhinm1 ) );
	II= ( II+( dt*dIInm1 ) );
	
	ddelta=( -(sin(xx)*cos(xx)*( ((Phi*Phi)+(II*II)).sub(0,n-1) ) ) );
	delta=((d|ddelta)&0);
	//Ap=(dd-ddelta);Ap.LUdcmp();a=((Ap|(ddelta))&0);A=( a + 1.0 );//cout << ((DoubleMatrix)x&A);

	t+=dt;nSave++;//cout << ((DoubleMatrix)x&Phi&II&delta&A&(ddelta&0));//cout << "The second step forward is done\n";
	
	do{
		
		// Adams-Bashforth 3rd
		dA=( (-2.0)*sin(x)*cos(x)*A*A*exp(-delta)*Phi*II );
		dPhi=( D*(A*exp(-delta)*II) );
		Adphi = ( A*exp(-delta)*Phi ); DAdphi=(D*Adphi);
		dII = (-DAdphi[0])&( ( 2.0*(Adphi.sub(1,n-1))/(sin(xxx)*cos(xxx)) ) + (DAdphi.sub(1,n-1)) )&(3.0*DAdphi[n]);
		
		dII[0]=0;dA[0]=0;dA[n]=0;dPhi[n]=0;//dPhi[0]=0;
				
		A= ( A+ ((dt/12.0)*(23.0*dA -16.0*dAnm1+5.0*dAnm2)) );
		Phi= ( Phi+ ( (dt/12.0)*(23.0*dPhi -16.0*dPhinm1+5.0*dPhinm2)) );
		II= ( II+ ((dt/12.0)*(23.0*dII -16.0*dIInm1+5.0*dIInm2)) );
		
		ddelta=( -(sin(xx)*cos(xx)*( ((Phi*Phi)+(II*II)).sub(0,n-1) ) ) );
		delta=((d|ddelta)&0);
				
		t+=dt;nSave++;dPhinm2=dPhinm1;dPhinm1=dPhi;dIInm2=dIInm1;dIInm1=dII;dAnm2=dAnm1;dAnm1=dA;
		
		for(int i=0;i<A.size();i++){
			if(A[i]<0.01){
				cout << "\n#t = " << t << "\n" << ((DoubleMatrix)x&Phi&II&delta&A);
				t=tmax;break;nSave=0;				
			}else if(isnan(A[i])){
				t=tmax;
			}

		}
		
		if(nSave==1000*NSAVE){
			//cout << ((DoubleMatrix)x&A&((A-AA)/A));
			cout << "\n#t = " << t << "\n" << ((DoubleMatrix)x&Phi&II&delta&A);
			nSave=0;
		}
		//getchar();
	}while(t<tmax);
	
	return 0;
}

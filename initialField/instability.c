using namespace std;
#include "Chebyshev.h"

#define SIG 1.0/16.0
#define EPS 45.0
//#define ZERO 1.0e-50
#define NSAVE 1000.0

int main(){
	int n=128;
	DoubleVector x=0.25*PI*(1.0+ChebPoints(n));
	DoubleVector xx=x.sub(0,n-1),xxx=x.sub(1,n-1);
	DoubleMatrix D=4.0*ChebDiff(n)/PI,d=D.sub(0,n-1,0,n-1),dd=d;
	for(int i=0;i<dd.col();i++) dd[0][i]*=-2.0;
	dd=dd+ (0&(( 1.0+ (2.0*(sin(xxx)^2)))/( sin(xxx)*cos(xxx) )));
	//DoubleMatrix ddd=( (D.sub(1,n-1,1,n-1))+( 1.0+ (2.0*(sin(xxx)^2)))/( sin(xxx)*cos(xxx) ));
	d.LUdcmp();

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
	
	//cout << Ap << A;
	
	//DoubleVector aold, intg = (0&(ddelta.sub(1,n-1)*sin(xxx)*exp(-delta.sub(1,n-1))/(cos(xxx)^3)) );
	//a= (((cos(xx)^3.0)*exp(delta.sub(0,n-1))*( (d|intg))/sin(xx))&0);
	//DoubleVector Amid=( a + 1.0 );//cout << ((DoubleMatrix)x&a&aold&( (((a-aold).sub(0,n-1))/aold.sub(0,n-1))&0 ));aold=a;
	
	//Ap=(ddd-(ddelta.sub(1,n-1)));Ap.LUdcmp();

	//DoubleVector Anew=( (0&(Ap|(ddelta.sub(1,n-1)))&0) + 1.0 );
	//cout << ((DoubleMatrix)x&A&Amid&Anew&((A-Amid)/A)&((A-Anew)/A)&((Anew-Amid)/Anew));
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
		
		/*
		// Runge-Kuta method
		DoubleVector dAk1,dPhik1,dIIk1,Ak1,Phik1,IIk1,dAk2,dPhik2,dIIk2,Ak2,Phik2,IIk2,dAk3,dPhik3,dIIk3,Ak3,Phik3,IIk3,dAk4,dPhik4,dIIk4,Ak4,Phik4,IIk4;
		
		dAk1=dt*( (-2.0)*sin(x)*cos(x)*A*A*exp(-delta)*Phi*II );
		dPhik1=dt*( D*(A*exp(-delta)*II) );
		Adphi = ( A*exp(-delta)*Phi ); DAdphi=(D*Adphi);
		dIIk1= dt*( (-DAdphi[0])&( ( 2.0*(Adphi.sub(1,n-1))/(sin(xxx)*cos(xxx)) ) + (DAdphi.sub(1,n-1)) )&(3.0*DAdphi[n]) );
		
		Ak1=A+0.5*dAk1;Phik1=Phi+0.5*dPhik1;IIk1=II+0.5*dIIk1;
		
		dAk2=dt*( (-2.0)*sin(x)*cos(x)*Ak1*Ak1*exp(-delta)*Phik1*IIk1 );
		dPhik2=dt*( D*(Ak1*exp(-delta)*IIk1) );
		Adphi = ( Ak1*exp(-delta)*Phik1 ); DAdphi=(D*Adphi);
		dIIk2= dt*( (-DAdphi[0])&( ( 2.0*(Adphi.sub(1,n-1))/(sin(xxx)*cos(xxx)) ) + (DAdphi.sub(1,n-1)) )&(3.0*DAdphi[n]) );		
		
		Ak2=A+0.5*dAk2;Phik2=Phi+0.5*dPhik2;IIk2=II+0.5*dIIk2;
		
		dAk3=dt*( (-2.0)*sin(x)*cos(x)*Ak2*Ak2*exp(-delta)*Phik2*IIk2 );
		dPhik3=dt*( D*(Ak2*exp(-delta)*IIk2) );
		Adphi = ( Ak2*exp(-delta)*Phik2 ); DAdphi=(D*Adphi);
		dIIk3= dt*((-DAdphi[0])&( ( 2.0*(Adphi.sub(1,n-1))/(sin(xxx)*cos(xxx)) ) + (DAdphi.sub(1,n-1)) )&(3.0*DAdphi[n]));		
		
		Ak3=A+dAk3;Phik3=Phi+dPhik3;IIk3=II+dIIk3;
		
		dAk4=dt*( (-2.0)*sin(x)*cos(x)*Ak3*Ak3*exp(-delta)*Phik3*IIk3 );
		dPhik4=dt*( D*(Ak3*exp(-delta)*IIk3) );
		Adphi = ( Ak3*exp(-delta)*Phik3 ); DAdphi=(D*Adphi);
		dIIk4= dt*((-DAdphi[0])&( ( 2.0*(Adphi.sub(1,n-1))/(sin(xxx)*cos(xxx)) ) + (DAdphi.sub(1,n-1)) )&(3.0*DAdphi[n]));		
		
		
		dA=((dAk1+dAk4)/6.0)+((dAk2+dAk3)/3.0);
		dPhi=((dPhik1+dPhik4)/6.0)+((dPhik2+dPhik3)/3.0);
		dII = (dIIk1+dIIk4)/6.0+(dIIk2+dIIk3)/3.0;
		
		dPhi[0]=0;dPhi[n]=0;dII[0]=0;dA[0]=0;dA[n]=0;
		
		A=A+dA;
		Phi=Phi+dPhi;
		II=II+dII;
		*/
		
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
		
		//dAnp1=( (-2.0)*sin(x)*cos(x)*A*A*exp(-delta)*Phi*II );
		//dPhinp1=( D*(A*exp(-delta)*II) );
		//Adphi = ( A*exp(-delta)*Phi ); DAdphi=(D*Adphi);
		//dIInp1 = (-DAdphi[0])&( ( 2.0*(Adphi.sub(1,n-1))/(sin(xxx)*cos(xxx)) ) + (DAdphi.sub(1,n-1)) )&(3.0*DAdphi[n]);

		//A= ( Anm1+ ((dt/12.0)*(5.0*dAnp1 -dAnm1+8.0*dA)) );//0.5*dt*( 3.0*dA - dAnm1 ) );//((dt/12.0)*(23.0*dA -16.0*dAnm1+5.0*dAnm2)) );
		//Phi= ( Phinm1+ ((dt/12.0)*(5.0*dPhinp1 -dPhinm1+8.0*dPhi)) );//0.5*dt*( 3.0*dPhi - dPhinm1 ));//( (dt/12.0)*(23.0*dPhi -16.0*dPhinm1+5.0*dPhinm2)) );
		//II= ( IInm1+ ((dt/12.0)*(5.0*dIInp1 -dIInm1+8.0*dII)) );//0.5*dt*( 3.0*dII - dIInm1 ));//((dt/12.0)*(23.0*dII -16.0*dIInm1+5.0*dIInm2)) );
		
		//ddelta=( -(sin(xx)*cos(xx)*( ((Phi*Phi)+(II*II)).sub(0,n-1) ) ) );
		//delta=((d|ddelta)&0);
		
		//Ap=(dd-ddelta);Ap.LUdcmp();a=((Ap|(ddelta))&0);A=( a + 1.0 );
		
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

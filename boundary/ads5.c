using namespace std;
#include "../Chebyshev.h"

#define L 1.5 
#define EPS 0.5
#define NSAVE 10.0
#define B 2.5
double phi0dot(double t){
		return -2.0*t*EPS*exp(-B*t*t);

}

int main(){
	int n=512;
	DoubleVector x=L*(1.0+ChebPoints(n));
	DoubleMatrix D=ChebDiff(n)/L,d=D.sub(0,n-1,0,n-1);
	DoubleVector xx=x.sub(0,n-1);
	d.LUdcmp();

	// Initial conditions
	DoubleVector Phi=0*x,II=0*x;
	
	DoubleVector delta,A,ddelta;DoubleMatrix Ap;
	delta=0*x;
	
	A=( delta + 1.0 );
	
	double dt=1.0e-5/NSAVE, tmax=100.0, t=-10.0;
	unsigned int nSave=0;
	
	cout << "#n = " << n << ", dt = " << dt << ", B = " << B << ", L = " << L << ", EPS = " << EPS << "\n";
	cout << "#t = " << t << "\n" << ((DoubleMatrix)x&Phi&II&delta&A);//a-aold);//

		
	// Time marching
	DoubleVector dII,dPhi,dA,dAnm1,dIInm1,dPhinm1,dAnm2,dIInm2,dPhinm2;//Time-marching
	DoubleVector Ad,AdII,Adphi,DAdphi;//temp vectors to speed up;
	
	dAnm2=Phi;dPhinm2=Phi;dIInm2=Phi;
	dAnm1=Phi;dPhinm1=Phi;dIInm1=Phi;
	
	double d4o3 = (4.0/3.0),ddto12=(dt/12.0),d2o3=(2.0/3.0);
	
	do{
		// Adams-Bashforth 3rd
		Ad=A*exp(-delta);AdII=Ad*II;
		dA=( d4o3*x*A*Phi*AdII );
		dPhi=( D*(AdII) );
		Adphi = ( Ad*Phi ); DAdphi=(D*Adphi);
		dII = (( (DAdphi.sub(0,n-1)) - 3.0*Adphi.sub(0,n-1)/xx )&0);
						
		A= ( A+ (ddto12*(23.0*dA -16.0*dAnm1+5.0*dAnm2)) );
		Phi= ( Phi+ ( ddto12*(23.0*dPhi -16.0*dPhinm1+5.0*dPhinm2)) );
		II= ( II+ (ddto12*(23.0*dII -16.0*dIInm1+5.0*dIInm2)) );
		
		//II[n]=phi0dot(t);//v3
		//II[n]=phi0dot(t);II[0]=0;//v1
		II[n]=phi0dot(t);Phi[0]=0;//v2
		
		ddelta=(d2o3*xx*( ((Phi*Phi)+(II*II)).sub(0,n-1) ) );
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
		
		//if(t<=0){	
			if(nSave==10000*NSAVE){
				cout << "\n#t = " << t << "\n" << ((DoubleMatrix)x&Phi&II&delta&A);
				nSave=0;
			}
		/*}else{
			if(nSave>=100*NSAVE){
				cout << "\n#t = " << t << "\n" << ((DoubleMatrix)x&Phi&II&delta&A);
				nSave=0;
			}
		}*/
	}while(t<tmax);
	
	return 0;
}

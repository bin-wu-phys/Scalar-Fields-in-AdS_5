using namespace std;
#include "../DoubleVector.h"
#define L 10.0
#define EPS 0.01
#define NSAVE 1.0
#define B 10 
/*
double phi0dot(double t){
	return EPS*sin(t);

}

double phi0dot(double t){
	return -2.0*t*EPS*exp(-B*t*t);

}

double phi0dot(double t){
	if(t<0)
		return EPS*exp(-B*pow(t,10));
	else
		return EPS;

}
*/

double phi0dot(double t){
	return EPS*exp(-pow(t/B,100));

}


int main(){
	int n=(int)(NSAVE*10000),i;double dx=L/n;
	DoubleVector x(n+1);for(i=0;i<=n;i++) x[i]=dx*i+0.005;

	// Initial conditions
	DoubleVector Phi=0*x,II=0*x;
	
	DoubleVector delta,A,ddelta;DoubleMatrix Ap;
	delta=0*x;
	
	A=( delta + 1.0 );
	
	double dt=0.1*dx, tmax=300.0, t=-1.1*B,tsave=0.5;
	unsigned int nSave=0,numsave=(int)(tsave/dt)+1;dt=tsave/numsave;
	
	cout << "#n = " << n << ", dt = " << dt << ", B = " << B << ", L = " << L << ", EPS = " << EPS << "\n";
	cout << "#t = " << t << "\n" << ((DoubleMatrix)x&Phi&II&delta&A);//a-aold);//

		
	// Time marching
	DoubleVector dII,dPhi(n+1),dA,dAnm1,dIInm1,dPhinm1,dAnm2,dIInm2,dPhinm2;
	dPhi[n]=0;
	DoubleVector Ad,AdII,Adphi,DAdphi(n+1);//temp vectors to speed up;
	
	dAnm2=Phi;dPhinm2=Phi;dIInm2=Phi;dAnm1=Phi;dPhinm1=Phi;dIInm1=Phi;

	double d4o3 = (4.0/3.0),ddto12=dt/12.0,d2o3=(2.0/3.0);

	double umin=L;
	
	do{
		// Adams-Bashforth 3rd
		Ad=A*exp(-delta);AdII=Ad*II;
		dA=( d4o3*x*A*Phi*AdII );
		for(i=0;i<n;i++)dPhi[i]=( AdII[i+1]-AdII[i] )/dx;
		Adphi = ( Ad*Phi ); for(i=1;i<=n;i++) DAdphi[i]=(Adphi[i]-Adphi[i-1])/dx;
		dII = (0&( (DAdphi.sub(1,n)) - 3.0*Adphi.sub(1,n)/x.sub(1,n) ));
		
		A= ( A+ (ddto12*(23.0*dA -16.0*dAnm1+5.0*dAnm2)) );
		Phi= ( Phi+ ( ddto12*(23.0*dPhi -16.0*dPhinm1+5.0*dPhinm2)) );
		II= ( II+ (ddto12*(23.0*dII -16.0*dIInm1+5.0*dIInm2)) );
		
		II[0]=phi0dot(t);Phi[n]=0;//v2
		
		
		ddelta=(d2o3*x.sub(1,n)*( ((Phi*Phi)+(II*II)).sub(1,n) ) );
		for(i=1;i<=n;i++)delta[i]=dx*ddelta[i-1]+delta[i-1];
		
		t+=dt;nSave++;dPhinm2=dPhinm1;dPhinm1=dPhi;dIInm2=dIInm1;dIInm1=dII;dAnm2=dAnm1;dAnm1=dA;
		
		for(int i=0;i<A.size();i++){
			if(A[i]<0.01&&x[i]<umin-0.1){
				umin=x[i];
				cout << "\n#t = " << t << "\n" << ((DoubleMatrix)x&Phi&II&delta&A);
				t=tmax;break;nSave=0;				
			}else if(isnan(A[i])){
				t=tmax;break;nSave=0;
			}

		}
		if(nSave==numsave){
			cout << "\n#t = " << t << "\n" << ((DoubleMatrix)x&Phi&II&delta&A);
			nSave=0;
		}
	}while(t<tmax);
	
	return 0;
}

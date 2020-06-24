using namespace std;
#include "Chebyshev.h"

#define L 5.0
#define EPS 10.0
//#define ZERO 1.0e-50
#define NSAVE 1.0

double phi0(double t){
        return EPS*exp(-B*t*t)/B;

}

double phi0dot(double t){
		return -2.0*t*EPS*exp(-10*t*t);
}

int main(){
	int n=10000;double dx=L/n;
	DoubleVector x(n+1);for(int i=0;i<=n;i++) x[i]=dx*i;
	
	// Initial conditions
	DoubleVector Phi=0*x,II=0*x,Phinm1=Phi,IInm1=Phi;//((2.0*EPS*exp(-4.0*(cot(x.sub(1,n))^2)/(PI*PI*SIG*SIG))/PI)&0);// Note: () for ^
	
	double dt=1.0e-5/NSAVE,mu=dt/dx, t=-10, tmax=10;//5000000000.0*NSAVE*dt;
	unsigned int nSave=0;
	
	cout << "#n = " << n << ", dt = " << dt << "\n\n";
	
	cout << "#t = " << t << "\n" << ((DoubleMatrix)x&Phi&II);//a-aold);//

		
	// Time marching
	do{
		II[0]=phi0dot(t);Phi[n]=0;
	
                phi[0]=phi0(t);

		for(int i=1;i<=n;i++){
			Phi[i]= ( Phinm1[i]+ mu*(IInm1[i]-IInm1[i-1]));
			II[i]= ( IInm1[i]+ mu*(Phinm1[i]-Phinm1[i-1]));
                        int im=i-1;
                        phi[i]=dx*Phi[im]+phi[im];
		}
	
		t+=dt;nSave++;Phinm1=Phi;IInm1=II;
		
		if(nSave==1000*NSAVE){
			//cout << ((DoubleMatrix)x&A&((A-AA)/A));
			cout << "\n#t = " << t << "\n" << ((DoubleMatrix)x&Phi&II);
			nSave=0;
		}
		//getchar();
	}while(t<tmax);

	return 0;
}

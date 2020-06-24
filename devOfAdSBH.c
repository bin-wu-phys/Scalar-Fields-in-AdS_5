/*
 *  a4fit.h
 *  
 *
 *  Created by Bin Wu on 7/8/12.
 *  Copyright 2012 Universitaet Bielefeld. All rights reserved.
 *
 */
#include <iostream>
#include <fstream>

using namespace std;
#include "Chebyshev.h"
#include "numericalAnalysis.h"

#define N		512
//#define N		256
#define EPS		0.5
#define ERR		1e-3
#define NT	
//#define NT		300
#define NT		50//400

//double a4=1/pow(0.987765,4.0);//4.188816127999566*EPS*EPS;
double a4=1/pow(0.988537,4.0);//400

double Abh(double u){
	return 1-a4*pow(u,4.0);
}

int main(){
	double pi,phi,delta,t=-0.5,dt=0.01,Amin,umin;
	//ifstream in("results/boundary/poincare/ads5L1n512e05a10.txt");
	ifstream in("results/boundary/poincare/ads5L1n512e05a400.txt");
	DoubleVector A(N+1),u(N+1);
	int lines=3 + NT*(N + 5);
	string line;
	
	for(int i=0;i<lines;i++){
		getline(in,line);
	}
	t+=NT*dt;
	
	do{
		Amin=1.0;umin=0;
		for(int i=0;i<=N;i++){
			in >> u[i] >> phi >> pi >> delta >> A[i];
			//cout << u[i] << "   " << phi << "   " << pi << "   " << delta << "   " << A[i] << endl;
		}
	
		for(int i=N;i>=0;i--){
			if(fabs(A[i]/Abh(u[i])-1)>ERR){
				if(i!=N){
					Amin=A[i+1];umin=u[i+1];
				}else {
					Amin = 1.0;umin=0;
				}
				cout <<t << "   " << umin << "   " << Amin << "\n";
				break;
			}/*
			if(Amin>A[i]){
				Amin=A[i];umin=u[i];
			}*/
		}
		cout <<t << "   " << umin << "   " << Amin << "\n";
		getline(in,line);getline(in,line);getline(in,line);getline(in,line);t+=dt;//cout << "\n\n";
		
	}while(t<2);
	
	return 0;
}
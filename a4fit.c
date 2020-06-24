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

#define NFIT	120
#define N		512
//#define NT		23//100
#define NT		1200//5

DoubleVector A(NFIT),u(NFIT);

double dE(double a4){
	return ((u^4)*(A-1+a4*(u^4))).sum();
}

int main(){
	double pi,phi,delta;
	//ifstream in("results/boundary/poincare/ads5L1n512e05a100.txt");
	//ifstream in("ads5L1p5n512e05b5.txt");
	ifstream in("ads5L2n512e05b1.txt");
	int lines=3 + NT*(N + 5);
	for(int i=0;i<lines;i++){
		string line;
		getline(in,line);
	}
	for(int i=0;i<(N-NFIT+1);i++){
		string line;
		getline(in,line);
	}
	
	for(int i=0;i<NFIT;i++){
		in >> u[i] >> phi >> pi >> delta >> A[i];
	}
	//cout << "\n\n";getchar();
	cout << u[0] << endl;
	cout << 1/pow(bisect(dE,0,10,1e-8),0.25) << endl;
	return 0;
}
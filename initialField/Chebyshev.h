/*
 *  Chebyshev.h
 *  
 *
 *  Created by Bin Wu on 5/27/12.
 *  Copyright 2012 Universitaet Bielefeld. All rights reserved.
 *
 */

#include "DoubleVector.h"

#define HX 1.0e-10

DoubleVector ChebPoints(const int n){
	DoubleVector v(n+1);
	
	for(int i=0;i<=n;i++) v[i]=cos(PI*i/n);//(1.0-2.0*i/n);//
	return v;
}

DoubleMatrix ChebDiff(const int n){
	DoubleMatrix d(n+1,n+1);DoubleVector x=ChebPoints(n);
	
	for(int i=0;i<=n;i++){
		for(int j=0;j<=n;j++){ 
			if(i==j){
				if (i==0) {
					d[0][0]=(2.0*n*n+1.0)/6.0;
				}else if(i==n){
					d[n][n]=-d[0][0];
				}
				else
					d[i][j]=-0.5*x[j]/(1.0-x[j]*x[j]);
			}else {
				double ci,cj;
				
				if(i==0||i==n)
					ci=2.0;
				else 
					ci=1.0;
				
				if(j==0||j==n)
					cj=2.0;
				else 
					cj=1.0;
				
				d[i][j]=ci*pow(-1.0,i+j)/(cj*(x[i]-x[j]));
			}
		}
	}

	return d;
}

DoubleMatrix Diff(const int n){
	DoubleMatrix d(n+1,n+1);DoubleVector x=ChebPoints(n);
	
	for(int i=0;i<=n;i++){
		for(int j=0;j<=n;j++){ 
			if(i==0){
				if(j==0){
					double h = x[i]-x[i+1];d[i][j]=1.0/h;
				}
				else if(j==1){ 
					double h = x[i]-x[i+1];d[i][j]=-1.0/h;
				}
				else d[i][j]=0;
			}
			else if(i==n){
				if(j==(n-1)){
					double h = x[i-1]-x[i];d[i][j]=1.0/h;
				}
				else if(j==n){
					double h = x[i-1]-x[i];d[i][j]=-1.0/h;
				}
				else 
					d[i][j]=0;
			}
			else {
				if(j==i-1){
					double h = x[i-1]-x[i+1];d[i][j]=1.0/h;
				}
				else if(j==i+1){
					double h = x[i-1]-x[i+1];d[i][j]=-1.0/h;
				}else {
					d[i][j]=0;
				}
			}
		}
	}
	return d;
}


double sfT(int n,double x){
	return cos(n*acos(x));
}
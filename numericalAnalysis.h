/*
 *  numericalAnalysis.h
 *  
 *
 *  Created by Bin Wu on 7/8/12.
 *  Copyright 2012 Universitaet Bielefeld. All rights reserved.
 *
 */

double bisect(double (*f)(double),double xbegin,double xend,double err){
	double fb,fe,xmid,fm;
	
	do{
		fb=(*f)(xbegin),fe=(*f)(xend),xmid=0.5*(xbegin+xend),fm=(*f)(xmid);
		
		if(fm*fb>0){
			xbegin=xmid;
		}else {
			xend=xmid;
		}
		
		cout << fb << "   " << fe << "   " << fm << endl;		
		
	}while (fabs(fm)>err);
	
	return xmid;
}	
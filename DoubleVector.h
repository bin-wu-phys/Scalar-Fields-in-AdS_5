/*
 *  DoubleVector.h
 *  
 *
 *  Created by Bin Wu on 5/25/12.
 *  Copyright 2012 Universitaet Bielefeld. All rights reserved.
 *
 */
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
//#define _MATLAB

//
//
// Vector
//
//
#define PI 3.141592653589793//23846 26433 83279 50288 41971 69399

class DoubleVector{
private:
	int _n;
	double *_v;

public:
	//Class definition
	DoubleVector();
	DoubleVector(int);
	DoubleVector(double v);
	DoubleVector(int n,double *v);
	DoubleVector(const DoubleVector& v);

	~DoubleVector();
	int size() const;

	DoubleVector sub(int rowb,int rowe); //Get subVector
	
	//Operators overloading
	double& operator[](int);
	const double& operator[](int) const;
	DoubleVector& operator=(const DoubleVector&);
	DoubleVector& operator=(const double);
	DoubleVector operator+(const DoubleVector&) const;
	DoubleVector operator-(const DoubleVector&) const;
	DoubleVector operator*(const DoubleVector&) const;
	DoubleVector operator/(const DoubleVector&) const;
	DoubleVector operator^(const double&) const;
	
	//Data analysis
	double sum();
};

DoubleVector::DoubleVector():_n(0),_v(NULL){}

DoubleVector::DoubleVector(int n){
	_v=new double[n];_n=n;
}

DoubleVector::DoubleVector(double v){
	_v=new double[1];_n=1;_v[0]=v;
}

DoubleVector::DoubleVector(int n,double *v){
	_v=new double[n];_n=n;for(int i=0;i<n;i++) _v[i]=v[i];
}

DoubleVector::DoubleVector(const DoubleVector& v){
	_n=v._n;
	_v=_n>0?new double[_n]:NULL;
	for(int i=0;i<_n;i++) _v[i]=v[i];
}

int DoubleVector::size() const{//const is also needed!
	return _n;
}

DoubleVector::~DoubleVector(){
	if(_v!=NULL) delete[] _v;
}

double& DoubleVector::operator[](int i){
	if(i<0||i>_n){
		cout << "Out of Bounds of DoubleVector\n";
		throw string("Out of Bounds of DoubleVector\n");
	}
	return _v[i];
}

const double& DoubleVector::operator[](int i) const{
	if(i<0||i>_n){
		cout << "Out of Bounds of DoubleVector\n";
		throw string("Out of Bounds of DoubleVector\n");
	}
	return _v[i];
}

DoubleVector& DoubleVector::operator=(const DoubleVector& v){
	if(this!=&v){
		if(_n!=v._n){
			if(_v!=NULL) delete[] _v;
			_n=v._n;
			_v=_n>0?new double[_n]:NULL;
		}
		for(int i=0;i<_n;i++)
			_v[i]=v[i];
	}
	return *this;
}

DoubleVector& DoubleVector::operator=(const double v){
	if(_n!=1){
		if(_v!=NULL) delete[] _v;
		_n=1;
		_v=_n>0?new double[_n]:NULL;
	}
		_v[0]=v;
	return *this;
}

DoubleVector DoubleVector::operator+(const DoubleVector& v) const{
	DoubleVector nrv(_n);
	if (_n!= v.size()) {
		cout << "***Adding two vectors with different dimensions!\n***";
		throw("Adding two vectors with different dimensions!\n");
	}
	else{
		for (int i=0; i<_n; i++) nrv[i]=_v[i]+v[i];
	}
	return nrv;
}

DoubleVector DoubleVector::operator-(const DoubleVector& v) const{
	DoubleVector nrv(_n);
	if (_n!= v.size()) {
		cout << "***Substraction of two vectors with different dimensions!\n***";
		throw("Substraction of two vectors with different dimensions!\n");
	}
	else{
		for (int i=0; i<_n; i++) nrv[i]=_v[i]-v[i];
	}
	return nrv;
}

DoubleVector DoubleVector::operator*(const DoubleVector& v) const{
	DoubleVector nrv(_n);
	if (_n!= v.size()) {
		cout << "***Multiplication of two vectors with different dimensions!\n***";
		throw("Multiplication of two vectors with different dimensions!\n");
	}
	else{
		for (int i=0; i<_n; i++) nrv[i]=_v[i]*v[i];
	}
	return nrv;
}

DoubleVector DoubleVector::operator/(const DoubleVector& v) const{
	DoubleVector nrv(_n);
	if (_n!= v.size()) {
		cout << "***Division of two vectors with different dimensions!\n***";
		throw("Division of two vectors with different dimensions!\n");
	}
	else{
		for (int i=0; i<_n; i++) nrv[i]=_v[i]/v[i];
	}
	return nrv;
}

DoubleVector DoubleVector::operator^(const double& v) const{
	DoubleVector nrv(_n);
	for (int i=0; i<_n; i++) nrv[i]=pow(_v[i],v);
	return nrv;
}


DoubleVector DoubleVector::sub(int rowb,int rowe){
	DoubleVector subVector(rowe-rowb+1);
	if(rowe<_n&&rowb>=0){
		int k=0;
		for(int i=rowb;i<=rowe;i++){
			subVector[k]=_v[i];k++;
		}
	}else{
		cout << "***Out of bounds in subVector***\n";
		throw("***Out of bounds in subVector***\n");
	}
	return subVector;	
}

//Data analysis
double DoubleVector::sum(){
	double sum=0;
	for(int i=0;i<_n;i++) sum+=_v[i];
	return sum;
}

//
//Top-level functions
//

//Output
ostream& operator<<(ostream& out,const DoubleVector& v){
	out.precision(std::numeric_limits<double>::digits10);
	out << setiosflags(ios::scientific | ios::showpoint);
	if(v.size()>0){
#ifdef _MATLAB
		int nm1=v.size()-1;
		out <<" [ ";
		for(int i=0;i<nm1;i++) out << v[i] <<"; ";
		out << v[nm1] << " ] ";
#else
		out << "\n";
		for(int i=0;i<v.size();i++){
			out << v[i] << "\n";
		}
		out << "\n";
#endif
	}else {
		out << "NULL";
	}
	
	return out;
}

DoubleVector operator&(const DoubleVector& v1,const DoubleVector& v2){
	
	if(v1.size()>0&&v2.size()>0){
		DoubleVector nrv(v1.size()+v2.size());
		for(int j=0;j<v1.size();j++) nrv[j]=v1[j];
		for(int k=0;k<v2.size();k++) nrv[v1.size()+k]=v2[k];
		return nrv;
	}
	else {
		if(v1.size()<=0) return v2;
		else if(v2.size()<=0) return v1;
		else{
			return DoubleVector(0);
		}
	}	
}


DoubleVector operator&(const double& v1,const DoubleVector& v2){
	
	if(v2.size()>0){
		DoubleVector nrv(1+v2.size());
		nrv[0]=v1;
		for(int k=0;k<v2.size();k++) nrv[1+k]=v2[k];
		return nrv;
	}
	else {
		return DoubleVector(v1);
	}	
}

DoubleVector operator&(const DoubleVector& v1,const double& v2){
	
	if(v1.size()>0){
		DoubleVector nrv(1+v1.size());
		for(int k=0;k<v1.size();k++) nrv[k]=v1[k];
		nrv[v1.size()]=v2;
		return nrv;
	}
	else {
		return DoubleVector(v2);
	}	
}

// Top-level operator overloading

DoubleVector operator+(const double& dv,const DoubleVector& v){
	DoubleVector nrv(v.size());
	for (int i=0; i<v.size(); i++) nrv[i]=dv+v[i];
	return nrv;	
}

DoubleVector operator+(const DoubleVector& v, const double& dv){
	DoubleVector nrv(v.size());
	for (int i=0; i<v.size(); i++) nrv[i]=dv+v[i];
	return nrv;	
}

DoubleVector operator-(const double& dv,const DoubleVector& v){
	DoubleVector nrv(v.size());
	for (int i=0; i<v.size(); i++) nrv[i]=dv-v[i];
	return nrv;	
}

DoubleVector operator-(const DoubleVector& v){
	DoubleVector nrv(v.size());
	for (int i=0; i<v.size(); i++) nrv[i]=-v[i];
	return nrv;	
}

DoubleVector operator-(const DoubleVector& v, const double& dv){
	DoubleVector nrv(v.size());
	for (int i=0; i<v.size(); i++) nrv[i]=v[i]-dv;
	return nrv;	
}

DoubleVector operator*(const double& dv,const DoubleVector& v){
	DoubleVector nrv(v.size());
	for (int i=0; i<v.size(); i++) nrv[i]=dv*v[i];
	return nrv;	
}

DoubleVector operator/(const double& dv,const DoubleVector& v){
	DoubleVector nrv(v.size());
	for (int i=0; i<v.size(); i++) nrv[i]=dv/v[i];
	return nrv;	
}

DoubleVector operator/(const DoubleVector& v, const double& dv){
	DoubleVector nrv(v.size());
	for (int i=0; i<v.size(); i++) nrv[i]=v[i]/dv;
	return nrv;	
}

// Mathematical functions overloading

DoubleVector sin(const DoubleVector& v){
	DoubleVector retV(v.size());
	for(int i=0;i<v.size();i++){
		if(v[i]==PI)
			retV[i]=0;
		else
			retV[i]=sin(v[i]);
	}
	return retV;
}

DoubleVector exp(const DoubleVector& v){
	DoubleVector retV(v.size());
	for(int i=0;i<v.size();i++){ 
		retV[i]=exp(v[i]);
	}
	return retV;
}

DoubleVector log(const DoubleVector& v){
	DoubleVector retV(v.size());
	for(int i=0;i<v.size();i++){ 
		retV[i]=log(v[i]);
	}
	return retV;
}

DoubleVector pow(const DoubleVector& v,const double &ex){
	DoubleVector retV(v.size());
	for(int i=0;i<v.size();i++){ 
		retV[i]=pow(v[i],ex);
	}
	return retV;
}

DoubleVector cos(const DoubleVector& v){
	DoubleVector retV(v.size());
	for(int i=0;i<v.size();i++){ 
		if(v[i]==0.5*PI)
			retV[i]=0;
		else
			retV[i]=cos(v[i]);
	}
	return retV;
}

DoubleVector tan(const DoubleVector& v){
	DoubleVector retV(v.size());
	for(int i=0;i<v.size();i++){ 
		retV[i]=tan(v[i]);
	}
	return retV;
}

DoubleVector cot(const DoubleVector& v){
	DoubleVector retV(v.size());
	for(int i=0;i<v.size();i++){ 
		if(v[i]==0.5*PI)
			retV[i]=0;
		else
			retV[i]=1/tan(v[i]);
	}
	return retV;
}

DoubleVector asin(const DoubleVector& v){
	DoubleVector retV(v.size());
	for(int i=0;i<v.size();i++){ 
		retV[i]=asin(v[i]);
	}
	return retV;
}

DoubleVector acos(const DoubleVector& v){
	DoubleVector retV(v.size());
	for(int i=0;i<v.size();i++){ 
		retV[i]=acos(v[i]);
	}
	return retV;
}

//
// Matrix
//

class DoubleMatrix{
private:
	int _row,_col,*_indx;
	double **_v,**_lu,_d;
public:
	DoubleMatrix();
	DoubleMatrix(int m, int n);
	DoubleMatrix(const DoubleVector&);
	DoubleMatrix(const DoubleMatrix&);
	inline int row() const;
	inline int col() const;
	~DoubleMatrix();
	
	DoubleMatrix sub(int rowb,int rowe,int colb,int cole);
	
	//Operators overloading
	inline double* operator[](const int i);
	inline const double* operator[](const int i) const;
	DoubleMatrix& operator=(const DoubleMatrix&);
	DoubleMatrix& operator=(const DoubleVector&);
	DoubleMatrix operator+(const DoubleMatrix&) const;
	DoubleMatrix operator-(const DoubleMatrix&) const;
	DoubleMatrix operator*(const DoubleMatrix&) const;
	DoubleVector operator*(const DoubleVector&) const;
	DoubleMatrix operator/(const DoubleMatrix&) const;
	DoubleMatrix operator^(const double&) const;
	
	//LU Decomposition
	void LUdcmp();
	double det();
	DoubleVector operator|(const DoubleVector&) const;
	DoubleMatrix operator~() const;
};

DoubleMatrix::DoubleMatrix():_row(0),_col(0),_d(1.0),_v(NULL),_lu(NULL),_indx(NULL){}

DoubleMatrix::DoubleMatrix(int m, int n):_row(m),_col(n),_d(1.0),_lu(NULL),_indx(NULL),_v(m>0?new double*[m]:NULL){
	for(int i=0;i<_row;i++)
		_v[i]=_col>0?new double[_col]:NULL;
}

DoubleMatrix::DoubleMatrix(const DoubleMatrix& dm):_d(1.0),_indx(NULL){
	_row=dm.row();_col=dm.col();_lu=NULL;
	
	_v=_row>0?new double*[_row]:NULL;
	for(int i=0;i<_row;i++)
		_v[i]=_col>0?new double[_col]:NULL;
	
	for(int i=0;i<_row;i++)
		for(int j=0;j<_col;j++)
			_v[i][j]=dm[i][j];
}


DoubleMatrix::DoubleMatrix(const DoubleVector& v):_d(1.0),_indx(NULL){
	// Allocate new memory
	_row=v.size();_col=1;_lu=NULL;
		
	_v=_row>0?new double*[_row]:NULL;
	for(int i=0;i<_row;i++)_v[i]=new double[1];
	
	for(int i=0;i<_row;i++)
		_v[i][0]=v[i];
}

inline int DoubleMatrix::row() const{
	return _row;
}

inline int DoubleMatrix::col() const{
	return _col;
}

DoubleMatrix::~DoubleMatrix(){
	
	if(_lu!=NULL){
		for(int i=0;i<_row;i++) if(_lu[i]!=NULL) delete [] _lu[i];
		delete [] _lu;
	}
	
	if(_indx!=NULL) delete [] _indx;
	
	if(_v!=NULL){
		for(int i=0;i<_row;i++) if(_v[i]!=NULL) delete [] _v[i];
		delete [] _v;
	}
}

inline double* DoubleMatrix::operator[](const int i){
	if(i<0||i>=_row){
		cout << "***DoubleMatrix row out of bounds***\n";
		throw("DoubleMatrix row out of bounds");
	}
	return _v[i];
}

inline const double* DoubleMatrix::operator[](const int i) const{
	if(i<0||i>=_row){
		cout << "***DoubleMatrix row out of bounds***\n";
		throw("DoubleMatrix row out of bounds");
	}
	return _v[i];
}

DoubleMatrix DoubleMatrix::sub(int rowb,int rowe,int colb,int cole){
	if(rowe<_row&&cole<_col&&rowb>=0&&colb>=0){
		DoubleMatrix submatrix(rowe-rowb+1,cole-colb+1);
		int k=0;
		for(int i=rowb;i<=rowe;i++){
			int l=0;
			for(int j=colb;j<=cole;j++){
				submatrix[k][l]=_v[i][j];l++;
			}
			k++;
		}
		return submatrix;
	}
	else {
		cout << "***Out of bounds in subMatrix***\n";	
		throw("***Out of bounds in subMatrix***\n");
		return DoubleMatrix(0,0);
	}

}

DoubleMatrix& DoubleMatrix::operator=(const DoubleMatrix& v){
	if(this!=&v){
		if(_row!=v._row){
			// Free old memory
			if(_v!=NULL){
				for(int i=0;i<_row;i++) if(_v[i]!=NULL) delete [] _v[i];
				delete [] _v;
			}
			
			// Allocate new memory
			_row=v.row();_col=v.col();
			
			_v=_row>0?new double*[_row]:NULL;
			for(int i=0;i<_row;i++)
				_v[i]=_col>0?new double[_col]:NULL;
		}
		else if(_col!=v._col) {
			for(int i=0;i<_row;i++) if(_v[i]!=NULL) delete [] _v[i];
			_col=v._col;for(int i=0;i<_row;i++) _v[i]=_col>0?new double[_col]:NULL;
		}
		
		//cout << "***Start to release LU***\n";
		if(_lu!=NULL){
			for(int i=0;i<_row;i++) if(_lu[i]!=NULL) delete [] _lu[i];
			delete [] _lu;
		}
		
		if(_indx!=NULL) delete [] _indx;	
		
		_lu=NULL;_indx=NULL;

		//cout << "***End to release LU***\n";
		
		for(int i=0;i<_row;i++)
			for(int j=0;j<_col;j++)
				_v[i][j]=v[i][j];
	}
	return *this;
}

DoubleMatrix& DoubleMatrix::operator=(const DoubleVector& v){
	if(_row!=v.size()){
		// Free old memory
		if(_v!=NULL){
			for(int i=0;i<_row;i++) if(_v[i]!=NULL) delete [] _v[i];
			delete [] _v;
		}

		// Allocate new memory
		_row=v.size();_col=1;
			
		_v=_row>0?new double*[_row]:NULL;
		for(int i=0;i<_row;i++)
			_v[i]=new double[1];
	}
	else if(_col!=1) {
		for(int i=0;i<_row;i++) if(_v[i]!=NULL) delete [] _v[i];
		_col=1;for(int i=0;i<_row;i++) _v[i]=new double[_col];
	}		
	
	if(_lu!=NULL){
		for(int i=0;i<_row;i++) if(_lu[i]!=NULL) delete [] _lu[i];
		delete [] _lu;
	}
	
	if(_indx!=NULL) delete [] _indx;	

	_lu=NULL;_indx=NULL;
	
	for(int i=0;i<_row;i++)
		_v[i][0]=v[i];
	return *this;
}

DoubleMatrix DoubleMatrix::operator+(const DoubleMatrix& v) const{
	DoubleMatrix nrv(_row,_col);
	if (_row!= v.row()||_col!= v.col()) {
		cout << "***Adding two Matrices with different dimensions***\n";
		throw("Adding two Matrices with different dimensions!\n");
	}
	else{
		for (int i=0; i<_row; i++) for (int j=0; j<_col; j++) nrv[i][j]=_v[i][j]+v[i][j];
	}
	return nrv;
}

DoubleMatrix DoubleMatrix::operator-(const DoubleMatrix& v) const{
	DoubleMatrix nrv(_row,_col);
	if (_row!= v.row()||_col!= v.col()) {
		cout << "***Substraction of two Matrices with different dimensions***\n";
		throw("Substraction of two Matrices with different dimensions!\n");
	}
	else{
		for (int i=0; i<_row; i++) for (int j=0; j<_col; j++) nrv[i][j]=_v[i][j]-v[i][j];
	}
	return nrv;
}

DoubleMatrix DoubleMatrix::operator*(const DoubleMatrix& v) const{
	DoubleMatrix nrv(_row,v.col());
	if (_col!= v.row()) {
		cout << "***Mutiplication of two Matrices with different dimensions***\n";
		throw("Mutiplication of two Matrices with different dimensions!\n");
	}
	else{
		for (int i=0; i<_row; i++){ 
			for (int j=0; j<v.col(); j++){
				nrv[i][j]=0;
				for(int k=0;k<_col;k++)
					nrv[i][j]+=(_v[i][k]*v[k][j]);
			}
		}
	}
	return nrv;
}

DoubleVector DoubleMatrix::operator*(const DoubleVector& v) const{
	DoubleVector nrv(_row);
	if (_col!= v.size()) {
		cout << "***Mutiplication of a Matrix and a Vector with different dimensions***\n";
		throw("Mutiplication of a Matrix and a Vector with different dimensions!\n");
	}
	else{
		for (int i=0; i<_row; i++){ 
				nrv[i]=0;
				for(int k=0;k<_col;k++)
					nrv[i]+=(_v[i][k]*v[k]);
		}
	}
	return nrv;
}

//Top-level functions
//Output
ostream& operator<<(ostream& out,const DoubleMatrix& m){
        out.precision(std::numeric_limits<double>::digits10);
        out << setiosflags(ios::scientific |
ios::showpoint);
	if(m.row()>0&&m.col()>0){
#ifdef _MATLAB
		int nrow=m.row()-1, ncol=m.col()-1;
		out <<" [";
		for(int i=0;i<nrow;i++){
			//cout << "{ ";
			for(int j=0;j<ncol;j++)
				out << m[i][j] <<", ";
			out << m[i][ncol] << " ; ";
		}
		//cout << "{ ";
		for(int j=0;j<ncol;j++)
			out << m[nrow][j] <<", ";
		out << m[nrow][ncol] << " ] ";
#else
		out<<"\n";
		for(int i=0;i<m.row();i++){
			for(int j=0;j<m.col();j++){
				out << m[i][j] << "   ";
			}
			out << "\n";		
		}
		out << "\n";
#endif
	}else {
		out << "NULL";
	}
	return out;
}

DoubleMatrix operator*(const double& dv,const DoubleMatrix& v){
	DoubleMatrix nrv(v.row(),v.col());
	for (int i=0; i<v.row(); i++) for(int j=0;j<v.col();j++) nrv[i][j]=dv*v[i][j];
	return nrv;	
}

DoubleMatrix operator*(const DoubleVector& dv,const DoubleMatrix& v){
	DoubleMatrix nrv(v.row(),v.col());
	for (int i=0; i<v.row(); i++) for(int j=0;j<v.col();j++) nrv[i][j]=dv[i]*v[i][j];
	return nrv;	
}

DoubleMatrix operator+(const DoubleMatrix& m,const DoubleVector& v){
	DoubleMatrix nrv=m;
	if(m.row()==m.col()&&m.row()==v.size()){
		for (int i=0; i<m.row(); i++) nrv[i][i]+=v[i];
	}else {
		cout << "***Addition of Matrix and Vector with different dimensions***\n";
		throw("***Addition of Matrix and Vector with different dimensions***\n");
	}

	return nrv;	
}

DoubleMatrix operator+(const DoubleVector& v,const DoubleMatrix& m){
	DoubleMatrix nrv=m;
	if(m.row()==m.col()&&m.row()==v.size()){
		for (int i=0; i<m.row(); i++) nrv[i][i]+=v[i];
	}else {
		cout << "***Addition of Matrix and Vector with different dimensions***\n";
		throw("***Addition of Matrix and Vector with different dimensions***\n");
	}
	
	return nrv;	
}


DoubleMatrix operator-(const DoubleMatrix& m){
	DoubleMatrix nrv(m.row(),m.col());
	for (int i=0; i<m.row(); i++) for(int j=0;j<m.col();j++) nrv[i][j]=-m[i][j];
	return nrv;	
}

DoubleMatrix operator-(const DoubleMatrix& m,const DoubleVector& v){
	DoubleMatrix nrv=m;
	if(m.row()==m.col()&&m.row()==v.size()){
		for (int i=0; i<m.row(); i++) nrv[i][i]-=v[i];
	}else {
		cout << "***Subtraction of Matrix and Vector with different dimensions***\n";
		throw("***Subtraction of Matrix and Vector with different dimensions***\n");
	}
	
	return nrv;	
}

DoubleMatrix operator-(const DoubleVector& v,const DoubleMatrix& m){
	DoubleMatrix nrv=(-m);
	if(m.row()==m.col()&&m.row()==v.size()){
		for (int i=0; i<m.row(); i++) nrv[i][i]+=v[i];
	}else {
		cout << "***Subtraction of Matrix and Vector with different dimensions***\n";
		throw("***Subtraction of Matrix and Vector with different dimensions***\n");
	}
	
	return nrv;	
}

DoubleMatrix operator/(const DoubleMatrix& v,const double& dv){
	DoubleMatrix nrv(v.row(),v.col());
	for (int i=0; i<v.row(); i++) for(int j=0;j<v.col();j++) nrv[i][j]=v[i][j]/dv;
	return nrv;	
}

DoubleMatrix operator&(const DoubleMatrix& v1,const DoubleMatrix& v2){
	
	if(v1.row()==v2.row()){
		DoubleMatrix nrv(v1.row(),v1.col()+v2.col());
		for(int i=0; i<v1.row(); i++){
			for(int j=0;j<v1.col();j++) 
				nrv[i][j]=v1[i][j];
			for(int k=0;k<v2.col();k++)
				nrv[i][v1.col()+k]=v2[i][k];
		}
		return nrv;
	}
	else {
		if(v1.row()<=0) return v2;
		else if(v2.row()<=0) return v1;
		else{
			cout << "***Trying to combine Matrices with different dimensions***\n";
			return DoubleMatrix(0,0);
		}
	}	
}

#ifdef _MATLAB
void OutputList(const DoubleMatrix& m){
	for(int i=0;i<m.row();i++){
		for(int j=0;j<m.col();j++){
			cout << m[i][j] << "   ";
		}
		cout << "\n";		
	}
}
#endif

// LU Decomposition from numerical recipies 3rd

void DoubleMatrix::LUdcmp(){
	if(_row==_col){
		if(_lu==NULL){
			
			const double TINY=1.0e-40;
			int i,imax,j,k,n=_row;
			double big,temp;
			DoubleVector vv(n);
			
			_lu=_row>0?new double*[_row]:NULL;
			for(i=0;i<_row;i++)
				_lu[i]=_col>0?new double[_col]:NULL;
			
			for(i=0;i<_row;i++)
				for(j=0;j<_col;j++)
					_lu[i][j]=_v[i][j];	
			
			_indx=new int[_row];

			_d=1.0;
			for (i=0;i<n;i++) {
				big=0.0;
				for (j=0;j<n;j++)
					if ((temp=abs(_lu[i][j])) > big) big=temp;
				if (big == 0.0){
					cout << "***Singular matrix in LUdcmp***\n";
					throw("Singular matrix in LUdcmp");
				}
				vv[i]=1.0/big;
			}
			for (k=0;k<n;k++) {
				big=0.0;
				for (i=k;i<n;i++) {
					temp=vv[i]*abs(_lu[i][k]);
					if (temp > big) {
						big=temp;
						imax=i;
					}
				}
				if (k != imax) {
					for (j=0;j<n;j++) {
						temp=_lu[imax][j];
						_lu[imax][j]=_lu[k][j];
						_lu[k][j]=temp;
					}
					_d = -_d;
					vv[imax]=vv[k];
				}
				_indx[k]=imax;
				if (_lu[k][k] == 0.0) _lu[k][k]=TINY;
				for (i=k+1;i<n;i++) {
					temp=_lu[i][k] /= _lu[k][k];
					for (j=k+1;j<n;j++)
						_lu[i][j] -= temp*_lu[k][j];
				}
			}
		}
	}
}

DoubleVector DoubleMatrix::operator|(const DoubleVector& b) const
{
	int i,ii=0,ip,j;
	double sum;int n=_row;
	DoubleVector x(_row);
	
	if (b.size() != _row){
		cout << "***LUdcmp::solve bad sizes***\n";
		throw("LUdcmp::solve bad sizes");
	}
	
	for (i=0;i<n;i++) x[i] = b[i];
	for (i=0;i<n;i++) {
		ip=_indx[i];
		sum=x[ip];
		x[ip]=x[i];
		if (ii != 0)
			for (j=ii-1;j<i;j++) sum -= _lu[i][j]*x[j];
		else if (sum != 0.0)
			ii=i+1;
		x[i]=sum;
	}
	for (i=n-1;i>=0;i--) {
		sum=x[i];
		for (j=i+1;j<n;j++) sum -= _lu[i][j]*x[j];
		x[i]=sum/_lu[i][i];
	}
	
	return x;
}

double DoubleMatrix::det()
{
	double dd = _d;
	for (int i=0;i<_row;i++) dd *= _lu[i][i];
	return dd;
}

DoubleMatrix DoubleMatrix::operator~() const{
	DoubleVector x(_row),y;
	DoubleMatrix invm(_row,_col);
	
	for(int i=0;i<_col;i++){
		for(int j=0;j<_row;j++) x[j]=0;x[i]=1;
		y=(this->operator|(x));
		for(int j=0;j<_row;j++) invm[j][i]=y[j];
	}
	
	return invm;
}

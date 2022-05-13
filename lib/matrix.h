// **************************************************************************************
// matrix.h
//		This file describes certain matri operations for the GRPulse code.  
// **************************************************************************************


#ifndef MATRIX
#define MATRIX

#include <stdio.h>

//calculate the determinant of an NxN matrix
//uses LU decomposition, based on GSL function
//will destroy the matrix m
template <size_t N>
double determinant(double (&m)[N][N]){
	double det=1.0, temp = 0.0, big = 0.0, ajj=0.0, coeff=0.0;
	int i=0,j=0,k=0,ipiv=0;
	
	for(j=0; j<N-1;j++){
		//find the max in this column
		big = 0.0;
		for(i=j;i<N;i++){
			temp = fabs(m[i][j]);
			if(temp>big){
				big=temp;
				ipiv=i;
			}
		}		
		//if max is not on diagonal, swap rows
		if(ipiv != j){
			for(k=0;k<N;k++){
				temp = m[j][k];
				m[j][k] = m[ipiv][k];
				m[ipiv][k]=temp;
			}
			//swapping rows swaps sign of det
			det = -det;
		}
		ajj = m[j][j];
		if(ajj!=0.0){
			for(i=j+1;i<N;i++){
				coeff = m[i][j]/ajj;
				for(k=j+1;k<N;k++){
					m[i][k] -= coeff*m[j][k];
				}
			}	
		}
		det *= m[j][j];	
	}
	det *= m[N-1][N-1];
	return det;
}

//in equation Ax=b, invert NxN matrix A to solve for x
//will destroy the matrix A
//can handle equations Ax=0
template <size_t N>
int invertMatrix(double (&A)[N][N], double(&b)[N], double(&x)[N]){
//void invertMatrix(int N, double** A, double *b, double *x){
	//a flag, in case we are soving homoegenous problem
	bool HOMOGENEOUS = true;
	double dummy = 0.0;
	for(int j=0; j<N; j++) HOMOGENEOUS &= (b[j]==0.0);
	//reduce to row-echelon form
	int L=0, indxBottom = N-1;	

	for(int i=0;i<N;i++){
		L=0;
		//make sure we are not on a zero element
		while(fabs(A[i][L])==0.0){
			L++;
			//if the entire row is zero
			if(L>=N){//swap to bottom
				L=0;
				dummy = b[indxBottom];
				b[indxBottom] = b[i];
				b[i] = dummy;
				for(int j=0; j<N; j++){
					dummy = A[indxBottom][j];
					A[indxBottom][j] = A[i][j];
					A[i][j] = dummy;
				}
				//switch more rows to bottom as time goes on
				indxBottom--;
				if(indxBottom <= 0) {printf("I think your whole matrix is zero...\n"); return 1;}
				break;
			}
		}
		//if pivot is too small, can cause a problem
		//so find a row where the pivot isn't too small and swap
		if(fabs(A[i][L])<1.0e-10){
			int K=i;
			while(fabs(A[K][L])<1.0e-10){
				K++;
				if(K>=N){
					K=i;
					break;	
				}
			}
			dummy = b[K];
			b[K]  = b[i];
			b[i]  = dummy;
			for(int j=0;j<N;j++){
				dummy   = A[K][j];
				A[K][j] = A[i][j];
				A[i][j] = dummy;
			}
		}
		//now that we have nonzero row, divide everthing by front element
		double a = 1.0/A[i][L];
		b[i] *= a;
		for(int j=0; j<N; j++){
			A[i][j] *= a;
		}
		A[i][L] = 1.0;//set to one, for good measure
		//now that we have a leading one, remove column elements below
		double lead = 0.0;
		for(int j=i+1; j<N; j++){
			lead = A[j][L];
			for(int k=0; k<N;k++){
				A[j][k] -= lead*A[i][k]; //subtract row i * lead
			}
			b[j] -= lead*b[i];
			A[j][L] = 0.0;//set to zero, for good measure
		}
	}
	//rearrange rows to make sure in proper order
	int roworder[N];
	for(int i=0;i<N;i++){
		int L=0;
		while(A[i][L]==0.0 && L<N){
			L++;
		}
		roworder[i] = L;
	}
	for(int i=0; i<N; i++){
		for(int j=i; j<N; j++){
			if(roworder[i]>roworder[j] && roworder[j]!=N){
				dummy = b[i];
				b[i] = b[j];
				b[j] = dummy;
				for(int k=0;k<N;k++){
					dummy   = A[i][k];
					A[i][k] = A[j][k];
					A[j][k] = dummy;
				}		
			}
		}
	}
	
	//resubstitute
	x[N-1] = b[N-1];
	if(HOMOGENEOUS) x[N-1] = 1.0;
	for(int i=N-2;i>=0; i--){
		L=0;
		while(A[i][L]==0 && L<N){	
			L++;
		}
		x[i] = (HOMOGENEOUS ? 0.0 : b[i]);
		for(int j=L+1; j<N; j++){
			x[i] -= A[i][j]*x[j];
		}
	}	
	return 0;
}

#endif
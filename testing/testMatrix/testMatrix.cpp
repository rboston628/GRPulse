// This program is meant to test the matrix inversion method.
// First, it will populate a matrix with random integers, based on a very (very very) simple
//   psuedo-random number  generator, as well as a vector, and then solve the system
//    A*x = b for the variable x
// The output for MATRIX and VECTOR are meant to be entered into the accompanying Mathematica notebook, 
// 		so that the resultrs can be checked
#include <stdio.h>
#include <math.h>
#include "../../lib/matrix.h"

int main(){
	const int a=53, b=122, r=17737;
	int rando = 39;
	
	printf("SEED = %d\n", rando); 
	const int N = 16;


	double myMatrix[N][N];
	printf("MATRIX:\n");
	for(int i=0; i<N; i++){
		printf("\t{");
		rando = (a*rando+b)%r;
		myMatrix[i][0] = double(rando);
		printf("%0.0lf", myMatrix[i][0]);
		for(int j=1; j<N; j++){
			rando = (a*rando+b)%r;
			myMatrix[i][j] = double(rando);
			printf(", %0.0lf", myMatrix[i][j]);
		}
		printf(" }%c\n", (i<N-1?',':' '));
	}

	double myVector[N];
	printf("VECTOR:\n");
	printf("\t{");
	for(int i=0; i<N; i++){
		rando = (a*rando+b)%r;
		myVector[i] = double(rando);
		printf("%0.0lf%c", myVector[i], (i<N-1?',':' '));
	}
	printf("}\n");

	double mySolution[N];
	invertMatrix(myMatrix, myVector, mySolution);

	printf("SOLUTION:\n");
	printf("\t[");
	for(int i=0; i<N; i++){
		printf(" %lf", mySolution[i]);
	}
	printf(" ]\n");


	rando = 0xDEADBEEF;
	printf("SEED = %X\n", rando); 
	printf("MATRIX:\n");
	for(int i=0; i<N; i++){
		printf("\t{");
		rando = (a*rando+b)%r;
		myMatrix[i][0] = double(rando);
		printf(" %0.0lf", myMatrix[i][0]);
		for(int j=1; j<N; j++){
			rando = (a*rando+b)%r;
			myMatrix[i][j] = double(rando);
			printf(", %0.0lf", myMatrix[i][j]);
		}
		printf(" }%c\n", (i<N-1?',':' '));
	}

	printf("VECTOR:\n");
	printf("\t{");
	for(int i=0; i<N; i++){
		rando = (a*rando+b)%r;
		myVector[i] = double(rando);
		printf("%0.0lf%c", myVector[i], (i<N-1?',':' '));
	}
	printf(" }\n");
	invertMatrix(myMatrix, myVector, mySolution);

	printf("SOLUTION:\n");
	printf("\t[");
	for(int i=0; i<N; i++){
		printf(" %lf", mySolution[i]);
	}
	printf(" ]\n");


	return 0;
}
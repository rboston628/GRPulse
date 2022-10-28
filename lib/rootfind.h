#ifndef ROOTFINDINGMETHODS
#define ROOTFINDINGMETHODS

#include <functional>

double pseudo_unif();

// *************************************************************************************
//					BISECTION METHODS
//  These methods are for implementing bisection searches in one parameter
// *************************************************************************************

// this will find brackets using a simple expansion/movement search
// one bracket is moved to one side until the function changes sign, to find brackets
void bisection_find_brackets_move(
	std::function<double(double)> func,		//the function to find zero of
	double x0,								//an initial guess for the zero
	double& xmin,							//the lower bracket -- will be returned
	double& xmax							//the upper bracket -- will be returned
);

// this will find brackets using newton's method to look for the next zero, then a bit beyond it
// why use one zero-finding method to prepare another zero-finding method? 
//    because Newton's method can fail, but bisection searches can go to nearly arbitrary accuracy
void bisection_find_brackets_newton(
	std::function<double(double)> func,		//the function to find zero of
	double x,								//an initial guess for the zero
	double& xmin,							//the lower bracket -- will be returned
	double& xmax							//the upper bracket -- will be returned
);

//given brackets bounding a single zero, find the zero
// the brackets xmin, xmax MUST bound a single zero
double bisection_search(
	std::function<double(double)> func,		//the function to find zero of
	double &x,								//the location of zero -- will be returned
	double xmin,							//the lower bracket
	double xmax								//the upper bracket
);


// *************************************************************************************
//					NEWTON METHODS
//  These methods are for gradient-descent methods
// *************************************************************************************

template <typename T>
T gen_abs(T x);

template<size_t np>
std::function<bool(double[np])> no_limit = [](double x[np])->bool{return true;};

//single-variable search, searching for ZERO of func
//intended to work with real-valued or complex-valued numbers
template <typename T>
T newton_search(
	std::function<T(T)> func,				//the function to find zero of
	T &x,									//an initial guess for the zero
	T dx,									//the step to use in numerical derivatives
	T const tol,							//tolerance of search, to be this accurate
	std::size_t const max_iter=0			//maximum number of iterations in search
);

//single-variable search, searching for func = target
//intended to work with real-valued or complex-valued numbers
template <typename T>
T newton_search(
	std::function<T(T)> func,				//the function to match to target
	T const target,							//the target to be matched to
	T &x,									//an initial guess for x
	T dx,									//the step to use in numerical derivatives
	T const tol,							//tolerance of search, to be this accurate
	std::size_t const max_iter=0			//maximum number of iterations in search
);

template <size_t np>
void newton_search(
	std::function<void(double f[np],double x[np])> func,	//f=vector function to zero, x=input array
	double (&x1)[np],							//an initial guess for x
	double (&dx)[np],							//the step to use in numerical derivatives
	double const tol=0.0,					//tolerance of search, to be this accurate
	std::size_t const max_iter=0,			//maximum number of iterations in search
	std::function<bool(double[np])> var_limit = no_limit<np> //a function limiting values of x1
);

template <size_t np>
void newton_search(
	std::function<void(double f[np],double x[np])> func,	//f=vector function, x=input array
	double (&target)[np], 						//the target to be matched to
	double (&x1)[np], 							//an initial guess for x
	double (&dx)[np],							//the step to use in numerical derivatives
	double const tol=0.0,					//tolerance of search, to be this accurate
	std::size_t const max_iter=0,			//maximum number of iterations in search
	std::function<bool(double[np])> var_limit = no_limit<np> //a function limiting values of x1
);

#include "rootfind.hxx"

#endif
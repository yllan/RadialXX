#include "radial.h"

#include <iostream>
using namespace std;

template<typename T>
T mq_dxx(int beta,  T x, T y, T xj, T yj, T c) 
{
  return x*x;
}


template<class RBF, typename T>
void fillM(T (*func)(int beta, T x, T y, T xj, T yj, T c), 
             RBF& rbf, const Polinomio<double>& pol, T c,
             const Vector<T>& x,const Vector<T>& y, int ini, int fin, Matrix<T>& A)
{
  int beta = rbf.get_beta();

}

int main(void)
{
 MQ<double>     rbf;
 Polinomio<double> pol;
 Vector<double> x,y;
 Matrix<double> A;
 double c;
 int  ni; 	
	   //fillMatrix(   "dx"     , rbf , c , pol , x , y, 0 ,  ni ,  A);;
	   
 fillM(rbf.eval, rbf, pol, c, x, y,  0, ni, A);
 	
}
//-------------------------------------------------------------------------
//
//  Copyright (C) 2009   Jose Antonio Munoz Gomez
//
//  This file is part of Radial++
//  http://sourceforge.net/projects/radial/
//
//  Radial++ is free software;  you can redistribute it and/or it under the
//  terms of the GNU Lesser General Public License as published by the Free 
//  Software Foundation; either version 3 of the License, or (at your option)
//  any later version.
//
//  Radial++ is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public 
//  License for more details.
//
//-------------------------------------------------------------------------

#include "radial.h"

#include "2d-comun.hpp"

//---------------------------------------------------------
template<typename RBF>
void test_rbf(RBF rbf)
{
// rbf:  RBF as kernel we can use : MQ, IMQ, GAU, TPS, POT
  Polinomio<double>   pol; 
  double              c;   
  Matrix<double>      A;   
  Vector<double>      x,y,f,b,lambda;
  Vector<double>      x_new,y_new,f_new;
  int                 n,m;
   
//define the number of nodes
   n = 20;   
   
//define the shape parameter for MQ,MQM, IMQ and GAU kernel
   c = 1.0;
   
//make the 2d data, see 2d-comun.hpp
   make_data(n,x,y,f);
   n = x.GetSize();
   
//configure the associate polynomial
// pol.make( data_dimension, degree_pol)  
   pol.make(2 , rbf.get_degree_pol());
   
//show rbf and pol info
   cout<<rbf;
   cout<<pol;      
   
//get the number of elements in the polynomial base  
   m = pol.get_M();
   
//make aditional data for:  Ax=b   
   b.Reallocate( n + m ); b = 0.0;
   
//fill the expanded vector b = [f 0]
   for( int i=0;  i<n;  i++)   
       b(i) = f(i);
       
//fill the Gramm matrix   
   fill_gram(rbf, pol, c, x, y, A);
   
//solve the linear system by LU factorization   
   lambda = solver::gauss(A,b);
   
//make a new grid to interpolate the data
   make_data(40, x_new,y_new, f_new);
   
//interpolate the data   
   f_new = interpolate(rbf,pol,c,lambda,x,y,x_new,y_new);
   
   cout<<endl;
   cout<<"Nnew = "<<x.GetSize()<<endl;
   cout<<"Nnew = "<<x_new.GetSize()<<endl;
   
//determine the maximum error
   double tmp,emax=0.0;
   for(int i=0; i<x_new.GetSize();i++)
   {
     tmp = fabs(f_new(i)-myf(x_new(i),y_new(i)));
     emax = tmp>emax? tmp : emax;   
   }
   
   printf("e_max =  %e\n",emax);    
   cout<<"----------------------------------------------------------------"<<endl;
}

//---------------------------------------------------------
int main(void)
{
   TPS<double>  tps;
   POT<double>  pot;   
   GAU<double>  gau;
   IMQ<double>  imq;
   MQ<double>    mq;
   MQM<double>  mqm;
   
   
   printf("Testing 2-d interpolation with different kernels\n");
   printf("we use the same data and c++ code and only pass a different kernel\n");
   printf("to perform the interpolation. This is powerful, because you can numerically\n");
   printf("experiment with differents radial basis functions without write more code.\n\n");
  
   test_rbf(tps);
   test_rbf(pot);
   test_rbf(imq);
   test_rbf(mq);
   test_rbf(gau);
   test_rbf(mqm);
      
  return 0;   
}   


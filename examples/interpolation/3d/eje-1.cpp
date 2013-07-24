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

//----------------------------------------------------------

double myf(double x, double y, double z)
{
  return sin(x)*exp(-y) + cos(y)*exp(-x) * sin(z)*exp(-x*y);
}
//----------------------------------------------------------
template<typename Vec>
void make_data(int N, Vec &x, Vec &y, Vec &z, Vec &f)
{
  double a,b,h;
  int    cont;   

   cont = 0;
   a    = 0.0;
   b    = 1.0;
   
   x.Reallocate(N*N*N);
   y.Reallocate(N*N*N);
   z.Reallocate(N*N*N);
   f.Reallocate(N*N*N);   
   
   h = (b-a)/double(N-1);   
   
   
   for(int i=0; i<N; i++ )
   {
     for(int j=0;j<N;j++)
     {
        for(int k=0; k<N; k++)
       {
          x(cont) = a + h*i;
          y(cont) = a + h*j;
         z(cont) = a + h*k;         
          cont++;
        }
     }
   }   

   
   for(int i=0; i<cont; i++ )
     f(i) = myf(x(i),y(i),z(i));   
   
}
//---------------------------------------------------------
int main(void)
{
  TPS<double>       rbf; // RBF as kernel we can use : MQ, IMQ, GAU, TPS, POT
  Polinomio<double> pol;   
  double            c;   
  Matrix<double>    A;   
  Vector<double>    x,y,z,f,b,lambda;
  Vector<double>    x_new,y_new,z_new,f_new;
  int               n,m;
 
//define the number of nodes
   n = 5;
      
//define the shape parameter for MQ, IMQ and GAU kernel
//in this example we use TPS, so, it is not necessary to give a value
   c = 1;
   
//make the 3d data
   make_data(n,x,y,z,f);
   n = x.GetSize();   
    
//set the exponent in  r^beta log(r)
   rbf.set_beta(4);
   
//configure the associate polynomial
// pol.make( data_dimension, degree_pol)  
   pol.make( 3 , rbf.get_degree_pol());
   
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
   fill_gram(rbf, pol, c, x, y, z, A);

//solve the linear system by LU factorization   
   lambda = solver::gauss(A,b);
   
//make a new grid to interpolate the data
   make_data(10,x_new,y_new, z_new, f_new);
   
//interpolate the data   
   f_new = interpolate(rbf,pol,c, lambda, x,y,z,x_new,y_new,z_new);
   
//store all the data in  ASCII format
   save_gnu_data("data_in.dat",x,y,z,f);   
   save_gnu_data("data_out.dat",x_new,y_new,z_new,f_new);

   
//determine the maximum error
   double tmp,emax=0.0;
   for(int i=0; i<x_new.GetSize();i++)
   {
      tmp = fabs(f_new(i)-myf(x_new(i),y_new(i),z_new(i)));
     emax = tmp>emax? tmp : emax;   
   }

   cout<<endl;
   cout<<"Testing 3-d interpolation"<<endl;
   cout<<" f(x,y,z) = sin(x)*exp(-y) + cos(y)*exp(-x) * sin(z)*exp(-x*y)"<<endl;
   cout<<"we build the interpolation with n-centers = "<<x.GetSize()<<endl;
   cout<<"determine the error in m centers = "<<x_new.GetSize()<<endl;
   cout<<"e_max = "<<emax<<endl<<endl;
   
 return 0;
}


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

//
//  Solve the Helmholtz equation in two-dimensions 
//      
//             Uxx + Uyy +k^2*U = f(x,y)     in (-1,1) x (-1,1)
//      
//  with boundary condition u=0 and    
//
//        f(x,y) = exp(-10*( (y-1)^2 + (x-1/2)^2  ))
// 
//  In particular, we employ k = 9 and rbf = r^7. This example
//  was taken from the book Meshfree Approximation Methods with
//  Matlab, page 413.
//  


#include "radial.h"

//---------------------------------------------------------
template <typename T>
T right_side(T xc, T yc)
{
   return exp(-10.0*( (yc-1.0)*(yc-1.0) + (xc-0.5)*(xc-0.5)   ));   
}
//---------------------------------------------------------
template <typename T>
T boundary_condition(T xc, T yc)
{
   return 0.0;   
}
//---------------------------------------------------------
int main(int argc, char **argv)
{
  POT<double>         rbf;
  Polinomio<double>   pol;
  Matrix<double>      A,B,C;
  Vector<double>      x,y,f;
  Vector<double>      lambda,b;
  Vector<double>      xnew,ynew,fnew; 
  double              c=0.01;
  int                 n,ni,m;
  double              k;
 
 
   k = 9.0; //requiered in  Uxx + Uyy + k^2*U = f(x,xy)

//make the data in the square [-1,1] x [-1,x]   
   make_data(-1,1,-1,1, 25, 25, x, y, ni, n);

//set the exponent in the radial basis function r^beta
   rbf.set_beta(7);
   
//configure the polynomial   
//pol.make( data_dimension, degree_pol)
   pol.make(2 , rbf.get_degree_pol());
   
//show the rbf and pol info
   cout<<rbf;
   cout<<pol;
   
//show the number of nodes   
   cout<<endl;
   cout<<"total nodes    N  = "<<n<<endl;
   cout<<"interior nodes ni = "<<ni<<endl;
   cout<<"boundary nodes nf = "<<n-ni<<endl;
   cout<<endl;
   

//get the number of elements in the polynomial base
   m = pol.get_M();
    
//creo la matriz A para la EDP
   A.Resize(n+m,n+m);  A = 0.0;
   B.Resize(n+m,n+m);  B = 0.0;
   C.Resize(n+m,n+m);  C = 0.0;
   
//make the derivatives matrices   
   fill_matrix(   "dxx"     , rbf , pol , c , x , y, 0 ,  ni ,  A);
   fill_matrix(   "dyy"     , rbf , pol , c , x , y, 0 ,  ni ,  B);
   fill_matrix(   "normal"  , rbf , pol , c , x , y, 0 ,  ni ,  C);
   
   A = A + B + (k*k)*C;  //  Uxx + Uyy + k^2*U 
      
//apply the moment conditions to the polynomial, they apply only in the
//boundary nodes       
   fill_matrix(   "normal"  , rbf , pol , c , x , y, ni,   n ,  A);
   fill_matrix( "pol_trans" , rbf , pol , c , x , y, n ,  n+m,  A);
     
   b.Resize(n+m); b = 0.0;
    
   for(int i=0;i<ni;i++)
      b(i) = right_side(x(i), y(i));
      
   for(int i=ni;i<n;i++)
      b(i) = boundary_condition(x(i),y(i));
      
//solve the linear system of equations 
   lambda =  gauss(A,b);
    
//make the new data grid    
   int ni2,n2;
   make_data(-1,1,-1,1, 50, 50, xnew, ynew, ni2, n2);
    
//interpolate in this new data grid (xnew,ynew)     
   fnew = interpolate(rbf,pol,c,lambda,x,y,xnew,ynew);
     
//now, interpolate in the node (0,0) 
   double fc,xc,yc;
   xc = 0.0;
   yc = 0.0;   
   fc = interpolate(rbf,pol,c,lambda,x,y,xc,yc);
    
//show the error and info    
   printf("\n");
   printf(" This equation was previously solved in the book of Trefethen(2000) and Fasshauer(pp413)\n");
   printf(" Meshfree approximation methods with matlab, book Fasshauer\n");
   printf(" the error reported in Fasshauer book was: u(0,0) = 0.01172256911\n");
   printf(" we determine the error as: u(0,0) = %.10f\n",fc);
   printf("\n");
         
//store the interpolated data  
   save_gnu_data("data",xnew,ynew,fnew);
     
  return 0;
}



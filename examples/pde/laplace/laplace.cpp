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
//  Solve the Laplace problem 
//      
//             Uxx + Uyy = 0     in (-1,1) x (-1,1)
//      
//  with boundary conditions
//                   
//                  |   sin^4(pi*x)      y=1 and  -1<x<0
//        u(x,y) =  | 0.2*sin(3*pi*y)    x=1
//                  |     0              otherwisze.
//  
// 
//  In particular, we employ the thin-plate spline r^4 log (r).
//  
//  This example was taken from the book Meshfree Approximation Methods with
//  Matlab, page 416.
//

  
#include "radial.h"

//---------------------------------------------------------
template <typename T>
T right_side(T xc, T yc)
{
   return 0.0;   
}
//---------------------------------------------------------
template <typename T>
T boundary_condition(T xc, T yc)
{
   if(xc==1)
      return 0.2*sin(3.0*M_PI*yc);
   else if( (yc==1) & ( (xc>-1) & (xc<0)))
      return pow(sin(M_PI*xc) ,4.0);      
   else
      return 0.0;   
}
//---------------------------------------------------------
int main(int argc, char **argv)
{
  TPS<double>         rbf;
  Polinomio<double>   pol;
  Matrix<double>      A,B;
  Vector<double>      x,y,f;
  Vector<double>      lambda,b;
  Vector<double>      fnew; 
  double              c=0.01;
  int                 n,ni,m;

//make the data in the square [-1,1] x [-1,1]   
   make_data(-1,1,-1,1, 31, 31, x, y, ni, n);

//set the exponent in TPS
   rbf.set_beta(4);
   
//configure the associate polynomial
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
    
//resize the matrices to store the partial derivatives
   A.Resize(n+m,n+m);  A = 0.0;
   B.Resize(n+m,n+m);  B = 0.0;
   
 
//Recall that this problem has the general form
//   (Uxx+Uyy)   (Pxx+Pyy)   =  f     interior nodes  0..ni 
//     U           P_b       =  g     boundary nodes  ni..n 
//   P^transpose   0         =  0     momentun conditions in P  
//
// P is the polynomial wit size m x n
// P_b is the polynomial working in the boundary nodes, size   nf x m
// Pxx+Pyy  has size  ni x m

//make the matriz derivatives    
   fill_matrix(   "dxx"     , rbf , pol , c , x , y, 0 ,  ni ,  A);
   fill_matrix(   "dyy"     , rbf , pol , c , x , y, 0 ,  ni ,  B);
      
   A = A + B;   //  A <-  Uxx + Uyy
      
//Add the submatriz for the boundary nodes:   U , P_b           boundary nodes  ni..n        
   fill_matrix(   "normal"  , rbf , pol , c , x , y, ni,   n ,  A);
   
//Add the submatriz P^transpose at the end:    P^transpose 
   fill_matrix( "pol_trans" , rbf , pol , c , x , y, n ,  n+m,  A);
     
//resize the vector to store the right_size of the PDE     
   b.Resize(n+m); b = 0.0;
    
//fill with   f  
   for(int i=0;i<ni;i++)
      b(i) = right_side(x(i), y(i));
      
//fill with the boundary condition       
   for(int i=ni;i<n;i++)
      b(i) = boundary_condition(x(i),y(i));
      
//solve the linear system of equations 
   lambda =  gauss(A,b);
    
//interpolate in the same original grid     
   fnew = interpolate(rbf,pol,c,lambda,x,y,x,y);
     
//now, interpolate in the node (0,0)     
   double fc,xc,yc;
   xc = 0.0;
   yc = 0.0;   
   fc = interpolate(rbf,pol,c,lambda,x,y,xc,yc);
    
//show the error     
   printf("\n");
   printf(" This equation was previously solved in the book of Trefethen(2000) and Fasshauer(pp416)\n");
   printf(" the value reported in Fasshauer book was: u(0,0) = 0.0495940466\n");
   printf(" we determine the numerical approximation: u(0,0) = %.10f\n\n",fc);
   
//store the interpolated data 
   save_gnu_data("data",x,y,fnew);
     
  return 0;
}



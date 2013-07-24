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
//  Solve the Poisson problem 
//      
//             Uxx + Uyy = f(x,y)     in (0,1) x (0,1)
//      
//  with boundary conditions
//                   
//        u(x,y) =  sin(2*pi*x)*cos(2*pi*y)
//   
//  the rhs is determined by applying the partial derivatives over u(x,y).
//  The exact solution of this problem is u(x,y) depicted above.      
//  
//  In particular, we employ the five radial basis functions given in this
//  library. !This is the power of Radial++, solve the same proble with 
//  differents kernels  without coding new functions. The degree of the
//  polynomial associate to each RBF depend of each selected RBF.
//

#include "radial.h"

//---------------------------------------------------------
template <typename T>
T right_side(T xc, T yc)
{
   return -2.0*(2.0*M_PI)*(2.0*M_PI)*sin(2.0*M_PI*xc)*cos(2.0*M_PI*yc);
}
//---------------------------------------------------------
template <typename T>
T boundary_condition(T xc, T yc)
{
   return sin(2*M_PI*xc)*cos(2*M_PI*yc);   
}
//---------------------------------------------------------
template<typename RBF>
void poisson_all(RBF rbf, int beta, double c)
{
  Polinomio<double>   pol;
  Matrix<double>      A,B;
  Vector<double>      x,y,f;
  Vector<double>      lambda,b;
  Vector<double>      xnew,ynew,fnew; 
  int                 n,ni,m;
  double              error;
  
   cout<<"-----------------------------------"<<endl;   
   
//make the data in the square [0,1] x [0,1]   
   make_data(0,1,0,1, 21, 21, x, y, ni, n);

//stablish the exponent in the radial basis function
   rbf.set_beta(beta);
   
//configure the associate polynomial
// pol.make( data_dimension, degree_pol)
   pol.make(2 , rbf.get_degree_pol());

   
//show the rbf and pol info
   cout<<rbf;
   cout<<pol;   
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
// P is the polynomial wit size n x m
// P_b is the polynomial working in the boundary nodes, size   nf x m
// Pxx+Pyy  has size  ni x m   
      
//make the matriz derivatives   
   fill_matrix(   "dxx"     , rbf , pol , c , x , y, 0 ,  ni ,  A);
   fill_matrix(   "dyy"     , rbf , pol , c , x , y, 0 ,  ni ,  B);
 
   A = A + B; //  A <-  Uxx + Uyy
 
//Add the submatriz for the boundary nodes:   U , P_b           boundary nodes  ni..n   
   fill_matrix(   "normal"  , rbf , pol , c , x , y, ni,   n ,  A);
   
//Add the submatriz P^transpose at the end:    P^transpose 
   fill_matrix( "pol_trans" , rbf , pol , c , x , y, n ,  n+m,  A);
    
//resize the vector to store the right_size of the PDE   
   b.Resize(n+m); b = 0.0;
    
//fill with  rhs of the PDE     
   for(int i=0;i<ni;i++)
      b(i) = right_side(x(i), y(i));
      
//fill with the boundary condition        
   for(int i=ni;i<n;i++)
      b(i) = boundary_condition(x(i),y(i));
      
    
//solve the linear system of equations  
   lambda =  solver::gauss(A,b);
  
//make the new data grid        
   int ni2,n2;
   make_data(0,1,0,1, 41, 41, xnew, ynew, ni2, n2);
   
//interpolate on this new data grid (xnew,ynew)     
   fnew = interpolate(rbf,pol,c,lambda,x,y,xnew,ynew);
     
//determine the maximum error     
   error = 0.0;  
   for(int i=0;i<ni2;i++)
   {
     error = max(error, fabs(fnew(i) - sin(2*M_PI*xnew(i))*cos(2*M_PI*ynew(i)))  );
   }
     
//show the error     
   printf("\n");
   printf("The maximum error: e_max = %e\n",error);
   
}
//---------------------------------------------------------
int main(int argc, char **argv)
{
   TPS<double> tps;
   POT<double> pot;   
   GAU<double> gau;
   IMQ<double> imq;
   MQ<double>  mq;
  
   printf("\n");
   printf("In this example, we demostrate the power of Radial++, in particular\n");
   printf("we solve the same Poisson problem  with differents kernels.\n");
   printf("\n");
   
   poisson_all(tps,4,0);
   poisson_all(pot,5,0);
   poisson_all(imq,1,0.851);
   poisson_all(mq,1,0.746);
   poisson_all(gau,1,3.777);
   
  return 0;   
}

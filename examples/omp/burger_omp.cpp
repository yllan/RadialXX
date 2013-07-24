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
//  Solve the viscous Burger's equation in two-dimensions 
//      
//             Ut + grad f(u) = eps . (Uxx + Uyy) 
//      
// with zero at Dirichlet data and initia condition a
// Gaussian function truncated. We employ the Euleer method for time
// discretization, the numerical method is the unsymmetric collocation
// method with  r^4 log( r ) as kernel function.
// 
// To compile: 
//             c++ -O3  -I../../../src burger_omp.cpp -fopenmp -o bur_omp
// To run:
//           ./bur_omp

#ifdef WITH_LAPACK  
  #define SELDON_WITH_CBLAS
  #define SELDON_WITH_LAPACK
#endif  

#define  SELDON_WITH_OMP
#include <omp.h>
#include "radial.h"


double difusivo = 0.005;  //diffusive parameter = eps

//---------------------------------------------------------
double u_t0(double x,double y)
{
  double tmp1,tmp2,z;
  double xc = 0.35;
  double yc = 0.35;
  double R  = 0.2;
   
   tmp1 = (x-xc)*(x-xc) + (y-yc)*(y-yc);
   tmp2 = sqrt(tmp1);
   
   if(tmp2<=R){
      z = exp( tmp1/(tmp1-R*R));
   }
   else{
      z=0.0;
   }
   
 return z;
}
//---------------------------------------------------------
template<typename Vec>
void initial_condition(Vec &x, Vec &y, Vec &u, int ni)
{
   for(int i=0; i<ni; i++)
     u(i) = u_t0(x(i),y(i));

}
//---------------------------------------------------------
template<typename T>
T boundary(T t,const T x,const T y)
{
   return  0.0;
} 
//---------------------------------------------------------
int main(int argc, char **argv)
{
  TPS<double>         rbf;
  Polinomio<double>   pol;
  LU<double>          lu;
  Matrix<double>      A,B,Uxy,Uxxyy,H,U;
  Vector<double>      x,y,f;
  Vector<double>      u,lambda,b;
  double              c=0.02;
  int                 n,ni,m;
  char                filename[255];   
  timer               time;
       
       
   printf("Testing the viscous Burger's equation\n");
   printf("\n");
       
//make the grid points in [0,1] x [0,1]
//the data is sorted as: *_d + *_b d=0,1,..,ni-1 b=ni,ni+1,..,n-1
//where ni denote the interior nodes and b the boundary nodes.
   make_data(0,1,0,1, 30, 30, x, y, ni, n);

//set the exponent in TPS 
   rbf.set_beta(4);
   
//configure the associate polynomial
// pol.make( data_dimension, degree_pol)
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
   
//create the initial condition at t=0  
   b.Reallocate(n+m);  b = 0.0;
   initial_condition(x,y,b,ni);
   
//create the required matrices to solve the PDE problem
   A.Resize(n+m,n+m);      A = 0.0;
   B.Resize(n+m,n+m);      B = 0.0;
   Uxy.Resize(n+m,n+m);    Uxy = 0.0;
   Uxxyy.Resize(n+m,n+m);  Uxxyy = 0.0;    
   H.Resize(n+m,n+m);      H = 0.0;    
   U.Resize(n+m,n+m);      U = 0.0;   

      
//build the spatial derivative matrices according to:   Ut =  -U*(Ux + Uy) + eps * (Uxx + Uyy)
   fill_matrix(   "dx"     , rbf , pol , c , x , y, 0 ,  ni ,  A);;
   fill_matrix(   "dy"     , rbf , pol , c , x , y, 0 ,  ni ,  B);

    Uxy = A + B;  //  Ux + Uy
    
   fill_matrix(   "dxx"     , rbf , pol , c , x , y, 0 ,  ni ,  A);
   fill_matrix(   "dyy"     , rbf , pol , c , x , y, 0 ,  ni ,  B);
   
    Uxxyy =  A + B;  //  Uxx + Uyy    
    
   A.Resize(0,0);
   B.Resize(0,0);
     
   fill_matrix(   "normal"  , rbf , pol , c , x , y, 0,   ni ,   U);

//   U_n+1 = U_n + ---
//   H lambda_n+1 =  U_n + --- 
   fill_matrix(   "normal"  , rbf , pol , c , x , y, 0,   n ,   H);
   fill_matrix( "pol_trans" , rbf , pol , c , x , y, n ,  n+m,  H);
    

//LU factorization over the matriz H, overwrite en H
    cout<<"running the LU matrix factorization .... "<<endl;
    lu.factorize(H,1);

//H contain the LU factorization: LU <-H
//  solve  Lz=b, Ux=z   x=lambda
   lambda = lu.solve(H,b);
   
//get the initial approximation at t=0   
   u = U*lambda;
     
//store the initial data solution at t=0  
   save_gnu_data("data/data0",x,y,b);
 
//now, the iterative temporal scheme Euler 
   cout<<"running the simulation .... "<<endl;
    
   double  stept = 0.0001;
   double  t     = 0.0;
   double  tmax  = 1.0;  
   int     cont  = 0;

//validate the fixed stept time
   check_stept(tmax, stept);
     
   time.tic();
   
   double te = omp_get_wtime();
     
   while(t<tmax)
   {        
     //correct the stept time
         if( correct_stept(t, tmax, stept) == false )
           break;
            
    //Euler Method for  Ut =  -U*(Ux + Uy) + eps * (Uxx + Uyy)
          u = U*lambda;
        
          b = u - stept*( pdot(u,Uxy*lambda)  )   +  (stept*difusivo)*( Uxxyy * lambda );
      
        // vector <-- pdot(vector,vector)
        //            like '.*' in matlab
        
        
    //apply boundary conditions
        for(int i=ni; i<n; i++) 
            b(i) = boundary(t+stept,x(i),y(i));
          
    //set zero in the polynomial part
         for(int i=n; i<(n+m); i++) 
             b(i) = 0.0;

    //solve the linear system of equations, H contain the LU=A
     lambda = lu.solve(H,b);
        
    //advance the time step  
         t=t+stept;

        if( (++cont%1000)==0){
 
           u = U*lambda;      
           sprintf(filename,"data/data%d",cont);
           save_gnu_data(filename,x,y,u);      

           printf("t = %.4f",t);
           printf("  max(u) = %.4f",max(u));
           printf("  stept  = %.4f   %s\n",stept,filename);

         }
  }

   te  = omp_get_wtime() - te;
   cout<<endl;
   cout<<"Elapsed time (loop) = "<<time.toc()<<"  sec"<<endl;
   cout<<"Elapsed time (loop) = "<<te<<"  sec"<<endl;
   cout<<endl;
      
//store the final data solution at t=tmax   
   save_gnu_data("data/datan",x,y,u);

   return 0;
}

   



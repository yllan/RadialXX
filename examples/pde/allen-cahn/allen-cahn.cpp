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
//  Solve the one dimensional Allen-Cahn equation 
//      
//             Ut = mu * Uxx + U - U^3   in (-1,1)
//      
//  with parameter mu>0, and initial condition
//  
//             u(x,0) = 0.53x + 0.47sin(-3/2 pi*x)
//         
//  together with boundary conditions
//  
//             u(-1,t) = -1   and  u(1,t) =  1 + sin^2(t/5)         
//             
//  This example was taken from the book Meshfree Approximation Methods with
//  Matlab, page 409.           
//


#ifdef WITH_LAPACK  
  #define SELDON_WITH_CBLAS
  #define SELDON_WITH_LAPACK
#endif  
#include "radial.h"

double mu = 0.01;

//---------------------------------------------------------
double u_t0(double x)
{
   return 0.53*x + 0.47*sin(-1.5*M_PI*x);
}
//---------------------------------------------------------
template<typename Vec>
void initial_condition(Vec &x,  Vec &u, int ni)
{
   for(int i=0;i<ni;i++)
     u(i) = u_t0( x(i) );
}
//---------------------------------------------------------
template<typename T>
T boundary(T t,const T x)
{
   if(x==-1)
     return -1;
   else
     return 1.0 + sin(t*0.2) * sin(t*0.2);  
}

//---------------------------------------------------------
int main(int argc, char **argv)
{
  TPS<double>         rbf;
  Polinomio<double>   pol;
  LU<double>          lu;
  Matrix<double>      U,Uxx,H;
  Vector<double>      x,f;
  Vector<double>      u,u3,lambda,b;
  double              c=0.01;
  int                 n,ni,m;
  char                filename[255];
  

   cout<<"Testing the Allen-Cahn equation"<<endl;
   cout<<endl;    
      
//make the grid points
   make_data(-1,1, 51,  x, ni, n);
    
//configure the polynomial  
   rbf.set_beta(4);
 
   pol.make(1 , rbf.get_degree_pol());

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
   b.Resize(n+m); b = 0.0;
   initial_condition(x,b,ni);

//create the required matrices to solve the PDE problem
   U.Resize(n+m,n+m);    U   = 0.0;
   Uxx.Resize(n+m,n+m);  Uxx = 0.0; 
   H.Resize(n+m,n+m);    H   = 0.0;   

      
//fill the matrices with the derivatives  
   fill_matrix(   "dxx"    , rbf , pol , c , x , 0 ,  ni ,  Uxx);
   fill_matrix(   "normal" , rbf , pol , c , x , 0,   ni ,  U);

//fill the matrix H required in: 
//   U_n+1 = U_n + ---
//   H lambda_n+1 =  U_n + --- 
   fill_matrix(   "normal"  , rbf , pol , c , x ,  0,   n ,   H);
   fill_matrix( "pol_trans" , rbf , pol , c , x ,  n ,  n+m,  H);


//LU factorization over the matriz H, overwrite en H
    lu.factorize(H,1);
   
//H contain the LU factorization: LU <-H
//  solve  Lz=b, Ux=z   x=lambda
   lambda = lu.solve(H,b);
   
//get the initial approximation at t=0   
   u = U*lambda;
      
//store the final data solution at t=0      
   save_gnu_data("data/data0",x,b);
 
//Entro al esquema iterativo temporal    
   cout<<"running the simulation .... "<<endl;

   double  stept = 0.001;
   double  t     = 0;
   double  tmax  = 50;
   int     cont  = 0;

//validate the fixed stept time  
   check_stept(tmax, stept);
 
  while( t<tmax )
  {        
        //correct the stept time
         if( correct_stept(t, tmax, stept) == false )
            break;
          
        //Euler Method
            u  = U*lambda;
            b  = u  +  (stept*mu)*(Uxx*lambda) + stept*u - stept*(u^3);
            
            //it is convenient to use pow?(u) ?=2,3,4,5,6,7,8 similar .^ in matlab

             
        //apply boundary conditions
         for(int i=ni; i<n; i++) 
           b(i) = boundary( t+stept ,x(i));
          
       //set zero in the polynomial part
         for(int i=n; i<(n+m); i++) 
             b(i) = 0.0;

       //solve the linear system of equations, H contain the LU=A
          lambda = lu.solve(H,b);
        
      //advance the time step  
         t = t + stept;
         
      //store the current solution u  
        if( (++cont%10000)==0){
        
           sprintf(filename,"data/data%d",cont);
           save_gnu_data(filename,x,u);  
           
           printf("t = %.4f",t);
           printf("  max(u) = %.4f",max(u));
           printf("  stept  = %.4f   %s\n",stept,filename);
        
         }
         
 } //while(t<tmax)
  
   printf("t = %.4f",t);
   printf("  max(u) = %.4f",max(u));
   printf("  stept  = %.4e\n",stept);
   
//store the final data solution at t=tmax        
   save_gnu_data("data/datan",x,u);
   
  return 0;
}



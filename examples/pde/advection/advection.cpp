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
//  Solve the advection equation in two-dimensions 
//      
//             Ut + Ux + Uy = 0
//      
// with zero at Dirichlet data and initia condition a
// Gaussian function at the point (0.3,0.3).
// 
// We employ the MQ kernel with c = 0.1.
//


#ifdef WITH_LAPACK  
  #define SELDON_WITH_CBLAS
  #define SELDON_WITH_LAPACK
#endif  


#include "radial.h"


//---------------------------------------------------------
template<typename T>
T solution(T t, T x, T y)
{
  return  exp(-1.0*(( (x-t-0.3)*(x-t-0.3)+(y-t-0.3)*(y-t-0.3) ) / 0.005));

}//---------------------------------------------------------
template<typename Vec>
void initial_condition(Vec &x, Vec &y, int ni, Vec &u)
{
   for(int i=0;i<ni;i++)
     u(i) = exp(-1.0*(( (x(i)-0.3)*(x(i)-0.3)+(y(i)-0.3)*(y(i)-0.3) ) / 0.005));

}
//---------------------------------------------------------
template<typename T>
T boundary(T t,const T x,const T y)
{
   return  0.0;
}
//---------------------------------------------------------
template <typename RBF>
void advection(RBF rbf, int beta, int N)
{
  Polinomio<double>   pol;
  LU<double>          lu;
  Matrix<double>      U,Ux,Uy,H;
  Vector<double>      x,y,f;
  Vector<double>      u,lambda,b;
  double              c=0.01;
  int                 n,ni,m;
  char                filename[255];
  timer               time;
  
  
   cout<<"Testing the advection equation"<<endl;
   cout<<endl;    
   
//shape parameters
   c = 0.12;
   
//make the grid points in [0,1] x [0,1]
   make_data(0,1,0,1, N, N, x, y, ni, n);
    
//set the exponent in TPS  
   if(beta>0){
      rbf.set_beta(beta);
   }
   
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
    
//create the initial condition at t=0    
   b.Reallocate(n+m); b = 0.0;
   
   initial_condition(x,y,ni,b);
   

//create the required matrices to solve the PDE problem
   U.Reallocate(n+m,n+m);   U  = 0.0;
   Ux.Reallocate(n+m,n+m);  Ux = 0.0;
   Uy.Reallocate(n+m,n+m);  Uy = 0.0;    
   H.Reallocate(n+m,n+m);   H  = 0.0;   

      
//fill the matrices with the derivatives  
   fill_matrix(   "dx"     , rbf , pol , c , x , y, 0 ,  ni ,  Ux);
   fill_matrix(   "dy"     , rbf , pol , c , x , y, 0 ,  ni ,  Uy);
     
   fill_matrix(   "normal" , rbf , pol , c , x , y, 0,   ni ,  U);

//fill the matrix H required in: H lambda_n+1 =  H lambda_n - stept(Ux+Uy) lambda_n 
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
      
//store the final data solution at t=0
   save_gnu_data("data/data0",x,y,b);
 
//now, the Euler scheme    
   cout<<"running the simulation .... "<<endl;

   double  stept = 0.0001;
   double  t     = 0;
   double  tmax  = 0.4;
   int     cont  = 0;

//validate the fixed stept time
   check_stept(tmax, stept);

     
   Ux =  stept*(Ux + Uy);


   time.tic();
   
  while( t<tmax )
  {        
        //correct the stept time
         if( correct_stept(t, tmax, stept) == false )
            break;
          
        //Euler Method
           u = u -   Ux * lambda ;  //MltAdd(-1.0, Ux, lambda, 1.0, u);            
            
       //set zero in the boundary and the polynomial part
           for(int i=ni; i<(n+m); i++) u(i) = 0.0;

       //solve the linear system of equations, H contain the LU=A
          lambda = lu.solve(H,u);

        
       //advance the time step  
          t = t + stept;
         
      //store the current solution u  
        if( (++cont%100)==0){
        
           u = U*lambda;
           sprintf(filename,"data/data%d",cont);
           save_gnu_data(filename,x,y,u);  
           
           printf("t = %.4f",t);
           printf("  max(u) = %.4f",max(u));
           printf("  stept  = %.4f   %s\n",stept,filename);
        
         }

         
 } //while(t<tmax)
 
   cout<<endl;
   cout<<"Elapsed time (loop) = "<<time.toc()<<" sec"<<endl;
   cout<<endl;

//store the final data solution at t=tmax    
   u = U*lambda;
   save_gnu_data("data/datan",x,y,u);
   
//determine the maximum error
   double e=0.0;
   double tmp;
   for(int i=0;i<n;i++)
   {
        tmp  = u(i) - solution(t,x(i),y(i));
        e  = max(e, fabs(tmp) ) ;
   }
     
//show the error     
   cout<<endl;
   cout<<"The maximum error: e_max = "<<e<<endl;  

}
//---------------------------------------------------------
int main(int argc, char **argv)
{
   MQ<double>     rbf_mq;
   TPS<double>     rbf_tps;
   POT<double>    rbf_pot;
   int   op,n;
  
   if(argc<3){
      printf("./advection <op> n\n");
      printf("       op = 1  use MQ  multiquadric with c=0.12\n");
      printf("       op = 2  use TPS thin-plate-spline r^4 log r\n");    
      printf("       op = 3  use POT  r^5\n");
      printf("        n = number of nodes, n x n\n");
      return 1;
   }
  
   op = atoi(argv[1]);
   n  = atoi(argv[2]);
  
   switch(op)
   {
     case 1: advection(rbf_mq,-1,n);
            break;
     case 2: advection(rbf_tps,4,n);
            break;
     case 3: advection(rbf_pot,5,n);
            break;
   }
  
  return 0;
}




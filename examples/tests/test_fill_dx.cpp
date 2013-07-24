/***************************************************************************
 *
 * Copyright (C) 2009   Jose Antonio Munoz Gomez
 *        
 * This file is part of Radial++
 * http://sourceforge.net/projects/radial/
 *
 * Radial++ is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 3
 * of the License, or (at your option) any later version.
 *
 * Radial++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * For more information, please see the  Home Page:
 *    http://www.dci.dgsca.unam.mx/pderbf
 *
 ***************************************************************************/

#include "radial.h"

#include "2d-comun.hpp"



//----------------------------------------------------------
// Data en 2-D
//----------------------------------------------------------
template<typename RBF, typename Pol, typename Mat, typename Vec>
void fillMatrix2( RBF rbf,  double c, Pol &pol, Vec &x, Vec& y, int ini, int fin, Mat &A)
//
//Build the Gramm matrix
//
//  ! A    P !  lambda   = f
//  ! P^t  0 !    mu      = 0
//
{
 int    m;     // degre of the polynomial, recal that is  m-1
 int    M;     // number of elements in the base of the polynomial
 int    N;     // number of elements in the vector x
 Mat    P;     // Matrix to store the evaluations of the Polynomial 


//get the degree of the polynomial required in the RBF
//recall that in practice, we employ the degree  m-1   
   m = rbf.get_degree_pol();      
   

//get the number of elements in the polynomial base
   M = pol.get_M();

//get the size of the vector
    N = x.GetSize();
    



         
  //fill the matrix with the evaluations of RBFs
    double* y_ = y.GetData();
    double* x_ = x.GetData();
  //  double* z_ = A.GetData(); 
    double xc,yc,zc;
    
   for(int i=ini; i<fin; i++)
   {
       xc = x_[i];
       yc = y_[i];
     for(int j=0; j<N; j++)
     {
                                                              //options -finline and  -finline -O3
                                       //original             //MFLOPS = 17  65                       
        // A(i,j) = rbf.eval(x(i),y(i),x(j),y(j),c);         //MFLOPS = 28  108   TP1
        // A(i,j) = rbf.eval(xc,yc,x(j),y(j),c);             //MFLOPS = 30  100   TP2
        // A(i,j) = rbf.eval(x_[i], y_[i], x_[j], y_[j],c);  //MFLOPS = 34  116   TP3
        // A(i,j) = rbf.eval(xc,yc,x_[j],y_[j],c);           //MFLOPS = 35  118   TP4
           
       //A(i,j) = sqrt( (x_[i]-x_[j])*(x_[i]-x_[j]) + (y_[i]-y_[j])*(y_[i]-y_[j]) + c*c);  //MFLOPS = 35  120
        //A(i,j) = sqrt( (xc-x_[j])*(xc-x_[j]) + (yc-y_[j])*(yc-y_[j]) + c*c);             //MFLOPS = 43  121
         //z_[i*N + j] = (xc-x_[j]) /sqrt( (xc-x_[j])*(xc-x_[j]) + (yc-y_[j])*(yc-y_[j]) + c*c);        //MFLOPS = 72  121
        A(i,j) = (xc-x_[j]) /sqrt( (xc-x_[j])*(xc-x_[j]) + (yc-y_[j])*(yc-y_[j]) + c*c); 
     }
   }


  //if m==0 indicate that does not require a polynomial
   /*
   if(m>0){
      
      P = pol.build_tnt(x.GetData(), N,2);   //build the matrix P of the polynomial  evaluations 
                                              // on the points x
      
      //insert the polinomial into the Gramm matrix
      for(int i=0; i<N; i++)
      {
      for(int j=N; j<(N+M); j++)
      {   
         A(i,j) = P(i,j-N);
         A(j,i) = A(i,j);
      }   
   }
 }
 */
}

template<typename RBF,   typename Vec>
void test( RBF rbf,  double c,  Vec &x)
{
}
template<typename RBF, typename T>
void test(RBF rbf, Vector<T> val)
{
 T *data;   
}
//---------------------------------------------------------
int main(int argc, char **argv)
{
  MQ<double>          rbf; // RBF as kernel we can use : MQ, IMQ, GAU, TPS, POT
  Polinomio<double>   pol; 
  double              c;   
  Matrix<double>      A;   
  Vector<double>      x,y,f,b,lambda;
  Vector<double>      x_new,y_new,f_new;
  int                 ni,n,m;
  timer               time;
  double             wtime;
  double             mflops,ops;
  Vector<double>    myc;

   
   test( rbf, myc);
   
//define the shape parameter for MQ kernel
   c = 1.0;
   

//make the data in the square [0,1] x [0,x]   
   make_data(0,1,0,1, 32, 32, x, y, ni, n);   
   
//number of operations to fill the Gramm-Matrix
   ops = 8*n*n;

//configure the associate polynomial
// pol.make( data_dimension, degree_pol)
   rbf.set_degree_pol(5);
   pol.make(1 , rbf.get_degree_pol() );
   
//show rbf and pol info
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
   cout<<"mm = "<<m<<endl;
   
   
   A.Reallocate(n+m,n+m);
   myc.Reallocate(n);
   
   myc = c;
   

    
  for(int iter=0;  iter<10;  iter++)
  {  
     
     //fill the Gramm matrix   
  time.tic();
  
  for(int i=0;  i<10;  i++)
     fillMatrix(   "dx"     , rbf , myc , pol , x , y, 0 ,  ni ,  A);

  wtime  = time.toc()/10.0;
  mflops = ops/(wtime*1.0E+06); 
  

  printf("Elapsed time = %5f    ",wtime);
  cout<<"MFLOPS1       = "<< mflops <<endl;   
     
    
    

  

  }
  

   
 return 0;
}



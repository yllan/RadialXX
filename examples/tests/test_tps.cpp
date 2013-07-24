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

//---------------------------------------------------------
template<typename T>
void pol_1d_x(Polinomio<T> pol)
{
	pol.deriva("x",1);
}
template <typename T, class RBF>
Vector<T> interpola(T (*func)(int,T,T),
                   void (func_pol)(Polinomio<T>),
                   Vector<T>& x_new)
{
   Vector<T>   px;   
   Vector<T> out;
   int   m,dim_px,N,Nnew;   
   T     xnew[2]; 
   T     s;
   int   i,j;
   int   beta,dim;
   
   dim    = 1;
   
   
   Nnew    = x_new.GetSize();
   
   out.Reallocate(Nnew);
   
   beta = rbf.get_beta();
   
   cout<<"nnew = "<<Nnew<<endl;
   return out;	
      /*
   
   for(i=0; i<Nnew;i++)
   {
   s    = 0.0;
   for(int j=0; j<N; j++)
   {
    s = s +  lambda(j) * func( beta, x_new(i), x(j));
   }
    out(i) = s;
   }
   
      m = pol.get_M();
  

	
	  func_pol(pol);
	
  if(m>0) //require a polynomial in the interpolation
   {
      dim_px = pol.get_M();
      
      cout<<"dim_px = "<<dim_px<<endl;
          
      
      px.Reallocate(dim_px);
      
      for(int i=0; i<Nnew; i++)
      {
         pol.eval(&x_new(i*dim),dim,px.GetData(),dim_px);
         s = 0.0;
         for(int k=0; k<dim_px; k++)
         {
            s = s + lambda(k+N)*px(k);   
         }      
         out(i) = out(i) + s;   
      }         
   }
*/
 
}
//---------------------------------------------------------
template <typename T, class RBF>
void fill_gram( T (*func)(int,T, T),
                 RBF rbf,
                 Polinomio<T>& pol,
                 Vector<T>&x,
                 Matrix<T>& A )
{
   int        n,m,M;
   int        i,j;
   int        beta;
   Matrix<T>  P;     // Matrix to store the evaluations of the Polynomial 
    int    dim=1;
   
//get the degree of the polynomial required in the RBF
//recall that in practice, we employ the degree  m-1   
   m = rbf.get_degree_pol();      
   

//get the number of elements in the polynomial base
   M = pol.get_M();
   
   n = x.GetSize();
   
   beta = rbf.get_beta();
   
   for(i = 0; i < n; i++)
   {
      for(j = 0; j < n; j++)
      {
        A(i,j)	= func(beta,x(i),x(j));
      }
      	
   }

//if m==0 indicate that does not require a polynomial
   if(m>0)
   {

      P = pol.build_tnt(x.GetData(), n, dim);   //build the matrix P of the polynomial  evaluations 
                                              // on the points x
      
      //insert the polinomial into the Gramm matrix
      for(int i=0; i<n; i++)
      {
         for(int j=n; j<(n+M); j++)
        {   
              A(i,j) = P(i,j-n);
              A(j,i) = A(i,j);
        }   
       } // if(m>0)
   }   
   
}     
//---------------------------------------------------------
int main(int argc, char **argv)
{
   timer times;
   POT<double>       rbf;
   Polinomio<double> pol;
   Vector<double>    x,x_new,b,f,f_new;
   Vector<double>    lambda;
   Matrix<double>    A;
   double            e_max;
   int               m;
   int               n=100;

   x = linspace(0,1.0,n);
   
   f.Reallocate(n);
   
   f = sin(x);
     
   pol.make(1, rbf.get_degree_pol());
   
//get the number of elements in the polynomial base  
   m = pol.get_M();
   
  
   A.Reallocate(n+m,n+m); A = 0.0; 
   
//make aditional data for:  Ax=b   
   b.Reallocate( n + m ); b = 0.0;
   
//fill the expanded vector b = [f 0]
   for( int i=0;  i<n;  i++)   
       b(i) = f(i);
       
   fill_gram(pot_1d, rbf, pol, x, A);
   
   
//solve the linear system by LU factorization
   lambda = solver::gauss(A, b); 
   
   
//make a new grid to interpolate the data
   x_new = linspace( 0.0 , 1.0 , 2*n );
   
//interpolate the data        
     f_new = interpola(rbf,0.0,pol,lambda,x,x_new); 
   f_new = interpola(pot_1d_x, pol_1d_x,  x_new);
   
   return 0;
   
//determine the maximum error   
   e_max = 0.0;
   for(int i=0;i<2*n; i++)
   {   
     e_max = max(e_max, fabs( sin(x_new(i)) - f_new(i) ) );
   }
   
   cout<<"Testing 1-d interpolation"<<endl;
   cout<<" f(x) = sin(x)"<<endl;
   cout<<"we build the interpolation with n-centers = "<<n<<endl;
   cout<<"determine the error in 2*n centers = "<<2*n<<endl;
   cout<<"e_max = "<<e_max<<endl;

  return 0;
}
//---------------------------------------------------------






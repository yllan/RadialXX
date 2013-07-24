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


#ifndef _GRAMM_HPP_
#define _GRAMM_HPP_

namespace rbf{

template<typename T>
void FillZero(Matrix<T>& A, int ini, int fin)
{
   int i,j;
   
   for( i = ini; i <(ini+fin); i++ )
   {
      for( j=ini; j<(ini+fin); j++)      
      {   
        A(i,j)= 0.0;
      }
   }
}
//----------------------------------------------------------
// Data en N-D
//----------------------------------------------------------
template<typename RBF, typename Pol, typename Mat, typename Vec>
void fill_gram( RBF rbf,  Pol &pol, Vec &c, Vec &x, int dim, Mat &A)
//
//Build the Gramm matrix
//
//  ! A    P !  lambda   = f
//  ! P^t  0 !    mu     = 0
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
    
   if(dim>1){
     if((N%dim)>0){
       fprintf(stderr,"!!ERROR:  N = %d must be divisible by dim = %d.\n",N,dim);
       fprintf(stderr,"in '%s' near to line %d\n\n",__FILE__,__LINE__);
       exit(1);
     }
   }
    
    N = N/dim;
    
    assert( c.GetSize()==(x.GetSize()/dim) );


//resize the Gramm matrix
   A.Reallocate(N + M, N + M);

//fill with zeros   
   FillZero( A, N, M );

//fill the matrix with the evaluations of RBFs
   for(int i=0; i<N; i++)
   {
     for(int j=0; j<N; j++)
     {
         A(i,j) = rbf.eval(&x(i*dim),&x(j*dim),dim,c(j));
     }
   }


  //if m==0 indicate that does not require a polynomial
   
   if(m>0){
      
      P = pol.build_tnt(x.GetData(), N, dim);   //build the matrix P of the polynomial  evaluations 
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
   
}

//----------------------------------------------------------
// Data in N-D
//----------------------------------------------------------
template<typename RBF, typename T, typename Pol, typename Mat, typename Vec>
void fill_gram( RBF rbf,  Pol &pol, T c, Vec &x, int dim,Mat &A)
//
//Build the Gramm matrix
//
//  ! A    P !  lambda   = f
//  ! P^t  0 !    mu     = 0
//
{
   Vec  myc;

   myc.Reallocate(x.GetSize()/dim);

   myc = c;
  
   fill_gram(rbf,pol,myc,x, dim, A);
}

//----------------------------------------------------------
// Data 1D
//----------------------------------------------------------
template<typename RBF, typename T, typename Pol, typename Mat, typename Vec>
void fill_gram(RBF rbf, Pol &pol, T c, Vec &x, Mat &A)
{
//
//Build the Gramm matrix
//
//  ! A    P !  lambda   = f
//  ! P^t  0 !    mu     = 0
//
// old release  fillGramm(rbf,c,pol,x, 1, A);


 int    m;     // degre of the polynomial, recal that is  m-1
 int    M;     // number of elements in the base of the polynomial
 int    N;     // number of elements in the vector x
 Mat    P;     // Matrix to store the evaluations of the Polynomial 
 int    dim=1;

//get the degree of the polynomial required in the RBF
//recall that in practice, we employ the degree  m-1   
   m = rbf.get_degree_pol();      
   

//get the number of elements in the polynomial base
   M = pol.get_M();
   
//get the size of the vector
    N = x.GetSize();
    
//resize the Gramm matrix
   A.Reallocate(N + M, N + M);
   
//fill with zeros   
   FillZero( A, N, M );

         
//fill the matrix with the evaluations of RBFs
    T   xc;
    T*  x_ = x.GetData();
    
   for(int i=0; i<N; i++)
   {
      xc = x_[i];
      
     for(int j=0; j<N; j++)
     {
         A(i,j) = rbf.eval(xc,x_[j],c);
     }
   }
   
   
//if m==0 indicate that does not require a polynomial
   if(m>0)
   {

      P = pol.build_tnt(x.GetData(), N, dim);   //build the matrix P of the polynomial  evaluations 
                                              // on the points x
      
      //insert the polinomial into the Gramm matrix
      for(int i=0; i<N; i++)
      {
         for(int j=N; j<(N+M); j++)
        {   
              A(i,j) = P(i,j-N);
              A(j,i) = A(i,j);
        }   
       } // if(m>0)
   }   
   
}
//----------------------------------------------------------
// Data 1D
//----------------------------------------------------------
template<typename RBF, typename Mat, typename Pol, typename Vec>
void fill_gram(RBF rbf, Pol &pol, Vec &c, Vec &x, Mat &A)
{
//
//Build the Gramm matrix
//
//  ! A    P !  lambda   = f
//  ! P^t  0 !    mu     = 0
//

   fill_gram(rbf,pol,c,x, 1, A);

}
//----------------------------------------------------------
// Data 2D
//----------------------------------------------------------
template<typename RBF, typename T, typename Pol,typename Mat, typename Vec>
void fill_gram(RBF rbf, Pol &pol, T c,  Vec &x, Vec &y,  Mat &A)
{
//
//Build the Gramm matrix
//
//  ! A    P !  lambda   = f
//  ! P^t  0 !    mu     = 0
//

 int    m;     // degre of the polynomial, recal that is  m-1
 int    M;     // number of elements in the base of the polynomial
 int    N;     // number of elements in the vector x
 Mat    P;     // Matrix to store the evaluations of the Polynomial 
 int    dim=2;
 Vec    X;

//get the degree of the polynomial required in the RBF
//recall that in practice, we employ the degree  m-1   
   m = rbf.get_degree_pol();      
   

//get the number of elements in the polynomial base
   M = pol.get_M();
   
//get the size of the vector
   N = x.GetSize();
    
//resize the Gramm matrix
   A.Reallocate(N + M, N + M);

//fill with zeros   
   FillZero( A, N, M );
         
//fill the matrix with the evaluations of RBFs
    T   xc,yc;
    T*  y_ = y.GetData();
    T*  x_ = x.GetData();
    
   for(int i=0; i<N; i++)
   {
      xc = x_[i];
      yc = y_[i]; 
      
     for(int j=0; j<N; j++)
     {
         A(i,j) = rbf.eval(xc,yc,x_[j],y_[j],c);
     }
   }
   
   
//if m==0 indicate that does not require a polynomial

   if(m>0)
   {
      
      compose(x,y,X);   
      
      P = pol.build_tnt(X.GetData(), N,dim);   //build the matrix P of the polynomial  evaluations 
                                              // on the points x
      
      //insert the polinomial into the Gramm matrix
      for(int i=0; i<N; i++)
      {
         for(int j=N; j<(N+M); j++)
        {   
              A(i,j) = P(i,j-N);
              A(j,i) = A(i,j);
        }   
       } // if(m>0)
   }   
   

}
//----------------------------------------------------------
// Data 2D
//----------------------------------------------------------
template<typename RBF,  typename Pol, typename Mat, typename Vec>
void fill_gram(RBF rbf,  Pol &pol, Vec &c, Vec &x, Vec &y, Mat &A)
{
//
//Build the Gramm matrix
//
//  ! A    P !  lambda   = f
//  ! P^t  0 !   mu      = 0
//
   Vec X; 
   
   compose(x,y,X);
   
   fill_gram(rbf,pol,c,X, 2, A);
}
//----------------------------------------------------------
// Data 3D
//----------------------------------------------------------
template<typename RBF, typename T, typename Pol, typename Mat, typename Vec>
void fill_gram(RBF rbf, Pol &pol,  T c, Vec &x, Vec &y, Vec &z, Mat &A)
{
//
//Build the Gramm matrix
//
//  ! A    P !  lambda   = f
//  ! P^t  0 !   mu      = 0
//
  Vector<T> X;

  compose(x,y,z,X);
 
  fill_gram(rbf,pol, c, X, 3, A);
}
//----------------------------------------------------------
// Data 3D
//----------------------------------------------------------
template<typename RBF, typename Pol, typename Mat, typename Vec>
void fill_gram(RBF rbf, Pol &pol,  Vec &c, Vec &x, Vec &y, Vec &z,Mat &A)
{
//
//Build the Gramm matrix
//
//  ! A    P !  lambda   = f
//  ! P^t  0 !  mu       = 0
//
  Vec X;

  compose(x,y,z,X);

  fill_gram(rbf, pol , c, X, 3, A);
}


} // RBF namespace

#endif // _GRAMM_HPP_

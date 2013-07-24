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
 /***************************************************************************
 * This file has been written to Radial++ in spirit to incorporate
 * the classical binary operations between vectors and matrices
 * in the library called Seldon.
 *
 * The original idea of this functions was taken from Blitz and suggested 
 * by Vivien Mallet.
 *
 ***************************************************************************/

#ifndef _MAT_OPERATORS_HPP_
#define _MAT_OPERATORS_HPP_

namespace Seldon{

   

//! Check the dimensions between Matrix vs. Matrix
/*! 
    \param two matrices.
*/
template<class Mat>
void CheckDims(const Mat &A, const Mat &B, string op)
{
  if ( ( A.GetM() != B.GetM() ) | ( A.GetN() != B.GetN() ) )
   {
         cout<<"Operation " + op + " not permitted,  ";
         cout<<"the matriz sizes are differents."<<endl;
         exit(1);
   }   
}




// Operations:
//            Matrix op T
//              T    op Matrix
//
// op = +, -, *
// T  = scalar


//!  Matrix <--  Matriz op alpha  or  alpha op Matrix, op = +,-,*
/*! 
    \param alpha scalar.
*/
#define SELDON_DEFINE_BINARY_MAT_OP(op)              \
template<typename T, typename Tc>                    \
Matrix<T> operator op (const Matrix<T>& A, Tc alpha) \
{                                                    \
   Matrix<T> C( A );                                 \
   T* data   = C.GetData();                          \
   int m     = C.GetSize();                          \
                                                     \
   for( int i=0;  i<m;  i++ )                        \
      data[i] = data[i] op alpha;                    \
                                                     \
   return C;                                         \
}                                                    \
template<typename Tc, typename T>                    \
Matrix<T> operator op (Tc alpha, const Matrix<T>& A) \
{                                                    \
   Matrix<T> C( A );                                 \
   T* data   = C.GetData();                          \
   int m     = C.GetSize();                          \
                                                     \
   for( int i=0;  i<m;  i++ )                        \
      data[i]  = data[i] op alpha;                   \
                                                     \
   return C;                                         \
}                                                    \

SELDON_DEFINE_BINARY_MAT_OP(+)
SELDON_DEFINE_BINARY_MAT_OP(-)
SELDON_DEFINE_BINARY_MAT_OP(*)






//
//!  Matrix <--  Matrix + Matrix  
//
template<typename T>
Matrix<T> operator+(const Matrix<T>& A, const Matrix<T>& B)
{

//check the data dimension
   CheckDims(A,B,"A + B");   
   
 #ifdef RADIAL_WITH_C_CODE
   int m_ = A.GetM();
   int n_ = A.GetN();
   
   Matrix<T> C( m_ , n_ );
   
   T*  a_ = A.GetData();
   T*  b_ = B.GetData();
   T*  c_ = C.GetData();
   
   for( int i=0; i<m_*n_;  i++)
     c_[i] = a_[i] + b_[i];
     
 #else
 
    Matrix<T> C( A );
    
     
//perform the plus command     
// computation of Y = Y + alpha*X
// Add(alpha, X, Y);

   Add(+1.0,B,C);
 
 #endif
    
   return C;
}   

  


//
//!  Matrix  <--  Matrix - Matrix  
//
template<typename T>
Matrix<T> operator-(const Matrix<T>& A, const Matrix<T>& B)
{
   
//check the data dimension
   CheckDims(A,B,"A - B");   
   
 #ifdef RADIAL_WITH_C_CODE
   int m_ = A.GetM();
   int n_ = A.GetN();
   
   Matrix<T> C( m_ , n_ );
   
   T*  a_ = A.GetData();
   T*  b_ = B.GetData();
   T*  c_ = C.GetData();
   
   for( int i=0; i<m_*n_;  i++)
     c_[i] = a_[i] + b_[i];
     
 #else
 
    Matrix<T> C( A );
    
     
//perform the plus command     
// computation of Y = Y + alpha*X
// Add(alpha, X, Y);

   Add(-1.0,B,C);
 
 #endif
    
   return C;
}   
  
  

//!  Matrix  <--  Matrix * Matrix  
//!
/*! 
    Perform the matrix multiplication:  A * B, internally
    check num_col(A) == num_rows(B).
    
    This is a expensive task, like O(N^3), so to increase
    the througput it is convenient to compile this library
    with the ATLAS library. Actually, only work with full
    storage matrices in RowMajor.
*/
template<typename T>
Matrix<T> operator*(const Matrix<T>& A, const Matrix<T>& B)
{
   Matrix<T> C( A.GetM() , B.GetN() );
    
//check the data dimension
  if ( ( A.GetN() != B.GetM() )  )
   {
         cout<<"Operation A*B not permitted,  ";
         cout<<"the matriz sizes are differents."<<endl;
         exit(1);
   }
     
// H  = Ux*H;
// similar method for matrix-matrix product
// Matrix<double> M(3,4), N(4,3);
// computation of A = beta*A + alpha*M*N
// MltAdd(alpha, M, N, beta, A);

   MltAdd(1.0, A, B, 0.0, C);

    
   return C;
}   



#ifdef  SELDON_WITH_OMP
//  
//! Vector <-- Matrix * Vector  with omp
//
template <class T>
void Mlt_omp(Matrix<T>& A, Vector<T>& x, Vector<T>& tmp)
{
   T* a  = A.GetData();
   T* x_ = x.GetData();
   T* y  = tmp.GetData();
   T     s;
   int   i,j,k;
   
   int m = A.GetM();
   int n = A.GetN();
   
   

   if((n%2)==0)  // is 4,6,8,10,...etc
   {
      #pragma omp parallel for private(i,j,s,k) shared(m,n,y,a,x_)
       for( i=0;  i<m;  i++)
       {
          s = 0.0;
          k = i*n;
          for( j=0;  j<n;  j+=2)
          {
              s += (a[k + j]*x_[j]) + (a[k + j +1 ]*x_[j+1])  ;
          } 
           
          y[i] = s;
       }
   }
   else{
      #pragma omp parallel for private(i,j,s,k) shared(m,n,y,a,x_)
       for( i=0;  i<m;  i++)
       {
          s = 0.0;
          k = i*n;
          for( j=0;  j<(n-1);  j+=2)
          {
              s += (a[k + j]*x_[j]) + (a[k + j +1 ]*x_[j+1])  ;
          } 
           
          y[i] = s;
       }

      for( i=0;  i<m;  i++)
       y[i] +=  (a[i*n + n-1 ]*x_[ n-1 ]);


   }

}
#endif

//  
//! Vector <-- Matrix * Vector 
//
template <class T>
inline Vector<T> operator*(Matrix<T> &A, Vector<T> &x)
{
   Vector<T>  tmp( x.GetSize() );


//check the data dimension
  if ( ( A.GetN() != x.GetSize() )  )
   {
         cout<<"Operation A*x not permitted,  ";
         cout<<"the matrix vector sizes are differents."<<endl;
         exit(1);
   }


 #ifdef  SELDON_WITH_OMP
     Mlt_omp(A , x, tmp);
 #else
     Mlt(A, x, tmp);
 #endif

  return tmp;

}


} // namespace Seldon
  
#endif //_MAT_OPERATORS_HPP_



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



#ifndef _VEC_OPERATORS_HPP_
#define _VEC_OPERATORS_HPP_





//! Check the dimensions between Vector vs. Vector
/*! 
    \param two vectors.
*/
template<typename T>
void CheckDims(const Vector<T> &A, const Vector<T> &B, string op)
{
   if (A.GetSize() != B.GetSize() )
   {
        cout<<"Operation " + op + " not permitted,  ";
         cout<<"the vector sizes are differents."<<endl;
         exit(1);
   }   
}  
  

// Operations:
//            Vector op T
//              T    op Vector
//
// op = +, -, *
// T  = scalar


//!  Vector <--  Vector op alpha  or  alpha op Vector, op = +,-,*
/*! 
    \param alpha scalar.
*/
#define SELDON_DEFINE_BINARY_OP(op)                  \
template<typename T, typename Tc>                    \
Vector<T> operator op (const Vector<T>& X, Tc alpha) \
{                                                    \
   Vector<T> C( X );                                 \
   T* data   = C.GetData();                          \
   int m     = C.GetSize();                          \
                                                     \
   for( int i=0;  i<m;  i++ )                        \
        data[i] = data[i] op alpha;                  \
                                                     \
   return C;                                         \
}                                                    \
template<typename Tc, typename T>                    \
Vector<T> operator op (Tc alpha, const Vector<T>& X) \
{                                                    \
   Vector<T> C( X);                                  \
   T* data   = C.GetData();                          \
   int m     = C.GetSize();                          \
                                                     \
   for( int i=0;  i<m;  i++ )                        \
      data[i]  = data[i] op alpha;                   \
                                                     \
   return C;                                         \
}                                                    \

SELDON_DEFINE_BINARY_OP(+)
SELDON_DEFINE_BINARY_OP(-)
SELDON_DEFINE_BINARY_OP(*)



//
//! Vector <-- Vector + Vector
//
template<typename T>
inline Vector<T> operator+(const Vector<T>& X,const Vector<T>& Y)
{   
//check the data dimension
   CheckDims(X,Y,"X + Y");
   

 #ifdef RADIAL_WITH_C_CODE
   
   Vector<T> C( X.GetSize() );
    
   T* x_  = X.GetData();
   T* y_  = Y.GetData();
   T* c_  = C.GetData();
   int n  = X.GetSize();
   
   
   for(int i=0;  i<n;  i++)
       c_[i] = x_[i]+y_[i];      
  
 #else 
     Vector<T> C( X );
       
   //perform the plus command     
   // computation of Y = Y + alpha*X
   // Add(alpha, X, Y);

      Add(+1.0,Y,C);
 #endif  

   return C;
}   






//
//! Vector <-- Vector - Vector
//
template<typename T>
inline Vector<T> operator-(const Vector<T>& X,const Vector<T>& Y)
{   
//check the data dimension
   CheckDims(X,Y,"X - Y");
   

  #ifdef RADIAL_WITH_C_CODE
  
    Vector<T> C( X.GetSize() );
    
   T* x_  = X.GetData();
   T* y_  = Y.GetData();
   T* c_  = C.GetData();
   int n  = X.GetSize();
   
   
   for(int i=0;  i<n;  i++)
       c_[i] = x_[i]-y_[i];      
  
 #else
 
      Vector<T> C( X );
       
   //perform the plus command     
   // computation of Y = Y + alpha*X
   // Add(alpha, X, Y);

    Add(-1.0,Y,C);

 #endif  

   return C;
}   

 
 

//
//! Vector <-- Vector.*Vector
//
template<typename T>
Vector<T> pdot(const Vector<T> &X, const Vector<T> &Y)
{
   Vector<T>  C( X.GetSize() );
   T* ptr_x   = X.GetData();
   T* ptr_y   = Y.GetData();
   T* ptr_c   = C.GetData();
   int m      = C.GetSize();
   
//check the data dimension
  if ( ( X.GetSize() != Y.GetSize() )  )
   {
         cout<<"Operation pdot(X,Y) not permitted,  ";
         cout<<"the vector sizes are differents."<<endl;
         exit(1);
   }   
    
   for( int i=0;  i<m;  i++ )
      ptr_c[i] = ptr_x[i] * ptr_y[i];
      
   return C;  
}
 

  
#endif //_VEC_OPERATORS_HPP_



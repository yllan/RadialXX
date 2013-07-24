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


#ifndef _SOLVER_HPP_ 
#define _SOLVER_HPP_ 


#include "jag_cond_all.cpp"


namespace solver {


template <typename T>
class LU
{
   Matrix<T>  A; 
    Vector<int> pivot;
   bool initialized;   
   int  n_eq;
   int op;
   public:
      LU(void);
      ~LU(void);
   T          factorize(Matrix<T> &B, int op=0);
   Vector<T>  solve(Vector<T> &b);
   Vector<T>  solve(Matrix<T> &B, Vector<T>  &b);
   Vector<T>  gauss(Matrix<T> &B, Vector<T>  &b);

};
//----------------------------------------------------------
template <typename T> 
LU<T>::LU(void)
{
   pivot.Resize(0);
   n_eq        = 0;
   initialized = false;
}
//----------------------------------------------------------
template <typename T> 
LU<T>::~LU(void)
{
   n_eq        = 0;
   initialized = false;
   
   pivot.Resize(0);
   A.Resize(0,0);
}   
//----------------------------------------------------------
template <typename T> 
T LU<T>::factorize(Matrix<T>  &B, int op)
{
  T cond;

   
   if(op==0) //LU factorization over a copy of B
   {

      initialized = true;
      this->op    = 0;
      
      n_eq = B.GetM();
      
      A.Resize(n_eq,n_eq);
     
      pivot.Reallocate(n_eq);
      
      A.Copy(B);

      #ifdef SELDON_WITH_LAPACK
         GetLU(A, pivot);
         cond = -1;
      #else
         Factor_t(A, n_eq, &cond, pivot);
      #endif
      
      return cond;
   }
   else if(op==1){
      
      
      initialized = true;
      this->op    = 1;
      
      n_eq = B.GetM();
      
      pivot.Reallocate(n_eq);
      
      
      #ifdef SELDON_WITH_LAPACK
      
         GetLU(B, pivot);
         
         cond = -1;
         
      #else
         
         Factor_t(B, n_eq, &cond, pivot);
    
         printf("condition number = %e\n",cond);
         
      #endif
      
      return cond;      

   }   
   else{
    
      fprintf(stderr,"\nERROR: invalid op=%d in LU::factorize, op = {0,1}.\n\n",op);
      fflush(stderr);
      exit(1);
      
   }   

}
//----------------------------------------------------------
template <typename T> 
Vector<T> LU<T>::solve(Vector<T>  &b)
{
 Vector<T> lambda;   

   
  if(initialized==false)
  {    
   fprintf(stderr,"\nERROR: in LU::solve,  you must be call first LU::factorize.\n\n");
   fflush(stderr);
   exit(1);   
  }
   
 if(op==0)
 {    
    lambda.Resize(n_eq);
    
    lambda.Copy(b);

     #ifdef SELDON_WITH_LAPACK
           SolveLU(A, pivot, lambda);
     #else
           Solve_t(A, n_eq, pivot, lambda);
     #endif

    
    return lambda;
 }    
 else{
      fprintf(stderr,"\nERROR: you must be use the funcion LU::solve(Matrix,Vector);\n\n");
      fflush(stderr);
      exit(1);
 }   
    
}
//----------------------------------------------------------
template <typename T> 
Vector<T> LU<T>::solve(Matrix<T> &B, Vector<T>  &b)
{
   Vector<T> lambda;   
   
   
   if(initialized==false)
   {    
      fprintf(stderr,"\nERROR: in LU::solve,  you must be call first LU::factorize.\n\n");
      fflush(stderr);
      exit(1);   
   }   

//validate is B is square
  if(B.GetM() != B.GetN()){
     printf("\nERROR: in solve function the input matrix must be square.\n\n");
     exit(1);   
  }
  
//validate the dimension of the input vectos  
  if(B.GetM() != b.GetSize()){
     printf("\nERROR: in solve function the input matrix and vector have different sizes.\n");
     printf("          size(A) = %d x %d   size(b) = %d\n\n",B.GetM(),B.GetN(),b.GetSize());
     exit(1);   
  }
   
   if(op==1)
   {    
      
      n_eq = B.GetM();
      
      lambda.Reallocate(n_eq);
      
      lambda = b; //lambda.Copy(b);
      
        
        #ifdef SELDON_WITH_LAPACK
             SolveLU(B, pivot, lambda);  
       #else
            Solve_t2(B, n_eq, pivot, lambda);
       #endif
     
      return lambda;
   }       

   printf("\n!!ERROR in lu.solve(A,b), op = %d\n",op);
   printf("  in file '%s' near to line %d\n",__FILE__,__LINE__);
   exit(1);
   return lambda;
}
//----------------------------------------------------------
template <typename T> 
Vector<T> LU<T>::gauss(Matrix<T> &B, Vector<T>  &b)
{
  Matrix<T>   Aux(B);
  Vector<T>   lambda(b);
  Vector<int> piv;
  int         n;

   

//validate is B is square
  if(B.GetM() != B.GetN()){
     printf("\nERROR: in gauss function the input matrix must be square.\n\n");
     exit(1);   
  }
  
//validate the dimension of the input vectos  
  if(B.GetM() != b.GetSize()){
     printf("\nERROR: in gauss function the input matrix and vector have different sizes.\n");
     printf("          size(A) = %d x %d   size(b) = %d\n\n",B.GetM(),B.GetN(),b.GetSize());
     exit(1);   
  }

   n = Aux.GetM();
   
   piv.Reallocate(n);
   
   #ifdef SELDON_WITH_LAPACK
          GetLU(Aux,pivot);
          SolveLU(Aux, pivot, lambda); 
   #else
          T           cond;
          int         error;
         
        error = Factor_t(Aux, n, &cond, piv);   
        
        printf("condition number = %e\n",cond);
        
        if(error!=0){
          printf("\nERROR in the LU factorization, error = %d\n\n",error);
          exit(1);   
        }

        Solve_t(Aux, n, piv, lambda);
   #endif

   return lambda;
}
//----------------------------------------------------------
//  Auxiliar Functions
//----------------------------------------------------------
template <typename T> 
Vector<T> gauss(Matrix<T> &A, Vector<T>  &b)
{
  LU<T> mylu;
  
  return mylu.gauss(A,b);
}


}

#endif // _SOLVER_HPP_ 

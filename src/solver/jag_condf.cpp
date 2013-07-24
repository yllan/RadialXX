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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef _JAG_COND_CPP_ 2
#define _JAG_COND_CPP_ 2


void Solve_f(float **a,int  neq,int * pivot_index, float * b);
int Factor_f(float **a, int neq, float * cond, int *pivot_index);

/* FUNCTION:     This is the file containing C auxiliary math functions.
   AUTHORS:      Lawrence Shampine, Richard Allen, Steven Pruess for
                 the text  Fundamentals of Numerical Computing
   DATE:         December 8, 1995
   LAST CHANGE:  July 3, 1996                                      */

/*/////////////////////////////////////////////////////////////////////////////
 //    This file contains utilities for common mathematical operations.
 /////////////////////////////////////////////////////////////////////////////*/

/* Return the absolute value of the double value x.     */
inline float abs_f(float x)
{
   if (x >= 0)
      return (x);
   else
      return (-x);
}

inline float max_f(float x, float y)
/* Return the maximum of double values x and y.      */
{
   if (x >= y)
      return (x);
   else
      return (y);
}

inline float min_f(float x, float y)
/* Return the minimum of double values x and y.     */
{
   if (x <= y)
      return (x);
   else
      return (y);
}

//-------------------------------------------------------------------------------------------------
int Factor_f(float **a, int neq, float * cond, int *pivot_index)
{
   /*
    Function Factor decomposes the matrix A using Gaussian elimination and
    estimates its condition number;  Factor may be used in conjunction with
    function Solve to solve A*x=b.
    
    Input parameters:
     a       = matrix A to be triangularized.
     neq     = number of equations to be solved.
    Output parameters:
     a       = the upper triangular matrix U in its upper portion (by
    rows) and a permuted version of a lower triangular
    matrix I-L such that (permutation matrix)*A = L*U;
    a record of interchanges is kept in pivot_index.
     return_value = an integer variable that reports whether or not the
    matrix A has a zero pivot.  A value of zero means no zero
    pivot was encountered; if positive, the first zero pivot
    occurred at the return_value's equation (1 to neq) and 
    the decomposition could not be completed.  If the
    return_value is -1 then there is an input error (neq not 
    positive).  If the return_value is -2 then there was
    insufficient memory for the code to continue.
     cond    = an estimate of the condition number of A (unless the
    return_value is nonzero).
     pivot_index = the pivot vector which keeps track of row interchanges;
    also,   pivot_index[neq-1] = (-1)**(number of interchanges).
    
    The determinant of A can be obtained on output from
    det(A) = pivot_index[neq-1] * a[0][0] *  ... * a[neq-1][neq-1].
    
    Declare local variables and initialize:    */
    float anorm, dnorm, ek, t, *temp, ynorm;
    int i, j, k, m;
    static float zero = 0.0, one = 1.0;
   
   
   
    if (neq <= 0) return -1;
   
    *cond = zero;
   
    pivot_index[neq-1] = 1;
   
    if (neq == 1)
    {
      /*
       neq = 1 is a special case.
       */
      if (a[0][0] == zero) return 1;
      else
      {
         *cond = one;
         return 0;
      }
    }
   /*
    Compute infinity-norm of A for later condition number estimation.
    */
    anorm = zero;
    for (i = 0; i < neq; i++)
    {
      
      t = zero;
      for (j = 0; j < neq; j++)
         t += abs_f(a[i][j]);
      anorm = max_f(t, anorm);
    }
   /*
    Gaussian elimination with partial pivoting.
    */
    for (k = 0; k < neq-1; k++)
    {
      
      /*
       Determine the row m containing the largest element in
       magnitude to be used as a pivot.
       */
      m = k;
      for (i = k+1; i < neq; i++)
         if (abs_f(a[i][k]) > abs_f(a[m][k])) m = i;
      /*
       Check for a nonzero pivot; if all possible pivots are zero,
       matrix is numerically singular.
       */
      if (a[m][k] == zero)
         return k + 1;
      pivot_index[k] = m;
      if (m != k)
      {
         /*
          Interchange the current row k with the pivot row m.
          */
         pivot_index[neq-1] = -pivot_index[neq-1];
         for (j = k; j < neq; j++)
         {
            t = a[m][j];
            a[m][j] = a[k][j];
            a[k][j] = t;
         }
      }
      /*
       Eliminate subdiagonal entries of column k.
       */
      for (i = k+1; i < neq; i++)
      {
         t = a[i][k]/a[k][k];
         a[i][k] = -t;
         if (t != zero)
         {
            for (j = k+1; j < neq; j++)
               a[i][j] -= t*a[k][j];
         }
      }
    }
    if (a[neq-1][neq-1] == zero)
      return neq;
   /*
    Estimate the condition number of A by computing the infinity
    norm of A directly and a lower bound for the norm of A-inverse.
    A lower bound for the norm of A-inverse is provided by the ratio
    norm(Y)/norm(D) for any vectors such that A*Y = D and D != 0.
    A "large" ratio is obtained by computing Y as one iteration of
    inverse iteration for the smallest singular value of A, i.e.,
    by solving for Y such that (transpose(A)*A)*Y = E.  This exploits
    the fact that an LU decomposition of A can be used to solve the
    linear system transpose(A)*D = E as well as A*Y = D.  The entries
    of E are +1 or -1 with the sign chosen during the computation of
    D to increase the size of the entry of D and so make a "large"
    lower bound for the norm of A-inverse more likely.
    
    First, allocate some space for a temporary array.
    */
    temp = (float*) malloc(neq*sizeof(float));
    if (temp == NULL) return -2;
   
    temp[0] = -1.0/a[0][0];
    for (k = 1; k < neq; k++)
    {
        t = 0.0;
        for (i = 0; i < k; i++) t = t+a[i][k]*temp[i];
        if (t < 0)
         ek = -1.0;
        else
         ek = 1.0;
      
        temp[k] = -(ek+t)/a[k][k];
    }
    for (k = neq-2; k >= 0; k--)
    {
        t = 0.0;
        for (i = k+1; i < neq; i++) t = t+a[i][k]*temp[i];
        temp[k] += t;
        m = pivot_index[k];
        t = temp[m];
        temp[m] = temp[k];
        temp[k] = t;
    }
    
    dnorm = zero;
    for (i = 0; i < neq; i++)
      dnorm = max_f(dnorm, abs_f(temp[i]));
    Solve_f(a, neq, pivot_index, temp);
    ynorm = zero;
    for (i = 0; i < neq; i++)
      ynorm = max_f(ynorm, abs_f(temp[i]));
    *cond = anorm*ynorm/dnorm;
    free(temp);
   
   //  *cond=-1;
   return 0;
}


/* ------------------------------------------------------------------ */
void Solve_f(float **a,int  neq,int * pivot_index, float * b)
{
   /*
    Function Solve solves the linear system A*x=b using the factorization
    obtained from Factor.  Do not use Solve if a zero pivot has
    been detected in Factor.
    
    Input parameters:
     a           = an array returned from Factor containing the
    triangular decomposition of the coefficient matrix.
     neq         = number of equations to be solved.
     pivot_index = vector of information about row interchanges
    obtained from Factor.
     b           = right hand side vector b.
    Output parameters:
     b           = solution vector x.
    
    Local variables:    */
   int       m;
   register int i,j,k;
   float t;
   float aux;
   /*
    Forward elimination.
    */
   if (neq > 1)
   {
      for (k = 0; k < neq-1; k++)
      {
         m = pivot_index[k];
         t = b[m];
         b[m] = b[k];
         b[k] = t;
         for (i = k+1; i < neq; i++)
            b[i]+= a[i][k]*t;
         
      }
      /*
       Back substitution.
       */
      for (i = neq-1; i >= 0; i--)
      {
         aux=b[i];
         for (j = i+1; j < neq; j++)
                aux -= a[i][j]*b[j];
         b[i]=aux;
         b[i] = b[i]/a[i][i];
      }
   }
   else
      b[0] = b[0]/a[0][0];
}

#endif // _JAG_COND_CPP_ 2

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

#ifndef _INTERPOLATE_HPP_
#define _INTERPOLATE_HPP_

namespace rbf{
	
//----------------------------------------------------------
// Validate that pol.dimension equal to argin in the function
//----------------------------------------------------------
template<class Pol>
void check_pol_dimension(Pol& pol, int dim)
{
   if( pol.get_d() != dim ){
      fprintf(stderr,"\n");
      fprintf(stderr,"ERROR: in function interpolate the pol.dimension=%d must be %d.\n",pol.get_d(),dim);
      fprintf(stderr,"\n");
      exit(1);
   }	
}	


//----------------------------------------------------------
// Interpolate Data N-dimensional
//----------------------------------------------------------
template <typename RBF, typename Vec, typename T>
Vec interpolate(RBF rbf, Polinomio<T> &pol, Vec& c, Vec &lambda,  Vec &x, Vec &x_new, int dim)
{
 Vec   px;
 Vec   z_new;
 int   m,dim_px;
 int   N,Nnew;
 T     s;   
   
  if(dim>1){
     if((x.GetSize()%dim)>0){
       fprintf(stderr,"ERROR:  N = %d must be divisible by dim = %d.\n",x.GetSize(),dim);
       fprintf(stderr,"in '%s' near to line %d\n\n",__FILE__,__LINE__);
       exit(1);
     }
    if((x_new.GetSize()%dim)>0){
       fprintf(stderr,"ERROR:  Nnew = %d must be divisible by dim = %d.\n",x_new.GetSize(),dim);
       fprintf(stderr,"in '%s' near to line %d\n\n",__FILE__,__LINE__);
       exit(1);
     }        
   }   
   
   assert( c.GetSize()==(x.GetSize()/dim) );
   
   check_pol_dimension(pol,dim);	
   
   
   N    = x.GetSize()/dim;
   
   Nnew = x_new.GetSize()/dim;
   
   z_new.Reallocate(Nnew);
   
   for(int i=0; i<Nnew; i++)
   {
      s = 0.0;
      for(int j=0; j<N; j++)
      {
         s = s + lambda(j)*rbf.eval(&x_new(i*dim),&x(j*dim),dim,c(j));
      }
      
      z_new(i) = s;   
   }      
   
//get the degree of the associated polynomil to the RBF
//recall the degree is at most m-1   
   m = pol.get_M();
   
   if(m>0) //require a polynomial in the interpolation
   {
      dim_px = pol.get_M();
          
      assert( lambda.GetSize()==(x.GetSize()/dim + dim_px) );
      
      px.Reallocate(dim_px);
      
      for(int i=0; i<Nnew; i++)
      {
         pol.eval(&x_new(i*dim),dim,px.GetData(),dim_px);
         
         s = 0.0;
         
         for(int k=0; k<dim_px; k++)
         {
            s = s + lambda(k+N)*px(k);   
         }      
         
         z_new(i) = z_new(i) + s;   
      }         
   }
   
   return z_new;      
}
//----------------------------------------------------------
// Data 1D
//----------------------------------------------------------
template <typename RBF, typename Vec, typename T, typename Pol>
Vec interpolate(RBF rbf, Pol &pol, T c, Vec &lambda,  Vec &x, Vec &x_new)
{
   Vector<T>  myc;
  
   myc.Reallocate(x.GetSize());
  
   myc = c;

   return interpolate(rbf, pol, myc, lambda, x, x_new, 1);
  
}
//----------------------------------------------------------
// Data 1D
//----------------------------------------------------------
template <typename RBF, typename Vec, typename T, typename Pol>
Vec interpolate(RBF rbf, Pol &pol, Vec& c, Vec &lambda, Vec &x, Vec &x_new)
{
   return interpolate(rbf, pol, c, lambda, x, x_new, 1);
}
//----------------------------------------------------------
// Data 2D
//----------------------------------------------------------
template <typename RBF, typename Vec, typename T,typename Pol>
Vec interpolate(RBF rbf,  Pol &pol, T c, Vec &lambda, Vec &x, Vec &y, Vec &x_new, Vec &y_new)
{
   Vector<T> X,Y,Z,myc;

   compose(x,y,X);
   compose(x_new,y_new,Y);      

   myc.Reallocate(x.GetSize());
   myc = c;

 return interpolate(rbf, pol, myc, lambda, X, Y , 2);
}
//----------------------------------------------------------
// Data 2D
//----------------------------------------------------------
template <typename RBF, typename Vec, typename Pol>
Vec interpolate(RBF rbf, Pol &pol, Vec& c, Vec &lambda, Vec &x, Vec &y, Vec &x_new, Vec &y_new)
{
   Vec  X,Y;

   compose(x,y,X);
   compose(x_new,y_new,Y);   

 return interpolate(rbf, pol, c, lambda, X, Y , 2);
}
//----------------------------------------------------------
// Data 3D
//----------------------------------------------------------
template <typename RBF, typename Vec, typename T, typename Pol>
Vec interpolate(RBF rbf,  Pol &pol, T c, Vec &lambda,  Vec &x, Vec &y, Vec &z, Vec &x_new, Vec &y_new, Vec &z_new)
{
   Vector<T> X,Y,Z,myc;

   compose(x,y,z,X);
   compose(x_new,y_new,z_new,Y);      

   myc.Reallocate(x.GetSize());
   myc = c;

 return interpolate(rbf,  pol, myc, lambda,  X, Y , 3);
}
//----------------------------------------------------------
// Data 3D
//----------------------------------------------------------
template <typename RBF, typename Vec,typename Pol>
Vec interpolate(RBF rbf, Pol &pol, Vec &c, Vec &lambda,Vec &x, Vec &y, Vec &z, Vec &x_new, Vec &y_new, Vec &z_new)
{
   Vec X,Y,Z;
   
   compose(x,y,z,X);
   
   compose(x_new,y_new,z_new,Y);   
   
 return interpolate(rbf, pol, c, lambda,  X, Y , 3);
}

//
// The next three functions are interpolators for one point of data in 1D,2D and 3D only.
// They are neccessary to provide a fast implementations for this cases.
//

//----------------------------------------------------------
// Data 1D
//----------------------------------------------------------
//
//   Interpolate only one point in 1D
//
template <typename RBF, typename Vec, typename T,typename Pol>
T interpolate(RBF rbf,  Pol &pol, T c, Vec &lambda, Vec &x,  T xc)
{
   Vec   px;   
   int   m,dim_px,N;   
   T     xnew[1]; 
   T     s;
   int   j;
   
   assert( lambda.GetSize()==(x.GetSize()+pol.get_M()) );
   
   check_pol_dimension(pol,1);	
   
   N    = x.GetSize();
   
   s    = 0.0;
   
   for(j = 0; j < N; j++)
   {
      s = s +  lambda(j)*rbf.eval(xc,x(j),c);
   }
   
   m = pol.get_M();
      
   if(m>0) //require a polynomial in the interpolation
   { 
      dim_px = pol.get_M();
      
      px.Reallocate(dim_px);
      
      xnew[0] = xc;
        
      pol.eval(&xnew[0],1,px.GetData(),dim_px);
         
      for(int k=0; k<dim_px; k++)
      {
         s = s + lambda(k+N)*px(k);   
      }      
   }
   
 return s;
}


//----------------------------------------------------------
// Data 1D
//----------------------------------------------------------
//
//   Interpolate only one point in 1D
//   with a c_vector as shape parameter
//
template <typename RBF, typename Vec, typename T,typename Pol>
T interpolate(RBF rbf, Pol &pol, Vec& c, Vec &lambda, Vec &x,  T xc)
{
   Vec   px;   
   int   m,dim_px,N;   
   T     xnew[1]; 
   T     s;
   
   
   assert( lambda.GetSize()==(x.GetSize()+pol.get_M()) );
   
   assert( c.GetSize()==x.GetSize() );
   
   
   check_pol_dimension(pol,1);
      
   N    = x.GetSize();
   
   s    = 0.0;
   for(int j=0; j<N; j++)
   {
     s = s +  lambda(j)*rbf.eval(xc,x(j),c(j));
   }
   
   m = pol.get_M();    
   
   if(m>0) //require a polynomial in the interpolation
   { 
      dim_px = pol.get_M();
      
      px.Reallocate(dim_px);
      
      xnew[0] = xc;
        
      pol.eval(&xnew[0],1,px.GetData(),dim_px);
         
      for(int k=0; k<dim_px; k++)
      {
         s = s + lambda(k+N)*px(k);   
      }      
   }
   
 return s;
}


//----------------------------------------------------------
// Data 2D
//----------------------------------------------------------
//
//   Interpolate only one point in 2D
//
template <typename RBF, typename Vec, typename T,typename Pol>
T interpolate(RBF rbf, Pol& pol, T c, Vec& lambda, Vec& x, Vec& y, T xc, T yc)
{
   Vec   px;   
   int   m,dim_px,N;   
   T     xnew[2]; 
   T     s;
   
   
   assert( lambda.GetSize()==(x.GetSize()+pol.get_M()) );
   
   assert( x.GetSize()==y.GetSize() );
   
   check_pol_dimension(pol,2);	
   
   N    = x.GetSize();
   
   s    = 0.0;
   
   for(int j=0; j<N; j++)
   {
       s = s +  lambda(j)*rbf.eval(xc,yc,x(j),y(j),c);
   }
   
   m = pol.get_M();
      
   if(m>0) //require a polynomial in the interpolation
   { 
       dim_px = pol.get_M();
      
       px.Reallocate(dim_px);
      
       xnew[0] = xc;
       xnew[1] = yc;
        
       pol.eval(&xnew[0],2,px.GetData(),dim_px);
         
       for(int k=0; k<dim_px; k++)
       {
          s = s + lambda(k+N)*px(k);   
       }      
   }
   
 return s;
}


//----------------------------------------------------------
// Data 2D
//----------------------------------------------------------
//
//   Interpolate only one point in 2D
//   with a c_vector as shape parameter
//
template <typename RBF, typename Vec, typename T,typename Pol>
T interpolate(RBF rbf,  Pol &pol, Vec &c, Vec &lambda,  Vec &x, Vec &y, T xc, T yc)
{
   Vec   px;   
   int   m,dim_px,N;   
   T     xnew[2]; 
   T     s;
   
   assert( lambda.GetSize()==(c.GetSize()+pol.get_M()) );
   
   assert( c.GetSize()==x.GetSize() );
  
   assert( x.GetSize()==y.GetSize() );
   
   check_pol_dimension(pol,2);	
    
   N    = x.GetSize();
   
   s    = 0.0;
   for(int j=0; j<N; j++)
   {
      s = s +  lambda(j)*rbf.eval(xc,yc,x(j),y(j),c(j));
   }
   
   m = pol.get_M();
      
   if(m>0) //require a polynomial in the interpolation
   { 
       dim_px = pol.get_M();
      
       px.Reallocate(dim_px);
      
       xnew[0] = xc;
       xnew[1] = yc;
        
       pol.eval(&xnew[0],2,px.GetData(),dim_px);
         
       for(int k=0; k<dim_px; k++)
       {
         s = s + lambda(k+N)*px(k);   
       }      
   }
   
 return s;
}


//----------------------------------------------------------
// Data 3D
//----------------------------------------------------------
//
//   Interpolate only one point in 3D
//
//
template <typename RBF, typename Vec, typename T,typename Pol>
T interpolate(RBF rbf,  Pol &pol, T c, Vec &lambda, Vec &x, Vec &y, Vec &z, T xc, T yc, T zc)
{
   Vec   px;   
   int   m,dim_px,N;   
   T     xnew[3]; 
   T     s;
   
   
   assert( lambda.GetSize()==(x.GetSize()+pol.get_M()) );
   
   assert( x.GetSize()==y.GetSize() );
   
   assert( x.GetSize()==z.GetSize() );
   
   assert( y.GetSize()==z.GetSize() );
  
   check_pol_dimension(pol,3);	
   
   N    = x.GetSize();
   
   s    = 0.0;
   
   for(int j=0; j<N; j++)
   {
      s = s +  lambda(j)*rbf.eval(xc,yc,zc,x(j),y(j),z(j),c);
   }
   
   m = pol.get_M();
     
   if(m>0) //require a polynomial in the interpolation
   { 
       dim_px = pol.get_M();
      
       px.Reallocate(dim_px);
      
       xnew[0] = xc;
       xnew[1] = yc;
       xnew[2] = zc;
        
       pol.eval(&xnew[0],3,px.GetData(),dim_px);
         
       for(int k=0; k<dim_px; k++)
       {
          s = s + lambda(k+N)*px(k);   
       }      
   }
   
 return s;
}

//----------------------------------------------------------
// Data 3D
//----------------------------------------------------------
//
//   Interpolate only one point in 3D
//   with a c_vector as shape parameter
//
template <typename RBF, typename Vec, typename T,typename Pol>
T interpolate(RBF rbf, Pol &pol, Vec &c, Vec &lambda, Vec &x, Vec &y, Vec &z, T xc, T yc, T zc)
{
   Vec   px;   
   int   m,dim_px,N;   
   T     xnew[3]; 
   T     s;
   
   assert( lambda.GetSize()==(c.GetSize()+pol.get_M()) );
   
   assert( c.GetSize()==x.GetSize() );
      
   assert( x.GetSize()==y.GetSize() );
   
   assert( x.GetSize()==z.GetSize() );
   
   assert( y.GetSize()==z.GetSize() );
   
   check_pol_dimension(pol,3);	
   
   N    = x.GetSize();
   
   s    = 0.0;
   for(int j=0; j<N; j++)
   {
      s  = s + lambda(j)*rbf.eval(xc,yc,zc,x(j),y(j),z(j),c(j));
   }
   
   m = pol.get_M();
     
   
   if(m>0) //require a polynomial in the interpolation
   { 
       dim_px = pol.get_M();
      
       px.Reallocate(dim_px);
      
       xnew[0] = xc;
       xnew[1] = yc;
       xnew[2] = zc;
        
       pol.eval(&xnew[0],3,px.GetData(),dim_px);
         
       for(int k=0; k<dim_px; k++)
       {
          s = s + lambda(k+N)*px(k);   
       }      
   }
   
 return s;
}



} // RBF namespace

#endif // _INTERPOLATE_HPP_

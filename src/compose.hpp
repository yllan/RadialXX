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

#ifndef _COMPOSE_HPP_
#define _COMPOSE_HPP_

namespace rbf{


//----------------------------------------------------------
//  Data 2D -->  1D 
//
//  (x,y)   -> X=(x1,y1,x2,y2,...)
//
//----------------------------------------------------------
template <typename Vec>
void   compose(Vec& x, Vec& y, Vec& X)
{
   int i,n;
   
   if(x.GetSize()!=y.GetSize())
   {
      fprintf(stderr,"ERROR: in compose the vectors have differents sizes.\n");
      fprintf(stderr,"\n");
      exit(1);
   }
  
   n = x.GetSize();

   X.Reallocate(2*n);

   for( i=0;  i<n;  i++ )
   {
      X(i*2+0) = x(i);
      X(i*2+1) = y(i);
   }

}
//----------------------------------------------------------
//  Data 3D -->  1D 
//
//  (x,y,z) -> X=(x1,y1,z1,x2,y2,z2,...)
//
//----------------------------------------------------------
template <typename Vec>
void   compose(Vec &x, Vec &y, Vec &z, Vec &X)
{
   int i,n;

   if( ( x.GetSize()!=y.GetSize() ) || ( x.GetSize()!=z.GetSize() ) || ( y.GetSize()!=z.GetSize() )  )
   {
      fprintf(stderr,"ERROR: in compose the vectors have differents sizes.\n");
      fprintf(stderr,"\n");
      exit(1);
   }
   
   n = x.GetSize();

   X.Reallocate(3*n);

   for( i=0;  i<n;  i++ )
   {
       X(i*3+0) = x(i);
       X(i*3+1) = y(i);
       X(i*3+2) = z(i);
   }
}

} // RBF namespace

#endif //_COMPOSE_HPP_



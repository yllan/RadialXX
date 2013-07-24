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

#ifndef _FILL_SELECT_HPP_
#define _FILL_SELECT_HPP_
 
namespace rbf{

#define  RBF_NORMAL    0
#define  RBF_DX        1
#define  RBF_DY        2
#define  RBF_DZ        3 
#define  RBF_DXX       4
#define  RBF_DYY       5
#define  RBF_DZZ       6
#define  RBF_DXY       7
#define  RBF_DYX       8
#define  RBF_DXXX      9
#define  RBF_DYYY      10
#define  RBF_DZZZ      11 
#define  POL_TRANSPUESTO 12
    
int select(string type)
{
   if(type=="normal")
      return RBF_NORMAL;
   
   if(type=="pol_trans")
      return POL_TRANSPUESTO;    
   
   if(type=="dx")
      return RBF_DX;
   
   if(type=="dy")
      return RBF_DY;
   
   if(type=="dz")
      return RBF_DZ;    
      
   if(type=="dxy")
      return RBF_DXY;   
      
   if(type=="dyx")
      return RBF_DYX;   
   
   if(type=="dxx")
      return RBF_DXX;   
   
   if(type=="dyy")
      return RBF_DYY; 
   
   if(type=="dzz")
      return RBF_DZZ;     
      
   if(type=="dxxx")
      return RBF_DXXX;   
   
   if(type=="dyyy")
      return RBF_DYYY; 
   
   if(type=="dzzz")
      return RBF_DZZZ;    
   
   return -1;  
}

} // RBF namespace

#endif //_FILL_SELECT_HPP_

//
// Future:
//       In 2d data we have the following possibilities:
//       dx   dy
//       dxx  dyy
//       dxy  dyx   ?equal,  not implemented
//
//       In 3d data we have more possibilities:
//       dx    dy   dz
//       dxx   dyy  dzz
//       dxxx  dyyy dzzz
//       dxy   dxz  dyz   
//       dxyy  dxxz    etc...
//       
//     In conclusion, we need a better form to incorporate the definition
//     of kernel derivatives.
//
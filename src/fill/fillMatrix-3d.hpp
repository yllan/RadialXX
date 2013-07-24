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

#ifndef _FILL_3D_HPP_
#define _FILL_3D_HPP_

namespace rbf{
   

//----------------------------------------------------------
// Data 3D
//---------------------------------------------------------- 
template <typename RBF, typename Vec, typename T, typename Pol, typename Mat>
void fill_matrix(string type,               // select the operation to perform
               RBF rbf, Pol &pol,           // defined the configuration of RBF + Pol 
               Vector<T>&  c, 
               Vec &x, Vec &y, Vec& z,     // three input vectors
               int ini, int fin,           // from row=ini to row<fin in A 
               Mat &A                      // output
               )
{
  int  n,dim,dim_px;
  int  i,j;
  T    xc[3];
  T*   y_ = y.GetData();
  T*   x_ = x.GetData();  
  T*   z_ = z.GetData();  
  T    xc_,yc_,zc_;
  
  
  assert( ini<=fin );
  
//stablish the data dimension
  dim = 3;
  
  
//get the vector dimension
  n      = x.GetSize();
  
//obtain the number of elements in the polynomial base
  dim_px = pol.get_M();

  switch(select(type))
  {
    case   RBF_NORMAL:
                       for( i=ini;  i<fin;  i++ )
                       { 
                          xc_ = x_[i];
                          yc_ = y_[i];
                          zc_ = z_[i];
                       
                          for( j=0; j<n;  j++ ) 
                             A(i,j) = rbf.eval(xc_,yc_,zc_,x_[j],y_[j],z_[j],c[j]);
                       }     
                   break;   
        
    case   POL_TRANSPUESTO:
                  break;
        
    case   RBF_DX: 
                   for( i=ini;  i<fin; i++ ) 
                   {
                       xc_ = x_[i];
                       yc_ = y_[i];
                       zc_ = z_[i];                      
                       
                       for( j=0;  j<n;  j++ ) 
                          A(i,j) = rbf.dx(xc_,yc_,zc_,x_[j],y_[j],z_[j],c[j]);
                   }
                   break;
        
    case   RBF_DY: 
                   for( i=ini;  i<fin; i++ ) 
                   {
                       xc_ = x_[i];
                       yc_ = y_[i];
                       zc_ = z_[i];  
                       
                       for( j=0;  j<n;  j++ ) 
                          A(i,j) = rbf.dy(xc_,yc_,zc_,x_[j],y_[j],z_[j],c[j]);
                   }                         
                         
                         
                   break;
                   
        
    case   RBF_DZ: 
                   for( i=ini;  i<fin; i++ ) 
                   {
                       xc_ = x_[i];
                       yc_ = y_[i];
                       zc_ = z_[i];  
                       
                       for( j=0;  j<n;  j++ ) 
                          A(i,j) = rbf.dz(xc_,yc_,zc_,x_[j],y_[j],z_[j],c[j]);
                   }                         
                         
                         
                   break;                   
    case   RBF_DXX: 
                   for( i=ini;  i<fin; i++ ) 
                   {
                       xc_ = x_[i];
                       yc_ = y_[i];
                       zc_ = z_[i];  
                       
                       for( j=0;  j<n;  j++ ) 
                          A(i,j) = rbf.dxx(xc_,yc_,zc_,x_[j],y_[j],z_[j],c[j]);
                   }                         
                        
                   break;
                   
    case   RBF_DYY: 
                   for( i=ini;  i<fin; i++ ) 
                   {
                       xc_ = x_[i];
                       yc_ = y_[i];
                       zc_ = z_[i];  
                       
                       for( j=0;  j<n;  j++ ) 
                         A(i,j) = rbf.dyy(xc_,yc_,zc_,x_[j],y_[j],z_[j],c[j]);
                   }                         
                        
                   break;
                   
                   
    case   RBF_DZZ: 
                   for( i=ini;  i<fin; i++ ) 
                   {
                       xc_ = x_[i];
                       yc_ = y_[i];
                       zc_ = z_[i];  
                       
                       for( j=0;  j<n;  j++ ) 
                          A(i,j) = rbf.dzz(xc_,yc_,zc_,x_[j],y_[j],z_[j],c[j]);
                   }                         
                        
                   break;
                   
     default:
                 fprintf(stderr,"\nERROR: undefined option '%s' in function fillMatrix\n",type.c_str());
                 fprintf(stderr,"in file '%s' near to line %d\n",__FILE__,__LINE__);
                 exit(1);
                 break;               
        
  } //switch
  
//now, insert the polynomial part


if(dim_px>0)
{
  Mat P;
  Vec X;
  switch(select(type))
  {
     case   POL_TRANSPUESTO:
                       //now evaluate  the polynomial   
                         compose(x,y,z,X);
                         
                         P = pol.build_tnt(X.GetData(),n,dim);

                        // P =  transpose(P);

                         for(int i=ini; i<fin; i++) 
                             for(int j=0; j<n; j++)
                                A(i,j) = P(j,i-ini);
                               // A(i,j) = P(i-ini,j);
          
                     //restore the info of the original polynomial 
                         pol.make(dim, pol.get_m());
                         return;
                  break; 
        
    case   RBF_NORMAL:
                   break;
        
    case   RBF_DX:  pol.deriva("x",1);  
                   break; 
        
    case   RBF_DY: pol.deriva("y",1);
                   break; 
                   
    case   RBF_DZ: pol.deriva("z",1);
                   break; 
        
    case   RBF_DXX: pol.deriva("x",2);
                   break;    
        
    case   RBF_DYY: pol.deriva("y",2);
                   break;  

    case   RBF_DZZ: pol.deriva("z",2);
                   break;  
        
    default:   fprintf(stderr,"\nERROR: undefined or invalid option '%s' in function fillMatrix\n",type.c_str());
               fprintf(stderr,"in file '%s' near to line %d\n",__FILE__,__LINE__);
               exit(1);
              break;
        
  } //switch
  


       //now evaluate the derivative of the polynomial
         for(int i=ini; i<fin; i++) 
         {               
              xc[0] = x(i);
              xc[1] = y(i);
              xc[2] = z(i);
              
              pol.eval(xc,dim, &A(i,n), dim_px );
          }   
          
          
      //restore the info of the original polynomial 
         pol.make(dim , pol.get_m());


 } //if(rbf.has_pol()) 
 

} 

//----------------------------------------------------------
// Data 3D
//---------------------------------------------------------- 
template <typename RBF, typename Vec, typename T, typename Pol, typename Mat>
void fill_matrix(string type,               // select the operation to perform
               RBF rbf, Pol &pol,           // defined the configuration of RBF + Pol 
               T  c, 
               Vec &x, Vec &y, Vec& z,     // three input vectors
               int ini, int fin,           // from row=ini to row<fin in A 
               Mat &A                      // output
               )
{
   Vec myc(x.GetSize());
 
   myc = c;
 
   fill_matrix(rbf,pol,myc,x,y,z,ini,fin,A);
}      


} // RBF namespace

#endif //  _FILL_3D_HPP_



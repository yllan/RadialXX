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


#ifndef _FILL_1D_HPP_
#define _FILL_1D_HPP_


namespace rbf{


//----------------------------------------------------------
// Data 1D
//---------------------------------------------------------- 
template <typename RBF, typename Vec, typename T, typename Pol, typename Mat>
void fill_matrix(string type, 
                RBF rbf, Pol &pol, 
                T  c,
                Vec &x,  
                int ini, int fin, 
                Mat &A
               )
{
   int  n,dim_px;
   int dim;
   int  i,j;
   T    xc[1];
   T*   x_ = x.GetData();  
   T    xc_;
  
    
   assert( ini<=fin );
  
  
//stablish the data dimension
   dim = 1;
   
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
                            
                          for( j=0;  j<n;  j++ ) 
                             A(i,j) = rbf.eval(xc_,x_[j],c);
                             
                       }      
                   break;           
    case   POL_TRANSPUESTO:
                  break;
        
    case   RBF_DX: 
                       for( i=ini;  i<fin;  i++ ) 
                       {
                            xc_ = x_[i];
                            
                          for( j=0;  j<n;  j++ ) 
                             A(i,j) = rbf.dx(xc_,x_[j],c);
                             
                       }    
                   break;
        
    case   RBF_DXX:  
                       for( i=ini;  i<fin;  i++ ) 
                       {
                            xc_ = x_[i];
                            
                          for( j=0;  j<n;  j++ ) 
                             A(i,j) = rbf.dxx(xc_,x_[j],c);
                             
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
  switch(select(type))
  {
     case   POL_TRANSPUESTO:
                       //now evaluate  the polynomial   
                         P = pol.build_tnt(x.GetData(),n,dim);

                         //P =  transpose(P);

                         for( i=ini;  i<fin;  i++ ) 
                             for( j=0;  j<n;  j++ )
                                A(i,j)  = P(j,i-ini);
          
                     //restore the info of the original polynomial 
                         pol.make(dim , pol.get_m());
     
                         return;
                  break; 
        
    case   RBF_NORMAL:
                   break;
        
    case   RBF_DX: pol.deriva("x",1);
                   break; 
        
    case   RBF_DXX: pol.deriva("x",2);
                   break;    
        
    default:   fprintf(stderr,"\nERROR: undefined or invalid option '%s' in function fillMatrix\n",type.c_str());
               fprintf(stderr,"in file '%s' near to line %d\n",__FILE__,__LINE__);
               exit(1);
              break;
        
  } //switch
  
       //now evaluate the derivative of the polynomial
         for( i=ini;  i<fin;  i++ ) 
         {               
              xc[0] = x(i);
              
              pol.eval(xc,dim, &A(i,n), dim_px );
          }   
          
      //restore the info of the original polynomial 
         pol.make(dim , pol.get_m());

 } //if(rbf.has_pol()) 
 
} 



//----------------------------------------------------------
// Data 1D
//---------------------------------------------------------- 
template <class RBF, typename T>
void fill_matrix(string type, 
                RBF rbf,   Polinomio<T> &pol, 
                Vector<T> &c,
                Vector<T> &x,  
                int ini, int fin, 
                Matrix<T> &A)

{
  int  n,dim_px;
  int  dim;
  int  i,j;
  T    xc[1];
  T*   x_ = x.GetData();  
  T*   c_ = c.GetData();  
  T    xc_;
  
    
  assert( ini<=fin );
  
//stablish the data dimension
   dim = 1;
  
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
                            
                          for( j=0;  j<n;  j++ ) 
                             A(i,j) = rbf.eval(xc_,x_[j],c_[j]);
                             
                       }      
                   break;           
    case   POL_TRANSPUESTO:
                  break;
        
    case   RBF_DX: 
                       for( i=ini;  i<fin;  i++ ) 
                       {
                            xc_ = x_[i];
                            
                          for( j=0;  j<n;  j++ ) 
                             A(i,j) = rbf.dx(xc_,x_[j],c_[j]);
                             
                       }    
                   break;
        
    case   RBF_DXX:  
                       for( i=ini;  i<fin;  i++ ) 
                       {
                            xc_ = x_[i];
                            
                          for( j=0;  j<n;  j++ ) 
                             A(i,j) = rbf.dxx(xc_,x_[j],c_[j]);
                             
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
  Matrix<T> P;
  switch(select(type))
  {
     case   POL_TRANSPUESTO:
                       //now evaluate  the polynomial   
                         P = pol.build_tnt(x.GetData(),n,dim);

                         //P =  transpose(P);

                         for( i=ini;  i<fin;  i++ ) 
                             for( j=0;  j<n;  j++ )
                                A(i,j)  = P(j,i-ini);
          
                     //restore the info of the original polynomial 
                         pol.make(dim , pol.get_m());
     
                         return;
                  break; 
        
    case   RBF_NORMAL:
                   break;
        
    case   RBF_DX: pol.deriva("x",1);
                   break; 
        
    case   RBF_DXX: pol.deriva("x",2);
                   break;    
        
    default:   fprintf(stderr,"\nERROR: undefined or invalid option '%s' in function fillMatrix\n",type.c_str());
               fprintf(stderr,"in file '%s' near to line %d\n",__FILE__,__LINE__);
               exit(1);
              break;
        
  } //switch
  
       //now evaluate the derivative of the polynomial
         for( i=ini;  i<fin;  i++ ) 
         {               
              xc[0] = x(i);
              
              pol.eval(xc, dim, &A(i,n), dim_px );
          }   
          
      //restore the info of the original polynomial 
         pol.make(dim , pol.get_m());

 } //if(rbf.has_pol()) 
 
} 



} // RBF namespace

#endif //  _FILL_2D_HPP_


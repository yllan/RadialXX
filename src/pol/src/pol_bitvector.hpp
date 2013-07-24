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


#ifndef _POL_BIT_VECTOR_HPP_
#define _POL_BIT_VECTOR_HPP_

#include <cmath>

class TVectorBits
{
protected:
   int*   data;
   int       n;
   int    base;	
   int    size;
   int    cont;	
  
private:
   int add(int &acum,int x, int y);

public:
   TVectorBits(void);
  ~TVectorBits(void);
   void  initialize(int n, int base);
   void  first(int *x);
   bool  next(int *x);
   int   get_size(void);
};
//----------------------------------------------------------
TVectorBits::TVectorBits(void)
{
  cont = 0;
  size = 0;
     n = 0;
 data  = NULL;
}
//----------------------------------------------------------
TVectorBits::~TVectorBits(void)
{
 if(n>0)	
   delete[] data;

  cont = 0;
  size = 0;
     n = 0;
}
//----------------------------------------------------------
int TVectorBits::get_size(void)
{
   return size;	
}
//----------------------------------------------------------
void TVectorBits::initialize(int n, int base)
{
   int    max_int;
   double size_tmp;
 
   if(n<=0)
   {
      fprintf(stderr,"\nERROR in TVectorBits.initialize, n = %d must be a positive value.\n\n",n);
      exit(1);	
   }
   
   if(n==1)
   {
      fprintf(stderr,"\nERROR in TVectorBits.initialize, n = %d must be greather than one.\n\n",n);
      exit(1);	
   }
   
   if(base<=0)
   {
      fprintf(stderr,"\nERROR in TVectorBits.initialize, base = %d must be a positive value.\n\n",n);
      exit(1);	
   }
   
   if(base==1)
   {
     // fprintf(stderr,"\nERROR in TVectorBits.initialize, beta = %d must be greather than one.\n\n",n);
     // exit(1);	
   }
 
   max_int  = numeric_limits<int>::max();
  
   size_tmp = pow((double)base,(double)n);
 
   if(size_tmp >= max_int )	
   {
 	
      fprintf(stderr,"\nERROR: in TVectorBits.initialize  with n=%d and base=%d the size of all combinations ...\n",n,base);
      fprintf(stderr,"       are (%f) which is upper of the max_int=%d\n\n",size_tmp,max_int);	
      exit(1);
   }	
   
   if(size_tmp<=0)
   {
      fprintf(stderr,"\nERROR: in TVectorBits.initialize the pow(base,n) --> %f \n\n",size_tmp);
      exit(1);	
   }
   
   
   this->base = base;	
	
   cont = 0;
  
   size = int(size_tmp);

   if(data!=NULL)
      delete[] data;

   data = new int[n];
	
   if(data==NULL)
   {
  	  fprintf(stderr,"\nERROR: in TVectorBits::initalize the dynamic memory can not be created.\n\n");
	  exit(1);
   }	 

   this->n = n;

}
//----------------------------------------------------------
void TVectorBits::first(int *x)
{
	for( int j=0; j<n; j++)
	{
       data[j] = 0;
	  	  x[j]  = 0;	
	}	
	cont = 0;
    cont = cont+1;
}
//----------------------------------------------------------
bool  TVectorBits::next(int *x)
{
   int acum;

 //add one to the vector with the bits
   acum=0;

   data[n-1] = add(acum,data[n-1],1);
	
   for( int j=(n-2); j>=0; j--)
       data[j] = add(acum,data[j],0);

   for( int j=0; j<n; j++)	
	   x[j] = data[j];
 
   cont=cont+1;

	
   if(cont>size)
	   return false;
   else
      return true;
}
//----------------------------------------------------------
int TVectorBits::add(int &acum,int x, int y)
{
   int z,out=0;	
   
   switch(acum)
   {
	case 0:
	       z = x + y;
		if(z<base){
			acum = 0;
			//return z;
			out = z;
		}else{
			acum = 1;
			//return 0;
			out = 0;
		}
		break;
	case 1:
		z = x + y + acum;
		if(z<base){
			acum = 0;
			//return z;
			out = z;
		}
		else if(z==base){
			acum = 1;
			//return 0;
			out = 0;
		}
		else if(z==(base+1)){
			acum = 1;
			//return 1;
			out = 1;
		}
		else{
			fprintf(stderr,"ERROR: int TVectorBits:add  acum = %d   x = %d   y= %d\n\n",acum,x,y); 
			exit(1);
		}
		
		break;
   } // switch
   
 return out;
 
}
#endif  //_POL_BIT_VECTOR_HPP_


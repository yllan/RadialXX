/***************************************************************************
 *
 * Copyright (C) 2009   Jose Antonio Munoz Gomez
 *        
 * This file is part of Radial++
 * http://sourceforge.net/projects/radial/
 *
 * Radial++ is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 3
 * of the License, or (at your option) any later version.
 *
 * Radial++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * For more information, please see the  Home Page:
 *    http://www.dci.dgsca.unam.mx/pderbf
 *
 ***************************************************************************/


#include "radial.h"

  
//---------------------------------------------------------
int main(int argc, char **argv)
{
   TVectorBits vbits;
   int  *vec;
   int   d,m;
   int   sum,cont;
   
   d = 10;
   m = 16;

           cont = 0;
          //create the dynamic memory
            vec = new int[d]; 
       
          //create the base vector of longitud d in the base number m
            vbits.initialize(d,m);
            
            printf("size vbits = %d cont_pol = %d\n",vbits.get_size(),vbits.cont);
            
            vbits.first(vec);
        
          do{
            sum = 0;
            for(int i=0; i<d; i++)  
              sum  +=  vec[i];
            
            if(sum<m){
            	cont++;
            	 printf("cont=%d size=%d , cont_pol=%d: ",cont,vbits.get_size(),vbits.cont);
                  for(int i=0; i<d; i++)  
                    printf("%d ",vec[i]);
                  printf("\n");  
               	}


               
         }while(vbits.next(vec));




  return 0;
}

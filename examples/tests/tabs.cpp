#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
  FILE *in,*out;
  char file_out[512];
  char  c;
  int   cont_line=1;


 if(argc!=2)
 {
    printf("ERROR:  tabs file.cpp  output->file.cpp.wt\n\n");
    return 1;   
  }
  in = fopen(argv[1],"r");
  
  sprintf(file_out,"%s.wt",argv[1]);
  
  out = fopen(file_out,"w");
  

  if(in==NULL)
  {
    fprintf(stderr,"ERROR: unabled to open '%s'\n\n",argv[1]);
    return 1;
  }

  if(out==NULL)
  {
    fprintf(stderr,"ERROR: unabled to open '%s'\n\n",file_out);
    return 1;
  }


  while(!feof(in))
  {
    c = fgetc(in);
    if(c=='\n')
    {
      cont_line++;
      fputc(c,out);
    }
    else if(c=='\t')
    {
      printf("WARNING:  detected tab near to line '%d' \n",cont_line);
       fputc(' ',out);
       fputc(' ',out);
       fputc(' ',out);
    }
    else{
    if(!feof(in))   
        fputc(c,out);
    }   

  }


  fclose(in);
  fclose(out);
  
 
 
 
   in = fopen(file_out,"r");
  out = fopen(argv[1],"w");
  
  
  while(!feof(in))
  {
    c = fgetc(in);
    if(!feof(in))   
        fputc(c,out);

  }


  fclose(in);
  fclose(out);
  
  
  system("rm -f *.wt");
  
  

 
  return 0;
}

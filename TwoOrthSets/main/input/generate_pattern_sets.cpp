#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <cstdlib>  // For srand() and rand()
//#include <mpi.h>
/*----------------------------------------------------------------------*
*		Generate Random Patterns written by Vizhe						*		
* ----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
*																		*
*							Enter Parameters							*
* 																		*		
* ----------------------------------------------------------------------*/

/*network parameters*/
int N = 5000;
int p = 60;
float a = 0.01;

FILE *fpattern;

typedef struct {
int index;
double field;
} mystruct;
  

int cmpf (const void *x, const void *y)
{
  double xx = (*(mystruct*)x).field;
  double yy = (*(mystruct*)y).field;
  /*sort by decreasing order*/
  if (xx > yy) return -1;
  if (xx < yy) return  1;
  return 0;
}

int main()
{
  srand48(time(0));	
 
  int **xi;

  xi= new int*[p];
  for(int i=0; i<p; i++)
    xi[i]=new int[N];
    
  char *buffer; 
  buffer = new char [1000*sizeof(char)];
  sprintf(buffer, "pattern_set0_N%d_a%.2f_p%d", N, a, p);
  fpattern = fopen(buffer, "w");    


  /*first initialize to zero*/
  for(int mu=0; mu<p; mu++)
  {
	  for(int i=0; i<N; i++)
	  {
	  	xi[mu][i] = 0;
	  }	
  }

  /* FIRST SET*/
  for(int mu=0; mu<int(p/2); mu++)
  {
	mystruct X[N];
    for(int i=0; i<int(N/2); i++)
    {
      X[i].index = i;
      X[i].field = drand48();
    }
    for(int i=int(N/2); i<N; i++)
    {    
      X[i].index = i;
      X[i].field = 0;
    }    	
    /*sort the patterns in order of decreasing field keeping the indices*/
    qsort(X, N, sizeof(mystruct), cmpf);  

    /*set N*a largest of them to become active*/
    for(int i=0; i<int(N*a); i++)
    {
	  xi[mu][X[i].index] = 1; 
    }  
  }

  /* SECOND SET*/

  for(int mu=int(p/2); mu<p; mu++)
  {
	mystruct X[N];
    for(int i=0; i<int(N/2); i++)
    {    
      X[i].index = i;
      X[i].field = 0;
    }    	
    for(int i=int(N/2); i<N; i++)
    {
      X[i].index = i;
      X[i].field = drand48();
    }
    /*sort the patterns in order of decreasing field keeping the indices*/
    qsort(X, N, sizeof(mystruct), cmpf);  

    /*set N*a largest of them to become active*/
    for(int i=0; i<int(N*a); i++)
    {
	  xi[mu][X[i].index] = 1; 
    }  
   }
  
  /*compressed writing, not ones and zeros, but only indices of non-zeros*/
  for(int mu=0; mu<p; mu++)
  {    
    for(int i=0; i<N; i++)
    {
    if(xi[mu][i] == 1)
    { 
      fprintf(fpattern, "%d ", i);
    }
   }
  fprintf(fpattern, "\n");
  }
}

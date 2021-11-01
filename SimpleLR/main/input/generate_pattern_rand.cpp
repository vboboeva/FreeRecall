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
int p = 1000;
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
  int i, mu;
  
  int **xi;

  xi= new int*[p];
  for(i=0; i<p; i++)
    xi[i]=new int[N];
    
  char *buffer; 
  buffer = new char [1000*sizeof(char)];
  sprintf(buffer, "pattern_N%d_a%.2f", N, a);
  fpattern = fopen(buffer, "w");    
  
  for(mu=0;mu<p;mu++)
  {
    mystruct X[N];
 
    for(i=0;i<N;i++)
    {
	  
      /*first initialize each pattern to null state*/
      xi[mu][i] = 0;		
      X[i].index = i;
      X[i].field = drand48();
    }
    /*sort the patterns in order of decreasing field keeping the indices*/
    qsort(X, N, sizeof(mystruct), cmpf);  

    /*set N*a largest of them to become active*/
    for(i=0;i<N*a;i++)
    {
	  xi[mu][X[i].index] = 1; 
    }  
  
    /*standard writing*/
    /*for(i=0;i<N;i++)
    {
      fprintf(fpattern, "%d ", xi[mu][i]);
    }
    fprintf(fpattern, "\n");*/
    /*compressed writing, not ones and zeros, but only indices of non-zeros*/
    for(i=0;i<N;i++)
    {
      if(xi[mu][i] == 1) 
      fprintf(fpattern, "%d ", i);
    }
    fprintf(fpattern, "\n");
  }
}

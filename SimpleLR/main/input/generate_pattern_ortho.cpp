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

  /*network parameters*/
  int N = 5000;
  int p = 100;
  float a = 0.02;  

  if (1./a < 1.*p) 
    {
      printf("not enough units! Making minimum number possible %d\n", int(1./a));
      p=int(1./a);
    }
  
  srand48(time(0));	
  int i, mu;
  
  int **xi;

  xi= new int*[p];
  for(i=0; i<p; i++)
    xi[i]=new int[N];
    
  char *buffer; 
  buffer = new char [1000*sizeof(char)];
  sprintf(buffer, "pattern_N%d_a%.2f_ortho", N, a);
  fpattern = fopen(buffer, "w");    
  
  for(mu=0;mu<p;mu++)
  {
    for(i=0;i<N*a;i++)
    {
      fprintf(fpattern, "%d ", int(mu*N*a)+i);
    }
    fprintf(fpattern, "\n");
  }
}


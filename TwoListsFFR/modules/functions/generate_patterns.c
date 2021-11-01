#define GENERATE_PATTERNS_C

#include <stdio.h> 
#include <stdlib.h>
#include <math.h> 
#include "random.h"
#include "global.h"
#include "const_potts.h"

/*----------------------------------------------------------------------*
*		New pattern generator written by Vizhe, version 4,				*
* 		that defines factors as random patterns acting on random		*
* 		subsets of patterns	and not units, as does Ale (version 0).		*
* 		The objective is to reverse the									*
* 		correlation statistics of Ale's algorithm, that is,				*
* 		skewed correlations distribution between patterns, 			 	*
* 		and symmetric between units. Difference with version 1: 		*
* 		parents are random patterns that send influence in a 			*
* 		different Potts state for each unit of a pattern.				*
* 		Difference with version 3: random input included to solve		*
* 		small sparsity for small values of a_pf.						*	
*		Limit cases:													*
* 		1) a_pf ---> 0 gives us randomly correlated patterns			*
* 		2) f*Num_fact = 1 gives us ultrametric patterns	 (f=p_fact/p)	*
* ----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
*																		*
*							Enter Parameters							*
* 																		*		
* ----------------------------------------------------------------------*/

/*No need, done in global.h! */


/*----------------------------------------------------------------------*
*																		*
*							FUNCTIONS									*
* 																		*		
* ----------------------------------------------------------------------*/

void GetMemory()
{
int i;

Factors= malloc(Num_fact*sizeof(int*));
for(i=0; i<Num_fact; i++)
	Factors[i]=malloc(N*sizeof(int));
	
PC= malloc(Num_fact*sizeof(int*));
for(i=0; i<Num_fact; i++)
	PC[i]=malloc(p_fact*sizeof(int));
 
hpat= malloc(N*sizeof(double*));
for(i=0; i<N; i++)
	hpat[i]=malloc((S+1)*sizeof(double));

hmax = malloc(N*sizeof(double));
smax = malloc(N*sizeof(int));

}

typedef struct{
int index;
double field;
int state;
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

/*Factors are just random patterns*/
void SetFactors()
{
  mystruct F[N];
  
  int n, i;
  for(n=0;n<Num_fact;n++)
  {
    for(i=0;i<N;i++)
    {
      /*first initialize each pattern to null state*/
      Factors[n][i] = S;		
      
      F[i].index = i;
      F[i].field = drand48();
      F[i].state = (int)((double)(S)*drand48());
    }
    /*sort the units in order of decreasing field keeping the indices: this is done to make patterns with fixed sparsity*/
    qsort(F, N, sizeof(mystruct), cmpf);  

    /*set N*a largest of them to become active*/
    for(i=0;i<N;i++)
    {
	  Factors[n][F[i].index] = F[i].state; 
    }  
  }
}

void AssignChildren()
{
  int n, m, i, patt_picked;
  int rand_num;

  for(n=0; n<Num_fact; n++)
  {
	m=0;
	while(m<p_fact)
	{
      rand_num = (int)((double)p*drand48());
      patt_picked = 0;
      /*check that it hasn't been assigned before*/
      for(i=0; i<m; i++)
      {
		if(PC[n][i] == rand_num) patt_picked = 1;
	  }
	  if (patt_picked == 0) 
	  {
		PC[n][m] = rand_num;
		m++;
	  }
    }
  }
   
}

void SumFields(int mu)
{

int i, n, m, k, rand_state;
double y, eigen_fact, expon;


/*set fields to 0*/
for(i=0;i<N;i++)
{
	for(k=0;k<S+1;k++)
	{
		hpat[i][k] = 0.0;
	}
}

/*each pattern sums field coming from all its parents*/
for(n=0; n<Num_fact; n++)
{
    expon = -fact_eigen_slope*n;
	for(m=0;m<p_fact; m++)
	{
		if(PC[n][m] == mu)
		{
			for(i=0;i<N;i++)
			{
			  /*component coming from parents*/	
		      y = (double)drand48();
		      if(y<=a_pf)
		      {	
		        eigen_fact = exp(expon)*y/a_pf;
				hpat[i][Factors[n][i]] += eigen_fact;
				/*printf("%.2f \t", hpat[i][Factors[n][i]]);*/
		      }
		    }
		}
	}
}
/*for children that have no parents, or when a_p is too small (this is epsilon in the paper + thesis)*/
for(i=0;i<N;i++)
{
  /*small input to a random state*/
  rand_state = (int)((double)(S)*drand48());
  hpat[i][rand_state] += eps*drand48();
}

/*find state of maximal field*/
double hM;
int SM;
for(i=0;i<N;i++)
{
	hM = 0.0;
	SM = S;
	for(k=0;k<S+1;k++)
	{
		/*printf("%d \t %d \t %.2f \n", i, k, hpat[i][k]);*/
		if(hpat[i][k] > hM) 
		{
			hM = hpat[i][k];
			SM = k;
		}
	}
	hmax[i] = hM;
	smax[i] = SM;
	/*printf("%.2f\t", hmax[i]);*/
	/*printf("%d\n", smax[i]);*/
}

}

void SetPatterns()
{
  int i, mu, count;
  for(mu=0;mu<p;mu++)
  {
	/*for each pattern sum the fields*/  
    SumFields(mu);
    
    mystruct X[N];
 
    for(i=0;i<N;i++)
    {
      /*first initialize each pattern to null state*/
      xi[mu][i] = S;		
      
      X[i].index = i;
      X[i].field = hmax[i];
      X[i].state = smax[i];
    }
    /*sort the patterns in order of decreasing field keeping the indices*/
    qsort(X, N, sizeof(mystruct), cmpf);  

    /*set N*a largest of them to become active*/
   
    i=0;
    count=0;
    while (count<N*a && i<N)
    {
	  if (X[i].state != S) 
	  {
		xi[mu][X[i].index] = X[i].state;
		count++;
      }
	i++;
    }
  }
}

void SavePatterns()
{
int mu, i;
	
char buffer[0x100];
snprintf(buffer, sizeof(buffer), "pattern_S%d_a%.2f_apf%.2f_pfact%d_p%d_Numfact%d", S, a, a_pf, p_fact, p, Num_fact);
fpattern = fopen(buffer, "w");	
	
for(mu=0;mu<p;mu++) 
{
    for(i=0;i<N;i++)
    {
      fprintf(fpattern, "%d ", xi[mu][i]);
    }
  fprintf(fpattern, "\n");
}
fclose(fpattern);
}


#define FUNCTIONS_C

#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <float.h> 
#include "random.h"
#include "params.h"
#include "global.h"
#include <limits.h>
/*--------------------------------------------*/

typedef struct {
int index;
double field;
} mystruct;

/*--------------------------------------------*/

void print_m(double t) 
{
	fprintf(mall, "%.2f\t", t);

	for(int mu=0;mu<p;mu++)
	{
	  fprintf(mall, "%.4f\t", m[mu]);
	}

	fprintf(mall, "\n");
	fflush(mall);
}

/*--------------------------------------------*/

int cmpf (const void *x, const void *y)
{
  double xx = (*(mystruct*)x).field;
  double yy = (*(mystruct*)y).field;
  /*sort by decreasing order*/
  if (xx > yy) return -1;
  if (xx < yy) return  1;
  return 0;
}

/*--------------------------------------------*/


void read_p_patterns(int lineNum) 
{
	int index;

    /*move through lines without scanning*/
	for (int j=0; j < lineNum; j++)
    {
	    for(int k=0;k<N*a;k++)
	    {    	
          fscanf(pat, "%*d ");
        }
    }

	for(int mu=0;mu<p;mu++)
	{
	  for(int i=0;i<N;i++)
	  {
	  	xi[mu*N+i]=0;
	  }
	  for(int i=0;i<N*a;i++)
	  {
	    fscanf(pat, "%d", &index);
	    /*printf("%d\t",index);*/
	    xi[mu*N+index]=1;
	  }
	  /*printf("\n");*/
	}
}


/*--------------------------------------------*/

void pick_Lrandom_patterns() 
{
	int *elements = malloc(p*sizeof(int));
	double uu[2*L];
	ranlxd(uu,2*L);

	/* List 2 */

	for (int i = 0; i < p; ++i)
		elements[i] = i;

	int temp, w;
	for (int i = 0; i < 2*L; i++)
	{
		w = i+ (int)(INT_MAX*uu[i])%(p-i);
		list_idx[i] = elements[w];
		temp = elements[i];
		elements[i] = elements[w];
		elements[w] = temp;
	}

	for (int i = 0; i < L; i++)
	{
		printf("%d ", list_idx[i]);
		fprintf(list_indices, "%d\t", list_idx[i]);
	}
	printf("\n");	
	fprintf(list_indices,"\n");

	for (int i = L; i < 2*L; i++)
	{
		printf("%d ", list_idx[i]);
		fprintf(list_indices, "%d\t", list_idx[i]);
	}
}

/*--------------------------------------------*/

void initialize_net() 
{	
	for (int i=0;i<N;i++)
	{
		u[i]=U;/*U+u[i]*(1.-U);*/
		x[i]=0.0;
		h[i]=0.0;
		r[i]=0.0;
		s[i]=0.0;
		g[i]=g0;

		for(int j=0;j<N;j++)
		{
			JSTP[i*N+j]=0.;
		}
	}
}

/*--------------------------------------------*/

void reinitialize_net() 
{	
	for (int i=0;i<N;i++)
	{
		u[i]=U;/*U+u[i]*(1.-U);*/
		x[i]=0.0;
		h[i]=0.0;
		r[i]=0.0;
		s[i]=0.0;
		g[i]=g0;
	}
}

/*--------------------------------------------*/


void construct_J()
{
	/*MAKE J*/
	for(int i=0; i<N; i++)
	{
		for(int j=0; j<N; j++)
		{
			J[i*N+j]=0.0;

			/*set self-connections to zero*/
			if (i != j)
			{
				for(int mu=0; mu<p; mu++)
				{
					J[i*N+j] += (1.*xi[mu*N+i]-a)*(1.*xi[mu*N+j]-a);
				}
				J[i*N+j] /= N*a*(1.-a);
			}
		/*printf("%f\t",J[i][j]);*/
	    }
	}
	printf("after J_ij \n");
}

/*--------------------------------------------*/


void update_state(int n, int i, double extfield, double backfield, int stoplearning)
{
	double y_EE, y_EI, z;
	double rand_num[1];
	y_EE=0.0;
	y_EI=0.0;
	
	/*first do the sums*/
	for(int j=0;j<N;j++)
	{	
		y_EE += (J[i*N+j] + JSTP[i*N+j])*s[j]*u[j]; 
	}
	for(int j=0;j<N;j++)
	{
		y_EI += s[j];
	}
	y_EI /= N*a;
	/*then update states*/
	/*extfield is specific to the cue*/
	/*readoutfield is nonspecific*/

	ranlxd(rand_num, 1);

	h[i] = y_EE - kappa*y_EI*y_EI + eps*(2*rand_num[0]-1.) + backfield + extfield; 
	r[i] += dt*(-r[i] + h[i] - x[i])/tau; 
	x[i] += dt*(s[i] - x[i])/tau_A;
	u[i] += dt*(U*(1./tau_F + s[i]) - (1./tau_F + U*s[i])*u[i]); 


	double denom=1./(N*a*a);

	for(int j=0;j<N;j++)
	{
		if(stoplearning==0 && JSTP[i*N+j]<= 0.01) 
		
			JSTP[i*N+j] += denom*dt*( alpha*(s[i]*s[j] - s[i]*s[i]*JSTP[i*N+j])/(1.*L) 
						+ beta*(s[i]*g[j] - s[i]*s[i]*JSTP[i*N+j])) - dt*zeta*JSTP[i*N+j] ;
		if(stoplearning==0 && JSTP[j*N+i]<= 0.01)
		{
			JSTP[j*N+i] += denom*dt*( gamma*(s[i]*g[j] - g[j]*g[j]*JSTP[j*N+i])) - dt*zeta*JSTP[j*N+i]; 
		}
	}
	/*set self-connections to zero*/

	JSTP[i*N+i]=0.;

	/*s[i] = 1./(1.+exp(-beta*(r[i])));*/

	/*Numerically stable sigmoid function*/
	if(r[i] >= 0)
	{
		z = exp(-invTemp*(r[i]-Th));
	    s[i] = 1./(1. + z);
	}
	else
	{
		z = exp(invTemp*(r[i]-Th));
	    s[i] = z/(1. + z);
	}

	g[i] += dt*(s[i]-g[i])/tau_pre;

}
/*--------------------------------------------*/
void update_aux()
{
	for (int i = 0; i < N; i++)
	{
		g[i] += dt*(s[i]-g[i])/tau_pre;
	}
}
/*--------------------------------------------*/
int* findmax(int* temp)
{
	float Mmax=-1.;
	float Mmax2=-2.;
	int Mumax=p+1;
	int Mumax2=p+2;

	for(int mu=0;mu<p;mu++) /*find pattern that has maximal overlap with net*/
	{
		if(m[mu]>Mmax)
		{
			Mmax2=Mmax;
			Mmax=m[mu];
			Mumax2=Mumax;
			Mumax=mu;
		} 
		else if(m[mu]>Mmax2)
		{
			Mmax2=m[mu];
			Mumax2=mu;
		}
	}

	temp[0]=Mumax;
	temp[1]=Mumax2;
	return(temp);
}
/*--------------------------------------------*/

void compute_m()
{
	double ma;

	for(int mu=0;mu<p;mu++)
	{
	  ma=0.;
	  for(int i=0;i<N;i++)
	  {
	    ma+=(1.*xi[mu*N+i]-a)*s[i];  
	  }
	m[mu]=ma/(1.*N*a*(1.-a));
	}
}

/*--------------------------------------------*/

void compute_s_h_u_x(double t)
{
	fprintf(sall, "%.2f\t", t);
	fprintf(hall, "%.2f\t", t);
	fprintf(uall, "%.2f\t", t);
	fprintf(xall, "%.2f\t", t);
	fprintf(Jall, "%.2f\t", t);
	for(int mu=0;mu<p;mu++)
	{
		smean[mu]=0.;
		hmean[mu]=0.;
		umean[mu]=0.;
		xmean[mu]=0.;
		Jmean[mu]=0.;
		for(int i=0;i<N;i++)
		{
			smean[mu]+=s[i];  
			if (xi[mu*N+i]==1)
			{
				hmean[mu]+=h[i];  
				umean[mu]+=u[i];  
				xmean[mu]+=x[i];  
				
				for(int j=0;j<N;j++)
				{
					if (xi[mu*N+j]==1 && j!=i)
					{
						Jmean[mu]+=JSTP[i*N+j];
					}			
				}				
			}
			
		}
		smean[mu]/=(1.*N);
		hmean[mu]/=(1.*N*a);
		umean[mu]/=(1.*N*a);
		xmean[mu]/=(1.*N*a);
		Jmean[mu]/=(1.*N*a*N*a);

		fprintf(sall, "%.4f\t", smean[mu]);
		fprintf(hall, "%.4f\t", hmean[mu]);
		fprintf(uall, "%.4f\t", umean[mu]);	
		fprintf(xall, "%.4f\t", xmean[mu]);	
		fprintf(Jall, "%.4f\t", Jmean[mu]);	
	}
	fprintf(sall, "\n");
	fprintf(hall, "\n");
	fprintf(uall, "\n");
	fprintf(xall, "\n");
	fprintf(Jall, "\n");
}

/*--------------------------------------------*/

void compute_u_J(double t)
{
	int mu;
	for(int nu=0;nu<L;nu++)
	{
		mu=list_idx[nu];
		umean[mu]=0.;
		Jmean[mu]=0.;
		for(int i=0;i<N;i++)
		{
			if (xi[mu*N+i]==1)
			{
				umean[mu]+=u[i];  
				
				for(int j=0;j<N;j++)
				{
					if (xi[mu*N+j]==1 && j!=i)
					{
						Jmean[mu]+=J[i*N+j]+JSTP[i*N+j];
					}			
				}				
			}
			
		}
		umean[mu]/=(1.*N*a);
		Jmean[mu]/=(1.*N*a*N*a);
		fprintf(UvsJ, "%f\t",umean[mu]);
	}
	for(int nu=0;nu<L;nu++)
	{
		mu=list_idx[nu];
		fprintf(UvsJ, "%f\t",Jmean[mu]);
	}	
}
/*--------------------------------------------*/

void print_Jij()
{
	/*now print to file*/
	for (int i=0; i<N; i++) 
	{ 
		for(int j=0; j<N; j++) 
		{
			fprintf(Jij, "%.2f ", J[i*N+j]+JSTP[i*N+j]);
		}
	fprintf(Jij, "\n");
	}
}
/*--------------------------------------------*/

void print_Jmean_before_after()
{
	for(int mu=0;mu<p;mu++)
	{
		Jmean[mu]=0.;
		for(int i=0;i<N;i++)
		{
			if (xi[mu*N+i]==1)
			{
				for(int j=0;j<N;j++)
				{
					if (xi[mu*N+j]==1)
					{
						Jmean[mu]+=J[i*N+j]+JSTP[i*N+j];
					}			
				}				
			}
		}
		Jmean[mu]/=(1.*N*a*N*a);
		fprintf(Jm, "%f\t",Jmean[mu]);
	}
	fprintf(Jm, "\n");
}

/*--------------------------------------------*/

/*produces NumSet sequences of random integers from 1 to 10 withOUT replacement.
this is the order of updating of units.
the idea is to be able to use these sequences randomly over runs,
but still be able to use a single sequence if needed.*/

void SetUpTables() 
{
	int item, randunit;
	int done;  
	double rand_num[1];

	for(int k=0; k<NumSet; k++)
	{
	  item = 0; 
	  while(item<N)
	  {
	    ranlxd(rand_num, 1);
	    randunit = (int)(N*rand_num[0]);
	    done=0;
	    while(done==0)
	    {
	      done=1;
	      for(int jtem=0; jtem<item; jtem++)
	      {

			if(Permut[k*N+jtem] == randunit)
			{
				randunit=(randunit +1)-(int)((randunit+1)/N)*N;
				jtem=item;
				done=0;
			}
	      }
	    }
	  Permut[k*N+item] = randunit;
	  item++;
	  }
	}
}
/*--------------------------------------------*/
void openfiles(int trial)
{
	char buffer[0x1000], buffer0[0x1000], buffer1[0x1000], buffer2[0x1000], buffer3[0x1000], buffer4[0x1000];
	
	snprintf(buffer, sizeof(buffer), "%s/Jmean_set%d", outdir, trial);
	Jm = fopen(buffer, "w");

	snprintf(buffer0, sizeof(buffer0), "%s/indicesL_set%d", outdir, trial);
	list_indices = fopen(buffer0, "w");

	snprintf(buffer1, sizeof(buffer1), "%s/UvsJ_set%d", outdir, trial);
	UvsJ = fopen(buffer1, "w");

	snprintf(buffer2, sizeof(buffer2), "%s/ld12_set%d", outdir, trial);
	ld12_file = fopen(buffer2, "w");

	snprintf(buffer3, sizeof(buffer3), "%s/recall_set%d", outdir, trial);
	recall = fopen(buffer3, "w");

	snprintf(buffer4, sizeof(buffer4), "%s/rehearse_set%d", outdir, trial);
	rehearse = fopen(buffer4, "w");
}

/*--------------------------------------------*/
void openfiles_detailed(int trial)
{
	char buffer5[0x1000], buffer6[0x1000], buffer7[0x1000], buffer8[0x1000], buffer9[0x1000], buffer10[0x1000];
	char buffer11[0x1000], buffer12[0x1000], buffer13[0x1000];

	snprintf(buffer5, sizeof(buffer5), "%s/mall_set%d", outdir, trial);
	mall = fopen(buffer5, "w");

	snprintf(buffer6, sizeof(buffer6), "%s/input_set%d", outdir, trial);
	input = fopen(buffer6, "w");

	snprintf(buffer7, sizeof(buffer7), "%s/sall_set%d", outdir, trial);
	sall = fopen(buffer7, "w");

	snprintf(buffer8, sizeof(buffer8), "%s/uall_set%d", outdir, trial);
	uall = fopen(buffer8, "w");

	snprintf(buffer9, sizeof(buffer9), "%s/xall_set%d", outdir, trial);
	xall = fopen(buffer9, "w");

	snprintf(buffer10, sizeof(buffer10), "%s/hall_set%d", outdir, trial);
	hall = fopen(buffer10, "w");

	snprintf(buffer11, sizeof(buffer11), "%s/Jall_set%d", outdir, trial);
	Jall = fopen(buffer11, "w");	

	snprintf(buffer12, sizeof(buffer12), "%s/Jij_set%d", outdir, trial);
	Jij = fopen(buffer12, "w");		

	snprintf(buffer13, sizeof(buffer13), "%s/test_set%d", outdir, trial);
	test = fopen(buffer13, "w");
}
/*--------------------------------------------*/

void closefiles()
{
	fclose(Jm);
	fclose(list_indices);
	fclose(UvsJ);
	fclose(ld12_file);
	fclose(recall);
	fclose(rehearse);
}

/*--------------------------------------------*/

void closefiles_detailed()
{
	fclose(mall);
	fclose(input);
	fclose(sall);
	fclose(hall);
	fclose(uall);
	fclose(xall);
	fclose(Jall);
	fclose(Jij);
	fclose(test);
}
/*--------------------------------------------*/

void getmemory()
{
	list_idx=malloc(2*L*sizeof(int));

	Permut= malloc(NumSet*N*sizeof(int*));

	xi = malloc(p*N*sizeof(int*));

	s = malloc(N*sizeof(double));
		
	h = malloc(N*sizeof(double));

	r = malloc(N*sizeof(double));

	u = malloc(N*sizeof(double));

	x = malloc(N*sizeof(double));

	g = malloc(N*sizeof(double));

	m = malloc(p*sizeof(double));

	smean = malloc(p*sizeof(double));

	hmean = malloc(p*sizeof(double));

	umean = malloc(p*sizeof(double));

	xmean = malloc(p*sizeof(double));

	Jmean = malloc(p*sizeof(double));

	J = malloc(N*N*sizeof(double*));

	JSTP = malloc(N*N*sizeof(double*));
}

/*--------------------------------------------*/

void deletememory()
{
	free(list_idx);

	free(Permut);

	free(xi);

	free(s);

	free(h);

	free(r);

	free(u);

	free(x);

	free(g);

	free(m);

	free(smean);

	free(hmean);

	free(umean);

	free(xmean);

	free(Jmean);

	free(J);

	free (JSTP);
}

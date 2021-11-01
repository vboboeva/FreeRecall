#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <string.h> 
#include <time.h>
/*#include <mpi.h>*/
#define _GNU_SOURCE  
#include <fenv.h> 
#include "random.h"
#include "global.h"
#include "params.h"
#include "functions.h"

void run_net(int set)
{		
	int k, n, n0, unit, Mumax, Mumaxold, Mumax2, dcounter, stoplearning;
	double t, rand_num[1], extfield, backfield, Mmax, Mmax2, d12, fL;

	time_t start_sim, end_sim;

	start_sim=time(0);

	/*to initialize the network at steady state*/
	initialize_net(); 

	printf("begin dynamics\n");

	t=0;
	n=0;

	int stimuli[2*L+2];
	int background[2*L+2];
	int Tupdates[2*L+2];
	int cues[2*L+2];

	stimuli[0]=1;
	background[0]=0;
	Tupdates[0]=T_start;
	cues[0]=list_idx[0];

	int c=0;
	for(int y=1;y<2*L;y++)
	{
		if (y%2==0) 
		{
			stimuli[y]=0;
			background[y]=0;
			Tupdates[y]=T_ISI;
			cues[y]=-1;
		}
		else
		{
			stimuli[y]=1;
			background[y]=0;
			Tupdates[y]=T_stim;
			cues[y]=list_idx[c];
			c=c+1;
		}			
	}
	stimuli[2*L]=0;
	background[2*L]=1; 
	Tupdates[2*L]=T_end; 
	cues[2*L]=-1;
	
	stimuli[2*L+1]=0;
	background[2*L+1]=0;
	Tupdates[2*L+1]=T_recall;
	cues[2*L+1]=-1;

	d12=0.;
	dcounter=0;
	fL=0.;
	Mumaxold=-1;

	/*printf("list presentation\n");*/

	/*------------------------START-----------------------------*/
	for(int i=0;i<2*L+2;i++)
	{
		printf("%d\n",stimuli[i]);
		n0=n;
		/*update whole net Runs times*/
		for(int ttt=0;ttt<Tupdates[i];ttt++) 
		{
			/*pick which random sequence to use for updating net*/
			ranlxd(rand_num, 1);
			k = (int)(NumSet*rand_num[0]);

			/*update whole net once*/
			for(int j=0;j<N;j++)
			{
				unit=Permut[k*N+j];
				
				/*external field used to cue the network with patterns*/
				extfield=0.;
				if (stimuli[i]==1) 
				{
					extfield=1.*(n>n0)*I_e*(xi[cues[i]*N+unit]==1);
				}	
				/*thermalize or else use to suppress activity through random input, either at pauses or at end*/
				ranlxd(rand_num, 1);
				backfield = (2.*rand_num[0]-1.)*(background[i]==1); 

				/*stop learning at the end of the list presentation*/
				stoplearning=0;
				if(Tupdates[i]==T_recall) stoplearning=1;

				if(i==0) backfield=I_b*rand_num[0]*(background[i]==1);
				
				/*update state of one unit for each state*/				
				update_state(n,unit,extfield,backfield,stoplearning); 

				/*Snapshot of the system at every dt step*/
				if((n%snap)==0)
				{
					t=dt*n/N;

					compute_m();

					if(set < NumSave)
					{		
						/*fprintf(input, "%f\t%d\t%f\t%f\n",t,cues[i],extfield,backfield);*/
						print_m(t);
						compute_s_h_u_x(t);
					}

					/*RECALL*/
					if(stimuli[i]==0 && background[i]==0 && Tupdates[i]==T_recall && cues[i]==-1)
					{
						/*printf("recall list\n");*/

						int temp[2];
						int* maxvals=findmax(temp);

						Mumax=maxvals[0];
						Mumax2=maxvals[1];
						if(Mumaxold!=Mumax && m[Mumax]>0.75) 
						{
							Mumaxold=Mumax;
							fprintf(recall, "%f\t%d\t%f\n", t, Mumax, m[Mumax]);
							
							printf("%f\t%d\t%f\n", t, Mumax, m[Mumax]);
							/*fflush(recall);*/ 
							d12 += m[Mumax]-m[Mumax2];
							dcounter++;
		   					for(int l = 0; l < L; l++)
						    {
						        if(list_idx[l] == Mumax) fL++;
						    }							
						}				

					}

				}
			n++;
			}
		}

	}	

	if (dcounter > 0)
	{
		d12 /= 1.*dcounter;
		fL /= 1.*dcounter;
	}

	fprintf(ld12_file, "%f\t%f\n", d12, fL);	

	end_sim=time(0);
	/*printf( "simulation end: %ld\n", end_sim);*/
	printf( "duration %ld seconds\n", end_sim-start_sim);	
}

/*---------------------------- START SIMULATION----------------------*/

int main(int argc, char **argv)
{

int seed = rlxd_seed();
/*rlxd_init(1,seed);*/
rlxd_init(1,1990);

int trial=atoi(argv[0]);

int pattset=atoi(argv[1]);

char command[0x1000];
snprintf(outdir, sizeof(outdir), "output_N%d_a%.2f_p%d_L%d_dJ%.4f_JEI%.1f", N, a, p, L, dJ, kappa);
snprintf(command, sizeof(command), "mkdir -p %s", outdir);
system(command);
printf("output directory --> %s\n", outdir);

getmemory();

/*Tables that contain lists of units to be updated randomly*/
SetUpTables();

/*Read p patterns from file*/
char filename[0x100];
snprintf(filename, sizeof(filename), "input/pattern_N%d_a%.2f", N, a);
pat=fopen(filename, "r");
if(pat==NULL) 
{
	printf("Error: could not open %S", filename);
}
else
{
	read_p_patterns(pattset*p);
}
fclose(pat);

/*Make p patterns*/
/*make_pattern();*/
construct_J();

/*if HPC=0 then the variable trial is ignored and the program enters a loop to collect stats*/
/*otherwise if HPC=1 then the variable trial is not ignored, but used to run many independent simulations*/

if (HPC == 0)
{
	for(int z=0;z<Runs;z++)
	{
		printf("list number %d\n",z);

		openfiles(z);

		if(z < NumSave)
		{
			openfiles_detailed(z);
		}			

		pick_Lrandom_patterns();

		construct_JSTP();

		run_net(z);

		closefiles();

		if(z < NumSave)
		{
			closefiles_detailed();
		}
	}
}

if (HPC == 1)
{
	openfiles(trial);

	if(trial < NumSave)
	{
		openfiles_detailed(trial);
	}			

	pick_Lrandom_patterns();

	construct_JSTP();

	run_net(trial);

	closefiles();
	if(trial < NumSave)
	{
		closefiles_detailed();
	}
	
}

deletememory();
}
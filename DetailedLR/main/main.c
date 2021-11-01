#include <stdio.h> /*input output*/
#include <stdlib.h> /*general purpose functions, including dynamic memory management, random number generation, communication with the environment, integer arithmetics, searching, sorting and converting.*/
#include <math.h> 
#include <string.h> /*non numeric characters*/
#include <time.h>
/*#include <mpi.h>*/
#define _GNU_SOURCE /* See feature_test_macros(7) */ 
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
	/*printf("simulation started: %ld \n", start_sim);*/
	
	/*print out weights across units*/
	if (set < NumSave) print_Jij();

	/*to initialize the network at steady state*/
	initialize_net(); 

	/*compute and print mean weights over all patterns before and after learning*/
	print_Jmean_before_after();

	printf("begin dynamics\n");

	t=0;
	n=0;

	int stimuli[2*L+2];
	int background[2*L+2];
	int Tupdates[2*L+2];
	int cues[2*L+2];

	stimuli[0]=0;
	background[0]=1;
	Tupdates[0]=T_start;
	cues[0]=-1;

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

			/*update auxilary variable g
			update_aux();*/

			/*update whole net once*/
			for(int j=0;j<N;j++)
			{
				unit=Permut[k*N+j];
				
				/*external field used to cue the network with patterns*/
				extfield=0.;
				if (stimuli[i]==1) 
				{
					extfield=1.*(n>n0)*I_e*(xi[cues[i]*N+unit]==1); /*exp(-(1.*(n-n0)/(5.*N)));*/
				}	
				/*thermalize or else use to suppress activity through random input, either at pauses or at end*/
				ranlxd(rand_num, 1);
				backfield = (2.*rand_num[0]-1.)*(background[i]==1); /*(2.*rand_num[0]-1.)*/

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
						fprintf(input, "%f\t%d\t%f\t%f\n",t,cues[i],extfield,backfield);
						print_m(t);
						compute_s_h_u_x(t);
					}

					if(i<2*L)
					{
						int temp[2];
						int* maxvals=findmax(temp);

						Mumax=maxvals[0];
						Mumax2=maxvals[1];

						if(Mumaxold!=Mumax && m[Mumax]>0.75) 
						{
							Mumaxold=Mumax;
							fprintf(rehearse, "%f\t%d\t%f\n", t, Mumax, m[Mumax]);
							/*fflush(rehearse);*/ 
						}										
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

					/*printf("%f\t%f\t%d\n", t, d12, dcounter);*/
				}
			n++;
			}
		}

		if (i == 2*L)
		{
			compute_u_J(t);
		}
	}	

	print_Jmean_before_after();


	if (set < NumSave) print_Jij();

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
rlxd_init(1,seed);
/*rlxd_init(1,1987);*/

int trial=atoi(argv[1]);

int pattset=atoi(argv[2]);

char command[0x1000];
snprintf(outdir, sizeof(outdir), "output_p%d_L%d_alpha%.4f_beta%.4f_gamma%.4f_taupre%.1f_kappa%.3f", p, L, alpha, beta, gamma, tau_pre, kappa);
snprintf(command, sizeof(command), "mkdir -p %s", outdir);
system(command);
printf("output directory --> %s\n", outdir);

getmemory();

/*Tables that contain lists of units to be updated randomly*/
SetUpTables();

/*Read connectivity from file
char buffer[0x100];
snprintf(buffer, sizeof(buffer), "input/Crand_N%d_C%d", N, Cm);
cij=fopen(buffer, "r");
read_C();
fclose(cij);*/

/*construct_C();*/

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

	run_net(trial);

	closefiles();
	if(trial < NumSave)
	{
		closefiles_detailed();
	}
	
}

deletememory();
}
#ifndef GLOBAL_H
#define GLOBAL_H

#include "params.h"

/*main and functions*/

char outdir[0x1000];

int *list_idx;
int *xi;	

double *JSTP;
double *J; 

double *s;	
double *m;
double *h;
double *r;
double *u;
double *x;
double *g;
double *umean;
double *xmean;
double *hmean;
double *Jmean;
int *Permut; 

FILE *pat;
FILE *mall;
FILE *sall;
FILE *xall;
FILE *uall;
FILE *hall;
FILE *Jall;
FILE *input;
FILE *Jij;
FILE *Jm;
FILE *list_indices;
FILE *recall;
FILE *rehearse;
FILE *ld12_file;
FILE *UvsJ;
FILE *test;

#endif

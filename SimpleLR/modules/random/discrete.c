#define DISCRETE_C

#include "random.h"

/*
 * Functions extracting random numbers with discrete values, either double
 * precision (discr_dble) or integers (discr_int). Arguments are:
 *  - r --> array to fill
 *  - x --> array containing the possible values
 *  - q --> array containing the corresponding probabilities
 *  - dim --> number of possible values ( = dimension of x and q)
 *  - n --> numbers to extract ( < dimension of r)
 */

void discr_dble(double r[], double x[], double q[], int dim, int n){
    int i,j;
    double qq;
    double u[n];

    /* normalize probability vector */
    qq = 0;
    for (i=0; i<dim; i++) qq += q[i];
    for (i=0; i<dim; i++) q[i] = q[i]/qq;

    ranlxd(u,n);   /* extract n random numbers in [0,1) */
    for (j=0; j<n; j++){
        /* "invert" discrete cumulative density function */
        qq = 0;
        i = 0;
        while (u[j]>=qq) qq += q[i++];
        r[j] = x[--i];
    
    }
}

/* Extracting random integers :
 *   - r --> vector to fill with random numbers ("output");
 *   - x --> possible values of the random numbers;
 *   - q --> probability vector (same dimension as x);
 *   - dim --> dimension of x and q (number of possible outcomes)
 *   - n --> number of random numbers (dimension of r)
 **/
void discr_int(int r[], int x[], double q[], int dim, int n){
    int i,j;
    double qq;
    double u[n];

    /* normalize probability vector */
    qq = 0;
    for (i=0; i<dim; i++) qq += q[i];
    for (i=0; i<dim; i++) q[i] = q[i]/qq;

    ranlxd(u,n);   /* extract n random numbers in [0,1) */
    for (j=0; j<n; j++){
        /* "invert" discrete cumulative density function */
        qq = 0;
        i = 0;
        while (qq<=u[j]){
             qq += q[i];
             i++;
        }
        r[j] = x[--i];
    
    }
}

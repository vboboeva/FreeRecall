#ifndef FUNCTIONS_H
#define FUNCTIONS_H

extern void print_m(double);
extern void print_Jmean_before_after();
extern void print_Jij();
extern int cmpf (const void *, const void *);
extern void make_pattern();
extern void read_p_patterns(int);
extern void read_pattern(int);
extern void read_all_patterns();
extern void pick_Lrandom_patterns();
extern void initialize_net();
extern void read_C();
extern void construct_C();
extern void construct_J();
extern void construct_JSTP();
extern void update_state(int, int, double, double, int);
extern void update_aux();
extern int* findmax(int*);
extern void compute_m();
extern void compute_s_h_u_x(double);
extern void compute_u_J();
extern void SetUpTables();
extern void getmemory();
extern void deletememory();
extern void openfiles(int);
extern void openfiles_detailed(int);
extern void closefiles();
extern void closefiles_detailed();
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sched.h>

#define EPS 1e-16

void print_matrix (double *a, int n, int m);
int  read_matrix  (const char *name, double *a, int n, int m);

void print_block_matrix (double *a, int n, int m);
int  block_read_matrix  (const char *name, double *a, int n, int m);

void init_matrix (double *a, int n);
void init_b      (double *a, double *b, int n);
int read_b (const char *name, double *b, int n);

void block_init_matrix (double *a, int n, int m);

double delta101  (double *x, int n);
double str_max   (double *a, int n);
double residual1 (double *a, int n, double *x, double *b);

void matrix_diff       (double *a, double *b, int n, int m);
void matrix_add_mult   (double *a, double *b, double *add, int n, int m, int k);
void minus_matrix_mult (double *a, double *b, double *res, int n, int m, int k);

int l_u_decomp  (double *a, int n, double *lu, double eps);
int l_u_solve   (double *lu, int n, double *b, double eps);
int l_u_inverse (double *lu, int n, double *inv, double eps);

double get_abs_time ();
double get_cpu_time ();

typedef struct args 
{ 
  double *matr;
  int n;
  int m;
  double *buf;
  int id;
  int p;
  pthread_t *thread_ids;
  pthread_barrier_t *bar;
  double *cpu_time;
} ARGS;

int block_l_u_decomp (ARGS *args);
int block_l_u_solve  (ARGS *args, double *b);
void *block_l_u_decomp_threaded (void *pa);

void *set_matrix_affinity (void *pa);

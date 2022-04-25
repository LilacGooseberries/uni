#include "matrix.h"

int main (int argc, char *argv[])
{
  int n, m, p, i;
  char *a_name = 0;
  char *b_name = 0;
  
  if ((argc != 4 && argc != 5 && argc != 6)  || 
      (n = atoi (argv[1])) == 0 || 
      (m = atoi (argv[2])) == 0 || 
      (p = atoi (argv[3])) == 0)
    {
      printf ("usage: %s <n> <m> <p> [matrix.txt] [rhs.txt]\n", argv[0]);
      return 1;
    }
  if (argc == 5)
    a_name = argv[4];
  else if (argc == 6)
    {
      a_name = argv[4];
      b_name = argv[5];
    }
  
  double *a = new double [n * n];
  double *b = new double [n];
  double *ans = new double [n];
  pthread_t *threads = new pthread_t [p];
  ARGS *args = new ARGS [p];
  double *cpu_time = new double [p + 1];
  if (!a || !b || !ans || !threads || !args || !cpu_time)
    {
      printf ("not enough mem!\n");
      delete [] a; delete [] b;  delete [] ans; delete [] threads; delete [] args; delete [] cpu_time;
      return 2;
    }

  pthread_barrier_t bar;
  pthread_barrier_init (&bar, 0, p);
  
  for (i = 0; i < p; i++)
    {
      args[i].matr = a;
      args[i].n = n;
      args[i].m = m;
      if (!(args[i].buf = new double [2 * m * m]))
        {
          printf ("not enough mem!\n");
          pthread_barrier_destroy (&bar);
          for (int j = 0; j < i; j++) delete [] args[j].buf;
          delete [] a; delete [] b;  delete [] ans; delete [] threads; delete [] args; delete [] cpu_time;
          return 2;
        }
      args[i].id = i;
      args[i].p = p;
      args[i].thread_ids = threads;
      args[i].bar = &bar;
      args[i].cpu_time = cpu_time;
    }
                                                              //распределение кусков матрицы по ядрам, занулением соотв. элементов
  for (i = 0; i < p; i++)
    {
      if (pthread_create (threads + i, 0, set_matrix_affinity, args + i))
        {
          printf ("cannot create thread %d\n", i);
          return -1;
        }
    }
  for (i = 0; i < p; i++)
    {
      if (pthread_join (threads[i], 0))
        {
          printf ("cannot wait thread %d\n", i);
          return -2;
        }
    }
                                                              
  if (a_name)
    {
      if (read_matrix (a_name, a, n, n) < 0)
        {
          printf ("read error\n");
          pthread_barrier_destroy (&bar);
          for (i = 0; i < p; i++) delete [] args[i].buf;
          delete [] a; delete [] b;  delete [] ans; delete [] threads; delete [] args; delete [] cpu_time;
          return 3;
        }
      
      if (b_name)
        {
          if (read_b (b_name, b, n) < 0)
            {
              printf ("read error\n");
              pthread_barrier_destroy (&bar);
              for (i = 0; i < p; i++) delete [] args[i].buf;
              delete [] a; delete [] b;  delete [] ans; delete [] threads; delete [] args; delete [] cpu_time;
              return 4;
            }
        }
      else
        init_b (a, b, n);
    }
  else
    {
      init_matrix (a, n);
      init_b (a, b, n);
    }
  
  printf ("A:");
  print_matrix (a, n, n);
  printf ("\nb:");
  print_matrix (b, n, 1);
  
  if (a_name) block_read_matrix (a_name, a, n, m);
  else        block_init_matrix (a, n, m);

  int ret = block_l_u_solve (args, b);
  
  if (!ret)
    {
      for (i = 0; i < n; i++) ans[i] = b[i];
      if (a_name)
        {
          read_matrix (a_name, a, n, n);
          if (b_name)
            {
              read_b (b_name, b, n);
            }
          else
            init_b (a, b, n);
        }
      else
        {
          init_matrix (a, n);
          init_b (a, b, n);
        }
      
      double res = residual1 (a, n, ans, b);
      double del = delta101 (ans, n);
      
      printf ("\nanswer:");
      print_matrix (ans, n, 1);
      printf ("\nresidual   = %e\n", res);
      if (!b_name)
        printf ("error norm = %e\n", del);
      printf ("\ncpu time per thread:\n");
      for (i = 0; i < p; i++) printf ("%.2f ", cpu_time[i]);
      printf ("\n\nastr time: %.2f\n", cpu_time[p]);
    }
  else
      printf ("\nalgorithm is not applicable\n");
    
  pthread_barrier_destroy (&bar);
  for (i = 0; i < p; i++) delete [] args[i].buf;
  delete [] a; delete [] b;  delete [] ans; delete [] threads; delete [] cpu_time; delete [] args;
  return 0;
}

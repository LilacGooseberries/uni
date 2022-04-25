#include "matrix.h"

#define ERR_OPEN -1
#define ERR_READ -2
#define ERR_DATA -3

double
get_cpu_time ()
{
  struct rusage u;
  getrusage (RUSAGE_THREAD, &u);
  return (double) u.ru_utime.tv_sec + (double) u.ru_utime.tv_usec / 1e6;
}

double
get_abs_time ()
{
  struct timeval tv;
  gettimeofday (&tv, 0);
  return (double) tv.tv_sec + (double) tv.tv_usec / 1e6;
}

int 
read_matrix (const char *name, double *a, int n, int m)
{
  FILE *fp; 
  int i, j;
  if (!(fp = fopen (name, "r")))
      return ERR_OPEN;
  
  for (i = 0; i < n; i++)
    {
      for (j = 0; j < m; j++)
        {
          if (fscanf (fp, "%lf", a+i*m+j) != 1)
            {
              if (feof (fp))
                {
                  fclose (fp);
                  return ERR_DATA;
                }
              
              fclose (fp);
              return ERR_READ;
            }
        }
    }
      
  fclose (fp);
  return 0;
}

int
read_b (const char *name, double *b, int n)
{
  FILE *fp; 
  int i;
  if (!(fp = fopen (name, "r")))
      return ERR_OPEN;
  
  for (i = 0; i < n; i++)
    {
      if (fscanf (fp, "%lf", b+i) != 1)
        {
          if (feof (fp))
            {
              fclose (fp);
              return ERR_DATA;
            }
          
          fclose (fp);
          return ERR_READ;
        }
    }
      
  fclose (fp);
  return 0;
}

void
init_b (double *a, double *b, int n)
{
  int i, j;
  double s;
  
  for (i = 0; i < n; i++)
    {
      s = 0;
      for (j = 0; j < n; j += 2)
          s += a[i*n+j];
      
      b[i] = s;
    }
}

int 
block_read_matrix (const char *name, double *a, int n, int m)
{
  FILE *fp;
  int k = n / m, l = n - k * m;
  if (!(fp = fopen (name, "r")))
      return ERR_OPEN;
  
  int i, j, block_i, block_j, block_pos, shift;
  for (int s = 0; s < n * n; s++)
    {
      i = s / n; j = s - i * n;
      block_i = i / m; block_j = j / m;
      block_pos = block_i * (k * m * m + m * l);
      if (block_i != k)
          block_pos += block_j * m * m;
      else
          block_pos += block_j * l * m;
      
      if (block_j != k)
          shift = (i % m) * m + j % m;
      else
          shift = (i % m) * l + j % m;
      //printf("\n(%d %d)    b(%d %d)    pos %d + %d\n", i, j, block_i, block_j, block_pos, shift);
      if (fscanf (fp, "%lf", a + block_pos + shift) != 1)
        { 
          if (feof (fp))
            {
              fclose (fp);
              return ERR_DATA;
            }
          
          fclose (fp);
          return ERR_READ;
        }
    }
      
  fclose (fp);
  return 0;
}

#define N_MAX 10
#define M_MAX 10

void
print_block_matrix (double *a, int n, int m)
{
  int i, j, s, t;
  int k = n / m, l = n - k * m;
  
  printf ("\n");
  
  for (t = 0; t < k; t++)
    {
      for (s = 0; s < m; s++)
        {
          for (i = 0; i < k; i++)
            {
              for (j = 0; j < m; j++)
                  printf ("%.3f ", a[i*m*m + j + s*m]);
            }
          for (j = 0; j < l; j++)
              printf ("%.3f ", a[k*m*m + j + s*l]);
          printf ("\n");
        }
      a += m*m*k + m*l;
    }
    for (s = 0; s < l; s++)
      {
        for (i = 0; i < k; i++)
          {
            for (j = 0; j < m; j++)
                printf ("%.3f ", a[i*l*m + j + s*m]);
          }
        for (j = 0; j < l; j++)
            printf ("%.3f ", a[k*l*m + j + s*l]);
        printf ("\n");
      }
}

void
print_matrix (double *a, int n, int m)
{
  int n_max = (n > N_MAX ? N_MAX : n), m_max = (m > M_MAX ? M_MAX : m), i, j;
  
  printf ("\n");
  for (i = 0; i < n_max; i++)
  {
      for (j = 0; j < m_max; j++)
          printf ("%.3f ", a[i*m+j]);
      
      printf ("\n");
  }
}

void 
init_matrix (double *a, int n)
{
  int i, j;
  
  for (i = 0; i < n; i++)
  {
      for (j = 0; j < n; j++)
          a[i * n + j] = n - ((i > j)? i : j)/*1. / (i + j + 1)*/;
  }
}

void 
block_init_matrix (double *a, int n, int m)
{
  int k = n / m, l = n - k * m;
  
  int i, j, block_i, block_j, block_pos, shift;
  for (int s = 0; s < n * n; s++)
    {
      i = s / n; j = s - i * n;
      block_i = i / m; block_j = j / m;
      block_pos = block_i * (k * m * m + m * l);
      if (block_i != k)
          block_pos += block_j * m * m;
      else
          block_pos += block_j * l * m;
      
      if (block_j != k)
          shift = (i % m) * m + j % m;
      else
          shift = (i % m) * l + j % m;

      a[block_pos + shift] = n - ((i > j)? i : j)/*1. / (i + j + 1)*/;
    }
}

double
delta101 (double *x, int n)
{
  int i;
  double max = -1, d;
  
  for (i = 0; i < n; i++)
    {
      d = fabs (x[i] - (i + 1) % 2);
      if (d > max) max = d;
    }
  
  return max;
}

double
str_max (double *a, int n)
{
  double max = -1, s;
  int i, j;
  
  for (i = 0; i < n; i++)
    {
      for (s = 0, j = 0; j < n; j++)
          s += fabs (a[n * i + j]);
      
      if (s > max)
          max = s;
    }
  
  return max;
}

double
residual1 (double *a, int n, double *x, double *b)
{
  int i, j; double max = -1, ax, d;
  
  for (i = 0; i < n; i++)
    {
      for (ax = 0, j = 0; j < n; j++)
          ax += a[n * i + j] * x[j];
      
      d = fabs (ax - b[i]);
      if (d > max) max = d;
    }
  
  return max;
}

void
matrix_diff (double *a, double *b, int n, int m)
{
  for (int i = 0; i < n * m; i++)
      a[i] -= b[i];
}

inline void
matrix_add_mult (double *a, double *b, double *add, int n, int m, int k)
{
  int i, j, t;
  
  double *p_a = a, *p_b = b, *p_c = add;
  double s00, s01, s02, s10, s11, s12, s20, s21, s22;
  for (i = 0; i + 3 < n; i += 3)
    { 
      for (j = 0; j + 3 < k; j += 3)
        {      
          p_b = b + j;
          s00 = s01 = s02 = s10 = s11 = s12 = s20 = s21 = s22 = 0;
          for (t = 0; t < m; t++)
            {  
              p_a = a + t;
              
              s00 += p_a[0] * p_b[0];
              s01 += p_a[0] * p_b[1];
              s02 += p_a[0] * p_b[2];
              p_a += m;
              
              s10 += p_a[0] * p_b[0];
              s11 += p_a[0] * p_b[1];
              s12 += p_a[0] * p_b[2];
              p_a += m;
            
              s20 += p_a[0] * p_b[0];
              s21 += p_a[0] * p_b[1];
              s22 += p_a[0] * p_b[2];
              
              p_b += k;
            }
          p_c = add + j;  
            
          p_c[0] += s00;
          p_c[1] += s01;
          p_c[2] += s02;
          
          p_c += k;
          p_c[0] += s10;
          p_c[1] += s11;
          p_c[2] += s12;
          
          p_c += k;
          p_c[0] += s20;
          p_c[1] += s21;
          p_c[2] += s22;
        }
      for ( ; j < k; j++)
        {
          p_b = b + j;
          
          s00 = s10 = s20 = 0;
          for (t = 0; t < m; t++)
            {
              p_a = a + t;
              
              s00 += p_a[0] * p_b[0];
              s10 += p_a[m] * p_b[0];
              s20 += p_a[2 * m] * p_b[0];
              
              p_b += k;
            }
          p_c = add + j; 
            
          p_c[0] += s00;
          p_c[k] += s10;
          p_c[2 * k] += s20;
        }
        
        a += 3 * m;
        add += 3 * k;
    }
    
  for ( ; i < n; i++)
    { 
      for (j = 0; j + 3 < k; j += 3)
        {
          p_b = b + j;
          
          s00 = s01 = s02 = 0;
          for (t = 0; t < m; t++)
            {
              p_a = a + t;
              
              s00 += p_a[0] * p_b[0];
              s01 += p_a[0] * p_b[1];
              s02 += p_a[0] * p_b[2];
              
              p_b += k;
            }
          p_c = add + j;
          
          p_c[0] += s00;
          p_c[1] += s01;
          p_c[2] += s02;
        }
      for ( ; j < k; j++)
        {
          p_b = b + j;
          
          s00 = 0;
          for (t = 0; t < m; t++)
            {
              s00 += a[t] * p_b[0];
              
              p_b += k;
            }
          
          p_c = add + j; 
          
          p_c[0] += s00;
        }
      
      a += m;
      add += k;
    }
}

inline void
minus_matrix_mult (double *a, double *b, double *res, int n, int m, int k)
{
  int i, j, t;
  
  double *p_a = a, *p_b = b, *p_c = res;
  double s00, s01, s02, s10, s11, s12, s20, s21, s22;
  for (i = 0; i + 3 < n; i += 3)
    { 
      for (j = 0; j + 3 < k; j += 3)
        {      
          p_b = b + j;
          s00 = s01 = s02 = s10 = s11 = s12 = s20 = s21 = s22 = 0;
          for (t = 0; t < m; t++)
            {  
              p_a = a + t;
              
              s00 += p_a[0] * p_b[0];
              s01 += p_a[0] * p_b[1];
              s02 += p_a[0] * p_b[2];
              p_a += m;
              
              s10 += p_a[0] * p_b[0];
              s11 += p_a[0] * p_b[1];
              s12 += p_a[0] * p_b[2];
              p_a += m;
            
              s20 += p_a[0] * p_b[0];
              s21 += p_a[0] * p_b[1];
              s22 += p_a[0] * p_b[2];
              
              p_b += k;
            }
          p_c = res + j;  
            
          p_c[0] = -s00;
          p_c[1] = -s01;
          p_c[2] = -s02;
          
          p_c += k;
          p_c[0] = -s10;
          p_c[1] = -s11;
          p_c[2] = -s12;
          
          p_c += k;
          p_c[0] = -s20;
          p_c[1] = -s21;
          p_c[2] = -s22;
        }
      for ( ; j < k; j++)
        {
          p_b = b + j;
          
          s00 = s10 = s20 = 0;
          for (t = 0; t < m; t++)
            {
              p_a = a + t;
              
              s00 += p_a[0] * p_b[0];
              s10 += p_a[m] * p_b[0];
              s20 += p_a[2 * m] * p_b[0];
              
              p_b += k;
            }
          p_c = res + j; 
            
          p_c[0] = -s00;
          p_c[k] = -s10;
          p_c[2 * k] = -s20;
        }
        
        a += 3 * m;
        res += 3 * k;
    }
    
  for ( ; i < n; i++)
    { 
      for (j = 0; j + 3 < k; j += 3)
        {
          p_b = b + j;
          
          s00 = s01 = s02 = 0;
          for (t = 0; t < m; t++)
            {
              p_a = a + t;
              
              s00 += p_a[0] * p_b[0];
              s01 += p_a[0] * p_b[1];
              s02 += p_a[0] * p_b[2];
              
              p_b += k;
            }
          p_c = res + j;
          
          p_c[0] = -s00;
          p_c[1] = -s01;
          p_c[2] = -s02;
        }
      for ( ; j < k; j++)
        {
          p_b = b + j;
          
          s00 = 0;
          for (t = 0; t < m; t++)
            {
              s00 += a[t] * p_b[0];
              
              p_b += k;
            }
          
          p_c = res + j; 
          
          p_c[0] = -s00;
        }
      
      a += m;
      res += k;
    }
}

int
block_l_u_decomp (ARGS *args)
{
  double *a = args->matr;
  int n = args->n;
  int m = args->m;
  double *buf = args->buf;
  int p = args->p;
  pthread_t *threads = args->thread_ids;
  int i;

  for (i = 0; i < p; i++)
    {
      if (pthread_create (threads + i, 0, block_l_u_decomp_threaded, args + i))
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
  
  int k = n / m, l = n - k * m;
  int x = k * m * m + m * l, y = m * m;
  if (!l)
    {
      l = m;
      k--;
    }

  if (args->id == -1)
      return -3;
                                                              //последний блок
  for (i = 0; i < l * l; i++) buf[i] = 0;
  for (i = 0; i < k; i++)
      matrix_add_mult (a + k * x + i * l * m,
                       a + i * x + k * y, buf, l, m, l);
    
  matrix_diff (a + k * x + k * l * m, buf, l, l);
  
  return 0;
}

inline int
l_u_decomp (double *a, int n, double *lu, double eps)
{
  int i, s, t;
  double sum;
  
  for (s = 0; s < n; s++)
    {
      for (i = s; i < n; i++)
        {
          for (t = 0, sum = 0; t < s; t++)
              sum += lu[i * n + t] * lu[t * n + s];
          
          lu[i * n + s] = a[i * n + s] - sum;
        }
        
      for (i = s + 1; i < n; i++)
        {
          for (t = 0, sum = 0; t < s; t++)
              sum += lu[s * n + t] * lu[t * n + i];
          
          if (lu[s * n + s] <= eps) return -1;
          lu[s * n + i] = (a[s * n + i] - sum) / lu[s * n + s];
        }
    }
    
  return 0;
}

int
block_l_u_solve (ARGS *args, double *b)
{
  double *a = args->matr;
  int n = args->n;
  int m = args->m;
  int p = args->p;
  double *buf = args->buf;
  
  int t, i, k = n / m, l = n - k * m;
  int x = k * m * m + m * l, y = m * m;
  double eps;
  
  if (!l)
    {
      l = m;
      k--;
    }
  
  if (block_l_u_decomp (args))
      return -1;
  
  double abs_time = get_abs_time ();
                                                     //LY = B, B <- Y
  for (t = 0; t < k; t++)
    {
      for (i = 0; i < m; i++) buf[i] = 0;
      for (i = 0; i < t; i++)
          matrix_add_mult (a + t * x + i * y,
                           b + i * m, buf, m, m, 1);
          
      matrix_diff (b + t * m, buf, m, 1);
      
      eps = str_max  (a + t * x + t * y, m) * EPS;
      if (l_u_decomp (a + t * x + t * y, m, buf, eps) ||
          l_u_solve (buf, m, b + t * m, eps)            )
          return -1;
    }
  for (i = 0; i < l; i++) buf[i] = 0;
  for (i = 0; i < k; i++)
      matrix_add_mult (a + k * x + i * l * m,
                       b + i * m, buf, l, m, 1);
      
  matrix_diff (b + k * m, buf, l, 1);
  
  eps = str_max  (a + k * x + k * l * m, l) * EPS;
  if (l_u_decomp (a + k * x + k * l * m, l, buf, eps) ||
      l_u_solve (buf, l, b + k * m, eps)                )
      return -1;
                                                     //UX = B, B <- X
  for (t = k - 1; t >= 0; t--)
    {
      for (i = 0; i < m; i++) buf[i] = 0;
      matrix_add_mult (a + t * x + k * y,
                       b + k * m, buf, m, l, 1);
      for (i = k - 1; i > t; i--)
          matrix_add_mult (a + t * x + i * y,
                           b + i * m, buf, m, m, 1);
          
      matrix_diff (b + t * m, buf, m, 1);
    }
    
  abs_time = get_abs_time () - abs_time;
  args->cpu_time[p] += abs_time;
    
  return 0;
}

inline int
l_u_solve (double *lu, int n, double *b, double eps)
{
  int i, t;
  double sum;
  
  for (i = 0; i < n; i++)
    {
      for (t = 0, sum = 0; t < i; t++)
          sum += b[t] * lu[i * n + t];
      
      if (lu[i * n + i] <= eps) return -1;
      b[i] = (b[i] - sum) / lu[i * n + i];
    }
    
  for (i = n - 1; i >= 0; i--)
    {
      for (t = n - 1, sum = 0; t > i; t--)
          sum += b[t] * lu[i * n + t];
      
      b[i] -= sum;
    }
    
    return 0;
}

inline int
l_u_inverse (double *lu, int n, double *inv, double eps)
{
  int i, t, col;
  double sum;
  
  for (col = 0; col < n; col++)
  {
    for (i = 0; i < col; i++) inv[i * n + col] = 0;
    if (lu[col * n + col] <= eps) return -1;
    inv[col * n + col] = 1 / lu[col * n + col];
    for (i = col + 1; i < n; i++)
      {
        for (t = col, sum = 0; t < i; t++)
            sum += inv[t * n + col] * lu[i * n + t];
        
        if (lu[i * n + i] <= eps) return -1;
        inv[i * n + col] = -sum / lu[i * n + i];
      }
      
    for (i = n - 1; i >= 0; i--)
      {
        for (t = n - 1, sum = 0; t > i; t--)
            sum += inv[t * n + col] * lu[i * n + t];
        
        inv[i * n + col] -= sum;
      }
  }
  
  return 0;
}

void *
block_l_u_decomp_threaded (void *pa)
{
  ARGS *pargs = (ARGS *) pa;
  double *a = pargs->matr;
  int n = pargs->n;
  int m = pargs->m;
  double *buf = pargs->buf;
  double *inv_buf = buf + m * m;
  int id = pargs->id;
  int p = pargs->p;
  pthread_barrier_t *bar = pargs->bar;
  int shift;
  
  double time = get_cpu_time ();
  double abs_time;
  if (id == 0) abs_time = get_abs_time ();
  int i, s, t;
  int k = n / m, l = n - k * m;
  int x = k * m * m + m * l, y = m * m;     // "координаты" для начала блоков (для первых k-1 строчек)
  if (!l)
    {
      l = m;
      k--;
    }
  double eps;

  cpu_set_t cpu;
  CPU_ZERO (&cpu);
  CPU_SET (id, &cpu);
  pthread_setaffinity_np (pthread_self (), sizeof (cpu), &cpu);

  for (s = 0; s < k; s++)
    {                                                         //L
      shift = (s + p - 1 - id) / p;
      for (i = id + p * shift; i < k; i += p)
        {
          for (t = 0; t < m * m; t++) buf[t] = 0;
          for (t = 0; t < s; t++)
              matrix_add_mult (a + i * x + t * y,
                               a + t * x + s * y, buf, m, m, m);
          
          matrix_diff (a + i * x + s * y, buf, m, m);
        }
                                                                        //L: i = k
      if (i == k)
      {
        for (t = 0; t < l * m; t++) buf[t] = 0;
        for (t = 0; t < s; t++)
            matrix_add_mult (a + k * x + t * l * m,
                             a + t * x + s * y, buf, l, m, m);
        
        matrix_diff (a + k * x + s * l * m, buf, l, m);
      }
                                                          //______________________
      pthread_barrier_wait (bar);
                                                              //U
      eps = str_max (a + s * x + s * y, m) * EPS;
      if (l_u_decomp (a + s * x + s * y, m, buf, eps) ||
          l_u_inverse (buf, m, inv_buf, eps)            )
        {
          pargs->id = -1;
          break;
        }
      
      shift = (s + p - id) / p;
      for (i = id + p * shift; i < k; i += p)
        {
          for (t = 0; t < m * m; t++) buf[t] = 0;
          for (t = 0; t < s; t++)
              matrix_add_mult (a + s * x + t * y,
                               a + t * x + i * y, buf, m, m, m);
          
          matrix_diff (buf, a + s * x + i * y, m, m);
          minus_matrix_mult (inv_buf, buf, a + s * x + i * y, m, m, m);
        }
                                                                        //U: i = k
        if (i == k)
          {
            for (t = 0; t < m * l; t++) buf[t] = 0;
            for (t = 0; t < s; t++)
                matrix_add_mult (a + s * x + t * y,
                                a + t * x + k * y, buf, m, m, l);
            
            matrix_diff (buf, a + s * x + k * y, m, l);
            minus_matrix_mult (inv_buf, buf, a + s * x + k * y, m, m, l);
          }
                                                          //______________________
      pthread_barrier_wait (bar);
    }
    
  if (id == 0)
    {
      abs_time = get_abs_time () - abs_time;
      pargs->cpu_time[p] = abs_time;
    }
  time = get_cpu_time () - time;
  pargs->cpu_time[id] = time;
  
  return 0;
}

void *
set_matrix_affinity (void *pa)
{
  ARGS *pargs = (ARGS *) pa;
  double *a = pargs->matr;
  int n = pargs->n;
  int m = pargs->m;
  int id = pargs->id;
  int p = pargs->p;
  
  int i, s, j, shift;
  int k = n / m, l = n - k * m;
  int x = k * m * m + m * l, y = m * m;    // "координаты" для начала блоков (для первых k-1 строчек)
  double *block_pos;
  
  cpu_set_t cpu;
  CPU_ZERO (&cpu);
  CPU_SET (id, &cpu);
  pthread_setaffinity_np (pthread_self (), sizeof (cpu), &cpu);
  
  for (s = 0; s < k; s++)
    {                                                         //L
      shift = (s + p - 1 - id) / p;
      for (i = id + p * shift; i < k; i += p)
        {     
          block_pos = a + i * x + s * y;
          for (j = 0; j < m * m; j++) block_pos[j] = 0;
        }
                                                                        //L: i = k
      if (i == k)
        {     
          block_pos = a + k * x + s * l * m;
          for (j = 0; j < l * m; j++) block_pos[j] = 0;
        }
                                                              //U
      shift = (s + p - id) / p;
      for (i = id + p * shift; i < k; i += p)
        {
          block_pos = a + s * x + i * y;
          for (j = 0; j < m * m; j++) block_pos[j] = 0;
        }
                                                                        //U: i = k
      if (i == k)
        {     
          block_pos = a + s * x + k * y;
          for (j = 0; j < m * l; j++) block_pos[j] = 0;
        }
    }
  
  return 0;
}



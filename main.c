#include <assert.h>
#include <getopt.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include "pthread_barrier.h"

#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)

#define DUMMY_5x5 {1, 1, 1, 1, 1, 1, 3, 7, 2, 1, 1, 8, 6, 5, 1, 1, 9, 0, 4, 1, 1, 1, 1, 1, 1}

struct args_t 
{
    int verbosity;    // -v
    int dimension;    // -d
    double precision; // -p
    int threads;      // -t
} args;

typedef struct thread_work
{
    int t_id;
    double** read_from;
    double** write_to;
    int start_row;
    int start_col;
    int finish_row;
    int ncells;
    int nrows;
    int ncols;
    double precision;
} thread_work_t;

pthread_barrier_t worker_threads_done;
pthread_barrier_t main_thread_done;

int global_continue = 0;
int global_done = 0;

int minimum(int a, int b) 
{ 
    return a < b ? a : b; 
}

void populate_random(double **one, double **two, int nrows, int ncols) 
{
    srand((unsigned)time(NULL));
    for (int i=0; i<nrows; ++i)
    {
        for (int j=0; j<ncols; ++j)
        {
            double r = 10*(double)rand()/(double)RAND_MAX;
            one[i][j] = r;
            two[i][j] = r;
        }
    }
}

void populate_values(double vals[],
                    int nvals,
                    double **one, 
                    double **two, 
                    int nrows, 
                    int ncols)
{
    // Make sure we can populate the whole array from vals[]
    assert(nvals == (nrows*ncols));
    int v = 0;
    for (int i=0; i<nrows; ++i)
    {
        for (int j=0; j<ncols; ++j)
        {
            one[i][j] = (double)vals[v];
            two[i][j] = (double)vals[v];
            ++v;
        }
    }
}

void print_matrix(double **arr, int nrows, int ncols) 
{
    for (int i=0; i<nrows; ++i)
    {
        for (int j=0; j<ncols; ++j)
        {
            printf("%f ", arr[i][j]);        
        }
        printf("\n");
    }
    printf("\n");
}

int relax_cell(double **in, double **out, int row, int col, double precision)
{
    double adjac = in[row-1][col]   //N
        + in[row+1][col]            //S
        + in[row][col-1]            //W
        + in[row][col+1];           //E
    double old = in[row][col];
    double new = adjac/(double)4;
    out[row][col] = new;
    /* If the the value differences are precise enough, return 0
     * Otherwise return 1
     */
    // if (args.verbosity)
    // {
    //     printf("Delta: %f \n", fabs(old-new));
    // }
    return (fabs(old-new) >= precision);
}

void swap(double ***one, double ***two)
{
    double **tmp = *one;
    *one = *two;
    *two = tmp;
}

void *relax_section(void *work_ptr)
{
    thread_work_t *work = (thread_work_t*)work_ptr;
    for (int iter=0;;++iter)
    {
        int local_continue = 0;
        int initial_jump = 1;
        int cells_done = 0;
        for (int i=1; i<work->nrows-1 && cells_done!=work->ncells; ++i)
        {
            for (int j=1; j<work->ncols-1 && cells_done!=work->ncells; ++j)
            {
                if (unlikely(initial_jump))
                {
                    i = work->start_row;
                    j = work->start_col;
                    initial_jump = 0;
                }
                local_continue += relax_cell(work->read_from, work->write_to, i, j, work->precision);
                ++cells_done;
            }
        }
        if (likely(local_continue))
        {
            ++global_continue;
        }
        int wait = pthread_barrier_wait(&worker_threads_done);
        if (unlikely(wait != 0 && wait != PTHREAD_BARRIER_SERIAL_THREAD))
        {
            exit(-1);
        }
        int main_wait = pthread_barrier_wait(&main_thread_done);
        if (unlikely(main_wait != 0 && main_wait != PTHREAD_BARRIER_SERIAL_THREAD))
        {
            exit(-1);
        }
        if (global_done)
        {
            return NULL;
        }
        if (work->t_id==0)print_matrix(work->write_to, work->ncols, work->nrows);
        swap(&work->read_from, &work->write_to);
    }
}

void assign_thread_work(thread_work_t *all_work,
                        double **from,
                        double **to,
                        int cell_split, 
                        int cell_remainder,
                        int nrows,
                        int ncols,
                        int nthreads, 
                        double precision)
{
    int next_col=0;
    int next_row=0;
    for (int t=0; t<nthreads; ++t)
    {
        thread_work_t work = 
        {
            .t_id = t,
            .read_from = from,
            .write_to = to,
            .start_row = next_row+1,
            .start_col = next_col+1,
            .ncells = cell_split + (t<cell_remainder ? 1:0),
            .nrows = nrows,
            .ncols = ncols,
            .precision = precision,
        };
        if (args.verbosity)
        {
            printf("Thread %d starting at (%d,%d) doing %d cells\n", 
            t, 
            work.start_row, 
            work.start_col, 
            work.ncells
        );
        int next = next_col + (work.ncells);
        next_row += ((next) / (ncols-2));
        next_col = ((next) % (ncols-2));
        all_work[t] = work;
        }
    }
}

void process_options(int argc, char **argv)
{
    const char *opt_string = "vd:p:t:";
    int opt = getopt(argc, argv, opt_string);
    while( opt != -1 ) 
    {
        switch( opt ) 
        {
            case 'v':
                args.verbosity = 1;
                break;
            case 'd':
                args.dimension = atoi(optarg);
                break;
            case 'p':
                args.precision = atof(optarg);
                break;
            case 't':
                args.threads = atoi(optarg);
                break;
            default:
                break;
        }
        opt = getopt(argc, argv, opt_string);
    }
}

int relax_matrix(thread_work_t *work, 
                 int nthreads)
{
    int iter_count = 0;
    pthread_t threads[nthreads];
    /* Initialise the barrier to nthreads+1 since the main thread
     * will hit both but does not count as a relaxation 'worker'
     */
    pthread_barrier_init(&worker_threads_done, NULL, nthreads+1);
    pthread_barrier_init(&main_thread_done, NULL, nthreads+1);
    
    for (int t=0; t<nthreads; ++t)
    {
        if (pthread_create(&threads[t], NULL, relax_section, &work[t]))
        {
            fprintf(stderr, "Failed creating thread %d\n", t);
            return 1;
        }
    }
    
    while (!global_done)
    {
        /* Wait for all the matrix workers to relax, then increment the
         * iteration count and calculate the precision.
         */
        pthread_barrier_wait(&worker_threads_done);
        ++iter_count;
        if (global_continue)
        {
            /* If any thread needs to continue it increments this, which needs
             * to be reset before they iterate again.*/
            global_continue = 0;
        }
        else
        {
            global_done = 1;
        }
        pthread_barrier_wait(&main_thread_done);
    }

    for (int t=0; t<nthreads; ++t)
    {
        if (pthread_join(threads[t], NULL))
        {
            fprintf(stderr, "Failed joining thread %d\n", t);
            return 1;
        }
    }
    
    if (args.verbosity)
    {
        printf("Reached in %d iterations.\n", iter_count);
    }
    
    //print_matrix(work[0].read_from, args.dimension, args.dimension);
    print_matrix(work[0].write_to, args.dimension, args.dimension);
    
    pthread_barrier_destroy(&worker_threads_done);
    pthread_barrier_destroy(&main_thread_done);
    return 0;
}

int main(int argc, char **argv)
{
    /* If any of the args are invalid, fall back to these defaults */
    args.verbosity = 0;
    args.dimension = 5;
    args.precision = 0.75;
    args.threads = 1;

    process_options(argc, argv);

    int nrows = args.dimension, ncols = args.dimension;
    double precision = args.precision;
    int nthreads = args.threads;

    thread_work_t* all_work = malloc(nthreads*sizeof(thread_work_t));
    if (unlikely(all_work==NULL))
    {
        return 1;
    }
    
    /* Allocate the read-from and write-to arrays */
    double **from = malloc(nrows*sizeof(double*));
    double **to   = malloc(nrows*sizeof(double*));
    double *from_buf = malloc(nrows*ncols*sizeof(double));
    double *to_buf   = malloc(nrows*ncols*sizeof(double));
    if (unlikely(from==NULL || from_buf== NULL || to==NULL || to_buf== NULL))
    {
        return 1;
    }
    for (int i=0; i<nrows; ++i)
    {
        from[i] = from_buf + ncols*i;
        to[i] = to_buf + ncols*i;
    }

    double dummy[25] = DUMMY_5x5;
    populate_values(dummy, 25, from, to, nrows, ncols);
    //populate_random(from, to, nrows, ncols);


    /* The total number of rows to be worked on is nrows minus the first 
     * and last. Each thread does a minimum of (nrows-2)/nthreads and a
     * maximum of ((nrows-2)/nthreads)+1 rows, depending on the remainder.
     */
    int cell_split = ((nrows-2)*(ncols-2))/nthreads;
    int cell_remainder = ((nrows-2)*(ncols-2))%nthreads;
    assign_thread_work(
        all_work, 
        from, to, 
        cell_split, cell_remainder, 
        nrows, ncols, 
        nthreads, precision
    );

    int err = relax_matrix(all_work, nthreads);
    
    free(from);
    free(to);
    free(from_buf);
    free(to_buf);
    free(all_work);
    return err;
}
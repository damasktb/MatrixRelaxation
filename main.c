#include <assert.h>
#include <getopt.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

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

pthread_barrier_t thread_work_barrier;

/* Set by threads which need to continue */
volatile int global_continue = 0;
/* Set by the main thread only, if all threads are done */
volatile int global_done = 0;

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

/* Populate initial arrays **one and **two identically with random floats 0-10.
 * To generate the same numbers each time, comment out the srand() call
 */
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

/* When passed references to two arrays, this will swap them */
void swap(double ***one, double ***two)
{
    double **tmp = *one;
    *one = *two;
    *two = tmp;
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
    /* If the the value delta is < precision, return 0, otherwise return 1.
     * Return values from this function can then be summed; any positive
     * value then indicates that the relaxation is incomplete.
     */
    return (precision < fabs(old-new));
}

/* This is the thread work function, which will iterate until the main thread
 * signals termination by setting global_done to 1
 */
void *relax_section(void *work_ptr)
{
    thread_work_t *work = (thread_work_t*)work_ptr;
    while (1)
    {
        /* 0 if all cells are relaxed 'enough', 1 otherwise */
        int local_continue = 0;
        int initial_jump = 1;
        /* We exit the loops if cells_done==work.ncells */
        int cells_done = 0;

        for (int i=1; i<work->nrows-1 && cells_done!=work->ncells; ++i)
        {
            for (int j=1; j<work->ncols-1 && cells_done!=work->ncells; ++j)
            {
                /* Threads don't all start at the beginning of a row, so if it's
                 * the first iteration, jump to the actual starting indices.
                 * The branch predictor will catch on to this very quickly.
                 */
                if (initial_jump)
                {
                    i = work->start_row;
                    j = work->start_col;
                    initial_jump = 0;
                }
                local_continue += relax_cell(work->read_from, work->write_to, 
                                            i, j, work->precision);
                ++cells_done;
            }
        }
        if (local_continue)
        {
            /* The data race here is deliberate; see report.pdf pages 4-5 */
            global_continue = 1;
        }
        int wait = pthread_barrier_wait(&thread_work_barrier);
        if (wait != 0 && wait != PTHREAD_BARRIER_SERIAL_THREAD)
        {
            exit(-1);
        }
        wait = pthread_barrier_wait(&thread_work_barrier);
        if (wait != 0 && wait != PTHREAD_BARRIER_SERIAL_THREAD)
        {
            exit(-1);
        }
        /* If the main thread signals termination, we can return */
        if (global_done)
        {
            return NULL;
        }
        /* Otherwise, another iteration is required */
        if (args.verbosity && work->t_id==0)
        {
            /* Print the current matrix state (from one thread only) */
            print_matrix(work->write_to, work->ncols, work->nrows);
        }
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
    /* Set the start location and number of cells to work on for each thread */
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
        }
        int next = next_col + (work.ncells);
        next_row += ((next) / (ncols-2));
        next_col = ((next) % (ncols-2));
        all_work[t] = work;
    }
}

int relax_matrix(thread_work_t *work, int nthreads)
{
    int iter_count = 0;
    pthread_t threads[nthreads];
    /* Initialise the barrier to nthreads+1 since the main thread will hit
     * both but does not do any of the relaxation work
     */
    pthread_barrier_init(&thread_work_barrier, NULL, nthreads+1);
    
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
        /* Wait for all the workers to relax once, then check if we can stop. */
        pthread_barrier_wait(&thread_work_barrier);
        ++iter_count;
        if (global_continue != 0)
        {
            /* If any thread needs to continue it sets this to 1. This needs
             * to be reset before they iterate again.*/
            global_continue = 0;
        }
        else
        {
            global_done = 1;
        }
        pthread_barrier_wait(&thread_work_barrier);
    }

    for (int t=0; t<nthreads; ++t)
    {
        if (pthread_join(threads[t], NULL))
        {
            fprintf(stderr, "Failed joining thread %d\n", t);
            return 1;
        }
    }
    
    printf("Reached in %d iterations.\n", iter_count);
    if (args.verbosity)
    {
        /* Print the finished matrix */
        print_matrix(work[0].write_to, args.dimension, args.dimension);
    }
    
    pthread_barrier_destroy(&thread_work_barrier);
    return 0;
}

void process_options(int argc, char **argv)
{
    const char *opt_string = "vd:p:t:";
    int opt = getopt(argc, argv, opt_string);
    while(opt != -1) 
    {
        switch(opt) 
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

int main(int argc, char **argv)
{
    /* If any of the args are invalid, fall back to these defaults */
    args.verbosity = 0;
    args.dimension = 100;
    args.precision = 0.5;
    args.threads = 1;
    process_options(argc, argv);

    /* This implmentation supports rectangular matrices */
    int nrows = args.dimension;
    int ncols = args.dimension;
    double precision = args.precision;
    int nthreads = args.threads;

    /* Allocate an empty thread_work_t struct for each thread  */
    thread_work_t* all_work = malloc(nthreads*sizeof(thread_work_t));
    if (all_work==NULL)
    {
        return 1;
    }
    
    /* Allocate the read-from write-to solution arrays */
    double **from = malloc(nrows*sizeof(double*));
    double **to   = malloc(nrows*sizeof(double*));
    double *from_buf = malloc(nrows*ncols*sizeof(double));
    double *to_buf   = malloc(nrows*ncols*sizeof(double));
    if (from==NULL || from_buf== NULL || to==NULL || to_buf== NULL)
    {
        return 1;
    }
    for (int i=0; i<nrows; ++i)
    {
        from[i] = from_buf + ncols*i;
        to[i] = to_buf + ncols*i;
    }

    populate_random(from, to, nrows, ncols);
    /* The total number of cells to be worked on is the full matrix minus
     * its outer values. Each thread does a minimum of cell_split and a
     * maximum of cell_split+1 cells, depending on the remainder.
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
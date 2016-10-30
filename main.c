#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)

#define DUMMY_3x3 {1, 1, 1, 1, 3, 1, 1, 1, 1}
#define DUMMY_4x4 {1, 1, 1, 1, 1, 4, 4, 1, 1, 4, 4, 1, 1, 1, 1, 1}

typedef struct thread_work
{
    double** read_from;
    double** write_to;
    int start_row;
    int finish_row;
    int nrows;
    int ncols;
} thread_work_t;

int within_precision(double **one, 
                    double **two, 
                    int rows, 
                    int cols, 
                    double prec)
{
    for (int i=1; i<rows-1; ++i)
    {
        for (int j=1; j<cols-1; ++j)
        {
            if (fabs(one[i][j]-two[i][j]) > prec)
            {
                return 0;
            }
        }
    }
    return 1;
}

void relax_cell(double **in, double **out, int row, int col)
{
    double adjac = in[row-1][col]   //N
        + in[row+1][col]            //S
        + in[row][col-1]            //W
        + in[row][col+1];           //E
    out[row][col] = adjac/(double)4;
}

// void swap(double ***one, double ***two)
// {
//     double **tmp = *one;
//     *one = *two;
//     *two = tmp;
// }

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

void populate_values(double vals[]
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

void *relax_section(void *work_ptr)
{
    thread_work_t *work = (thread_work_t*)work_ptr;
    for (int i=work->start_row; i<=work->finish_row; ++i)
    {
        for (int j=1; j<work->ncols-1; ++j)
        {
            relax_cell(work->read_from, work->write_to, i, j);
        }
    }
    return NULL;
}

int main(void)
{
    int nrows = 7;
    double precision = (double)0.05;
    int nthreads = 4;

    printf("Enter array dimension: ");
    scanf("%d", &nrows);

    printf("Enter precision: ");
    scanf("%lf", &precision);

    printf("Enter number of threads: ");
    scanf("%d", &nthreads);

    int ncols = nrows;
    int work_rows = (nrows-2)/nthreads;
    int work_remainder = (nrows-2)%nthreads;
    thread_work_t* all_work = malloc(nthreads*sizeof(thread_work_t));

    if (unlikely(all_work==NULL))
    {
        return 1;
    }
    
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

    populate_random(from, to, nrows, ncols);

    int next_start=1;
    for (int t=0; t<nthreads; ++t)
    {
        thread_work_t work = 
        {
            .read_from = from,
            .write_to = to,
            .start_row = next_start,
            .finish_row = next_start+work_rows-1 + (t<work_remainder ? 1:0),
            .nrows = nrows,
            .ncols = ncols
        };
        next_start = work.finish_row+1;
        printf("Thread %d doing rows %d-%d\n", t, work.start_row, work.finish_row);
        all_work[t] = work;
    }


    pthread_t threads[nthreads];
    int relaxed = 0, count = 0;
    while (!relaxed)
    {
        for (int t=0; t<nthreads; ++t)
        {   
            if (count%2)
            {
                all_work[t].read_from = to;
                all_work[t].write_to = from;
            }
            else
            {
                all_work[t].read_from = from;
                all_work[t].write_to = to;
            }
            if (pthread_create(&threads[t], NULL, relax_section, &all_work[t]))
            {
                fprintf(stderr, "Failed creating thread %d\n", t);
                return 1;
            }
        }
        for (int t=0; t<nthreads; ++t)
        {
            if (pthread_join(threads[t], NULL))
            {
                fprintf(stderr, "Failed joining thread %d\n", t);
                return 1;
            }
        }
        relaxed = within_precision(from, to, nrows, ncols, precision);
        ++count;
    }
    printf("Reached in %d iterations.\n", count);

    free(from);
    free(to);
    free(from_buf);
    free(to_buf);
    free(all_work);
    return 0;
}
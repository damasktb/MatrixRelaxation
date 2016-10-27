#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>


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

void swap(double ***one, double ***two)
{
    double **tmp = *one;
    *one = *two;
    *two = tmp;
}

void populate_random(double **one, double **two, int nrows, int ncols) 
{
    srand((unsigned)time(NULL));
    for (int i=0; i<nrows; ++i)
    {
        for (int j=0; j<ncols; ++j)
        {
            one[i][j] = 10*(double)rand()/(double)RAND_MAX;
            two[i][j] = one[i][j];
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

int main(void)
{
    int nrows = 7;
    int ncols = nrows;
    double precision = (double)0.05;
    
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

    int relaxed = 0, count = 0;
    while (!relaxed)
    {
        for (int i=1; i<nrows-1; ++i)
        {
            for (int j=1; j<ncols-1; ++j)
            {
                relax_cell(from, to, i, j);
            }
        }
        relaxed = within_precision(from, to, nrows, ncols, precision);
        swap(&from, &to);
        ++count;
    }

    printf("Reached in: %d iterations.\n", count);
    
    free(from);
    free(from_buf);
    return 0;
}
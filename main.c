#include <getopt.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define ROOT_PROC 0 // The root process

struct args_t 
{
    int verbosity;      // -v
    int dimension;      // -d
    double precision;   // -p
} args;

// Determine if index i lies on an edge of a matrix with dimension nrows*ncols
int is_edge(int nrows, int ncols, int i)
{
    return (i < ncols                // Top row
        ||  i+ncols >= (nrows*ncols) // Bottom row
        ||  i%ncols == 0             // Left edge
        ||  (i+1)%ncols == 0 );      // Right edge
}

// Print a 2d nrows*ncols matrix stored in a contiguous 1d array
void print_matrix(double *arr, int nrows, int ncols) 
{   
    for (int cell=0; cell<nrows*ncols; ++cell)
    {
        printf("%f ", arr[cell]);
        if ((cell+1) % ncols==0) // end of a row
        {
            printf("\n"); 
        }
    }
    printf("\n");
}

// Populate initial array with ones on the edges and zeroes in the centre
void populate_square(double *m, int nrows, int ncols) 
{
    for (int cell=0; cell<nrows*ncols; ++cell)
    {
        m[cell] = is_edge(nrows, ncols, cell) ? 1 : 0;
    }
}

int relax_cell(double *in, double *out, int cell, int ncols, double precision)
{
    double adjac = in[cell-ncols]   //N
        + in[cell+ncols]            //S
        + in[cell-1]                //W
        + in[cell+1];               //E
    double old = in[cell];
    double new = adjac * 0.25;
    out[cell-ncols] = new;
    /* If the the value delta is < precision, return 0, otherwise return 1.
     * Return values from this function can then be summed; any positive
     * value indicates that the relaxation is incomplete. */
    return (precision < fabs(old-new));
}

void process_options(int argc, char **argv)
{
    const char *opt_string = "vd:p:t:n:";
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
            default:
                break;
        }
        opt = getopt(argc, argv, opt_string);
    }
}

int main(int argc, char** argv)
{
    // If any of the args are invalid/unspecified, fall back to these defaults
    args.verbosity = 0;
    args.dimension = 10000;
    args.precision = 0.01;
    process_options(argc, argv);

    // This implmentation supports rectangular matrices if you override optarg
    const int nrows = args.dimension;
    const int ncols = args.dimension;
    const double precision = args.precision;

    MPI_Init(NULL, NULL);
    int world_rank; // The rank (id) of this process
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size; // The total number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // The full matrix will only be allocated on the root process
    double *global_matrix = NULL;
    if (world_rank == ROOT_PROC) {
        global_matrix = malloc(nrows*ncols*sizeof(double));
        if (global_matrix == NULL)
        {
            return 1;
        }
        populate_square(global_matrix, nrows, ncols);
    }

    // Each process will do at least this many cells
    const int cell_split = ncols * ((nrows-2)/world_size);  
    // These processes will do an extra row
    const int extra_row = (nrows-2)%world_size;

    // Arrays of row counts and offsets for scattering to each process
    int *send_counts = malloc(world_size * sizeof(int)); 
    int *send_offset = malloc(world_size * sizeof(int)); 
    // Arrays of counts and offsets for gathering at root from each process
    int *rcv_counts = malloc(world_size * sizeof(int));
    int *rcv_offset = malloc(world_size * sizeof(int));
    if (send_counts==NULL || send_offset==NULL || rcv_counts==NULL || rcv_offset==NULL)
    {
        return 1;
    }

    // Each process must read from the rows above and below the rows it relaxes
    const int read_only = 2*ncols;
    for (int proc=0; proc < world_size; ++proc)
    {
        /* Start at the row before where each process will start, and send the 
         * allocated number of rows plus the two read only rows (above and below) */
        send_counts[proc] = cell_split + (proc<extra_row? ncols : 0) + read_only;
        send_offset[proc] = proc==0? 0 : send_offset[proc-1] + rcv_counts[proc-1];

        // Receive the number of rows we sent minus the rows above and below
        rcv_counts[proc] = send_counts[proc] - read_only;
        // Offset starting at the row after the top boundary
        rcv_offset[proc] = send_offset[proc] + read_only/2;

        if (args.verbosity && world_rank == ROOT_PROC)
        {
            printf("Process %d doing %d rows at offset %d.\n", proc, 
                rcv_counts[proc]/ncols, rcv_offset[proc]/ncols);
        }
    }

    /* In each process allocate a buffer that will hold a read-only subset of
     * the matrix and a (smaller) buffer that will hold a write-only subset */
    double *from = malloc(sizeof(double) * send_counts[world_rank]);
    double *to   = malloc(sizeof(double) * rcv_counts[world_rank]);
    if (from==NULL || to==NULL)
    {
        return 1;
    }

    // If this process needs another iteration, it will set local_continue > 0.
    int local_continue = 0;
    /* Reducing MPI_SUM over local_continue for every process tells us whether
     * we need another iteration or not. */
    int global_continue = 0;

    /* The first and last rows are read-only and will be relaxed by another
     * process, so start at the second and finish after the penultimate. */
    const int start = ncols;
    const int end = start + rcv_counts[world_rank];
    
    for (int iters = 1; /* no exit condition - we'll break */ ;++iters)
    {
        /* Scatter chunks of the array from root to each process. Use MPI_Scatterv 
         * as the chunks will only be the same size when (nrows-2)%ntasks == 0 */
        MPI_Scatterv(global_matrix, send_counts, send_offset, MPI_DOUBLE, 
            from, send_counts[world_rank], MPI_DOUBLE, ROOT_PROC, MPI_COMM_WORLD);
        
        for (int cell=start; cell<end; ++cell)
        {
            /* If this is an edge cell, just copy the old fixed value.
             * The global index of a cell is calculated as the local
             * zero-based index plus the offset its row will be returned to. */
            if (is_edge(nrows, ncols, rcv_offset[world_rank] + cell-start))
            {
                to[cell-start] = from[cell];
            }
            // Otherwise, relax the cell and determine if it's within precision.
            else
            {
                local_continue += relax_cell(from, to, cell, ncols, precision);
            }
        }

        // Gather the chunks back into the complete array on the root.
        MPI_Gatherv(to, rcv_counts[world_rank], MPI_DOUBLE, global_matrix,
            rcv_counts, rcv_offset, MPI_DOUBLE, ROOT_PROC, MPI_COMM_WORLD);

        // Sum local_continue across every process to determine if we're done.
        MPI_Allreduce(&local_continue, &global_continue, 1, MPI_INTEGER, 
            MPI_SUM, MPI_COMM_WORLD);

        if (args.verbosity && world_rank == ROOT_PROC) 
        {
            print_matrix(global_matrix, nrows, ncols);
        }
        if (global_continue)
        {
            local_continue = 0; // If we're not done, reset for the next loop.
        }
        else 
        {
            if (world_rank == ROOT_PROC) 
            {
                printf("Reached in %d iterations.\n", iters);
            }
            break;
        }
    }

    // Tidy up and call MPI_Finalize from every process.
    free(from);
    free(to);
    if (world_rank == ROOT_PROC) 
    {
        free(global_matrix);
        free(send_counts);
        free(send_offset);
        free(rcv_counts);
        free(rcv_offset);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
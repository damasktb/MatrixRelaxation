+ Compile with: mpicc -std=gnu99 main.c -o mpi-relax
+ Run with: mpirun -np $tasks ./mpi-relax [optargs]

Where [optargs] are:

    -d (Integer, default=10000) The matrix dimensions    
    -p (Double, default=0.01) The precision threshold   
    -v (No specifier, default=false) Enable verbose mode, do not use while performance testing

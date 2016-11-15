+ Compile with: gcc -std=gnu99 -pthread main.c -o relax
+ Run with: ./relax [optargs]

Where [optargs] are:

    -d (Integer, default=100) The matrix dimensions    
    -t (Integer, default=1) The number of threads   
    -p (Double, default=0.5) The precision threshold   
    -v (No specifier, default=false) Enable verbose mode, do not use while performance testing
    

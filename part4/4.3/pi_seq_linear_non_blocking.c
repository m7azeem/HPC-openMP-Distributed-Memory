
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#define SEED     921
#define NUM_ITER 1000000000

int main(int argc, char* argv[])
{   
    int provided,size,rank;
    int count = 0;
    double x, y, z, pi;
    
    
    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status stat;
    MPI_Request requests[size - 1];
    
    double start_time, stop_time, elapsed_time;
    start_time = MPI_Wtime();

    int *allCount = calloc(size, sizeof(int));

    srand(SEED * rank); // Important: Multiply SEED by "rank" when you introduce MPI!
    
    int loc_iter = NUM_ITER/size;

    if(rank == 0)
    {
        for (int i = 1; i < size; i++) 
        {
            MPI_Irecv(&allCount[i-1], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &requests[i-1]);
        }
    }

    // Calculate PI following a Monte Carlo method
    for (int iter = 0; iter < loc_iter; iter++)
    {
        // Generate random (X,Y) points
        x = (double)random() / (double)RAND_MAX;
        y = (double)random() / (double)RAND_MAX;
        z = sqrt((x*x) + (y*y));
        
        // Check if point is in unit circle
        if (z <= 1.0)
        {
            count++;
        }
    }

    if(rank == 0)
    {
        MPI_Waitall(size - 1, requests, MPI_STATUSES_IGNORE);
    }
    else
    {
        MPI_Send(&count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }


    if(rank == 0)
    {
        for(int j = 0; j < size; j++)
        {
            count += allCount[j];
        }
    
        // Estimate Pi and display the result
        pi = ((double)count / (double)NUM_ITER) * 4.0;

        stop_time = MPI_Wtime();
	    elapsed_time = stop_time - start_time;
    	
        printf("It took %f seconds\n", elapsed_time);

        printf("The result is %f by rank: %d\n", pi, rank);
    }
    // printf("The result is %f\n", pi);
    
    MPI_Finalize();

    return 0;
}


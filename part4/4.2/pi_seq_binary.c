
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#define SEED     921
#define NUM_ITER 100//0000000


// TODO, NOT FINISHED YET. :'(
int main(int argc, char* argv[])
{   
    int provided,size,rank;
    int count, temp_count = 0;
    double x, y, z, pi;
    
    printf("log2 %f \n", log2(2));
    printf("log2 %f \n", log2(4));
    printf("log2 %f \n", log2(6));
    printf("log2 %f \n", log2(8));
    printf("log2 %f \n", log2(10));
    printf("log2 %f \n", log2(12));
    printf("log2 %f \n", log2(14));
    printf("log2 %f \n", log2(16));

    
    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status stat;
    // int *allCount = calloc(size, sizeof(int));

    srand(SEED * rank); // Important: Multiply SEED by "rank" when you introduce MPI!
    
    int loc_iter = NUM_ITER/size;

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

    if (rank != 0)
    {
        if (rank % 2 > 0) // all odd ranks
        {
            MPI_Send(&count, 1, MPI_INT, rank-1, 10, MPI_COMM_WORLD);
        } else { // all even ranks

            //get count of higher odd rank.
            MPI_Recv(&temp_count, 1, MPI_INT, rank + 1, 10, MPI_COMM_WORLD, &stat);
            count += temp_count;

            if (log2(rank) % 2 == 0) // all log2(rank) as whole numbers
            {
                // receive earlier counts 
                int next_power_of_2 = 2 ^((int)log2(rank) + 1);

                for (size_t i = rank +1; i < next_power_of_2; i++)
                {
                    //TODO: put condition to avoid it going out of size


                    /* code */
                }
                

                // send to zero
                MPI_Send(&count, 1, MPI_INT, 0, 10, MPI_COMM_WORLD);
            } else {
                int destination = 2^(int)log2(rank);
                MPI_Send(&count, 1, MPI_INT, destination, 10, MPI_COMM_WORLD);
            }

            if (sqrt(rank)/2 != 2)
            {
                MPI_Send(&count, 1, MPI_INT, rank-1, 10, MPI_COMM_WORLD);
            }
            MPI_Send(&count, 1, MPI_INT, rank-1, 10, MPI_COMM_WORLD);
        }
    }
    else if (rank == 0)
    {
        for(int j = 1; j < size; j^2)
        {
            MPI_Recv(&temp_count, 1, MPI_INT, j, 10, MPI_COMM_WORLD, &stat);
            count += temp_count;
        }
    }
    

    if(rank == 0)
    {
        for(int j = 0; j < size; j++)
        {
            // count += allCount[j];
        }
    
        // Estimate Pi and display the result
        pi = ((double)count / (double)NUM_ITER) * 4.0;

        printf("The result is %f by rank: %d\n", pi, rank);
    }
    // printf("The result is %f\n", pi);
    
    MPI_Finalize();

    return 0;
}


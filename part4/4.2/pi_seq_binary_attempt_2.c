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
      
    int rank, size, i, provided;
    
    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

    double start_time, stop_time, elapsed_time;

    start_time = MPI_Wtime();

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    int count_loc = 0;
    int number_reductions=log(size)/log(2); 

    
    double x, y, z, pi;
    int NUM_ITER_LOC=NUM_ITER/size;
    
    srand(SEED*rank); 
    
    
    for (int iter = 0; iter < NUM_ITER_LOC; iter++)
    {
        x = (double)random() / (double)RAND_MAX;
        y = (double)random() / (double)RAND_MAX;
        z = sqrt((x*x) + (y*y));
        
        if (z <= 1.0)
        {
            count_loc++;
        }
    }
    
    
    for (int j=1;j<=number_reductions;j++){
    	if (rank==0 || rank % (int)pow(2,j) ==0){
   		int count_temp=0;
     	if (j==1){
    		MPI_Recv(&count_temp, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    		} else {
    		MPI_Recv(&count_temp, 1, MPI_INT, rank+pow(2,(j-1)), 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    			}
    		
       	count_loc=count_loc+count_temp; 
    	} else {
    		if (j==1){
    		MPI_Send(&count_loc, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD);
    		} else if (rank % (int)pow(2,j-1)==0) {
    		MPI_Send(&count_loc, 1, MPI_INT, rank-pow(2,(j-1)), 0, MPI_COMM_WORLD);
    		}
    	}
    }
    if (rank==0){
    	pi = ((double)count_loc / (double)NUM_ITER) * 4.0;
    	
    	stop_time = MPI_Wtime();
    	elapsed_time = stop_time - start_time;
    	printf("The result is %f. It took %f seconds\n", pi, elapsed_time);	
    }
    MPI_Finalize();
    return 0;
}

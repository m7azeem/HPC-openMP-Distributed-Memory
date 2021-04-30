#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
	/* -------------------------------------------------------------------------------------------
		MPI Initialization 
	--------------------------------------------------------------------------------------------*/
	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
	// MPI_Init(&argc, &argv);

	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	MPI_Status stat;

	if(size != 2){
		if(rank == 0){
			printf("This program requires exactly 2 MPI ranks, but you are attempting to use %d! Exiting...\n", size);
		}
		MPI_Finalize();
		exit(0);
	}

	/* -------------------------------------------------------------------------------------------
		Loop from 8 B to 1 GB
	--------------------------------------------------------------------------------------------*/

	for(int i=0; i<=27; i++){

		long int N = 1 << i;
	
   	 	// Allocate memory for A on CPU
		double *A = (double*)malloc(N*sizeof(double));

		// Initialize all elements of A to 0.0
		for(int i=0; i<N; i++){
			A[i] = 0.0;
		}


		// int S = 4;
		// double* local_mem = (double*)malloc(S * sizeof(double));
		// double* shared_mem = (double*)malloc(S * sizeof(double));

		int S = 4;
		double* local_mem = (double*)malloc(S * sizeof(double));
		double* shared_mem = (double*)malloc(S * sizeof(double));

		for (int i = 0; i < S; i++)
    	{
			local_mem[i] = rank * S +i;
			shared_mem[i] = 0;
    	}

		// for (int j = 0; j < S; j++)
		// {
		// 	local_mem[j] = rank * S +j;
		// 	shared_mem[j] = 0;
		// }

		MPI_Win win;
		MPI_Win_create(shared_mem, S, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &win);

		



		//warmup

		for(int i=1; i<=5; i++){
			MPI_Win_fence(0, win);
			if(rank == 0){
				MPI_Put(local_mem, S, MPI_INT, 1, 1, S, MPI_INT, win);
			}
			MPI_Win_fence(0, win);
			if(rank == 0){
				MPI_Get(shared_mem, S, MPI_DOUBLE, 1, 0, S, MPI_DOUBLE, win);
			}
		}	

		MPI_Win_fence(0, win);
		
		//printf("Warmup done!-- rank=%d \n", rank);


		int loop_count = 50;

		// // Time ping-pong for loop_count iterations of data transfer size 8*N bytes
		 double start_time, stop_time, elapsed_time;
		 start_time = MPI_Wtime();
	
		//printf("rank %d \n", rank);
		for(int i=1; i<=loop_count; i++){
			MPI_Win_fence(0, win);
			if(rank == 0){
		
				
				//printf(" i=%d    putttt \n", i);
				MPI_Put(local_mem, S, MPI_INT, 1, 1, S, MPI_INT, win);
				

				//printf("i=%d   putttt done \n", i);
			}
			MPI_Win_fence(0, win);

			if(rank == 0){
				MPI_Get(shared_mem, S, MPI_DOUBLE, 1, 0, S, MPI_DOUBLE, win);
			}
			
			
			//printf("rank=%d  getttt done i=%d  \n", rank, i);
		}	

		MPI_Win_fence(0, win);
		
		//printf("done with loop, rank= %d \n", rank);
		
		

		int tag1 = 10;
		int tag2 = 20;
		

		// // Warm-up loop
		// for(int i=1; i<=5; i++){
		// 	if(rank == 0){
		// 		MPI_Send(A, N, MPI_DOUBLE, 1, tag1, MPI_COMM_WORLD);
		// 		MPI_Recv(A, N, MPI_DOUBLE, 1, tag2, MPI_COMM_WORLD, &stat);
		// 	}
		// 	else if(rank == 1){
		// 		MPI_Recv(A, N, MPI_DOUBLE, 0, tag1, MPI_COMM_WORLD, &stat);
		// 		MPI_Send(A, N, MPI_DOUBLE, 0, tag2, MPI_COMM_WORLD);
		// 	}
		// }

		
		// for(int i=1; i<=loop_count; i++){
		// 	if(rank == 0){
		// 		MPI_Send(A, N, MPI_DOUBLE, 1, tag1, MPI_COMM_WORLD);
		// 		MPI_Recv(A, N, MPI_DOUBLE, 1, tag2, MPI_COMM_WORLD, &stat);
		// 	}
		// 	else if(rank == 1){
		// 		MPI_Recv(A, N, MPI_DOUBLE, 0, tag1, MPI_COMM_WORLD, &stat);
		// 		MPI_Send(A, N, MPI_DOUBLE, 0, tag2, MPI_COMM_WORLD);
		// 	}
		// }

		stop_time = MPI_Wtime();
		elapsed_time = stop_time - start_time;

		//printf(" \n  rank=%d elaspsed time %f  \n",rank, elapsed_time);

		long int num_B = 8*N;
		long int B_in_GB = 1 << 30;
		double num_GB = (double)num_B / (double)B_in_GB;
		double avg_time_per_transfer = elapsed_time / (2.0*(double)loop_count);


		if(rank == 0) printf("%10li\t%15.9f\n", num_B, avg_time_per_transfer);

		free(A);
		MPI_Win_free(&win);
	}
	
	MPI_Finalize();

	return 0;
}

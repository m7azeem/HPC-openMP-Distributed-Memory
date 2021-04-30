#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
    int proc_id, num_procs, provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    int N = 4;
		
    double* local_mem = (double*)malloc(N * sizeof(double));
    double* shared_mem = (double*)malloc(N * sizeof(double));
    for (int i = 0; i < N; i++)
    {
        local_mem[i] = proc_id * N +i;
        shared_mem[i] = 0;
    }

    MPI_Win win;
    MPI_Win_create(shared_mem, N, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &win);

    MPI_Win_fence(0, win);
    MPI_Put(local_mem, N, MPI_INT, (proc_id + 1)%num_procs, 0, N, MPI_INT, win);
    MPI_Win_fence(0, win);

    printf("mpi rank %d got data: ", proc_id);
    for(int i=0; i<(N-1); ++i)
        printf("%f, ", shared_mem[i]);
    printf("%f \n", shared_mem[N-1]);

    MPI_Win_free(&win);
    MPI_Finalize();

    return 0;


}
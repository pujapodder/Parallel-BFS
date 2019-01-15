#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "math.h"
#include<assert.h>
#include <stdbool.h>
#define MAX_RAND 10

int *GF, *CF, *CP, *LF, *T;
char *adj,t; //adj adjacency matrix stored in each processor
int rno, cno, nor, nov, novpp, nv, d,nvrp,rowno,colno,o=0,dens=0;

void initialize()
{
	int i;
	if (colno >= rowno) {
		for (i = 0; i< nvrp*nvrp; i++) {
			int rand1 = rand() % MAX_RAND;
		//	printf("densrand=%d ", rand1);
			if (rand1<d) {
				adj[i] = 1;
				o++;
			}
		}
	}
	//for (i = 0; i<nvrp*nvrp; i++) {
	//	printf("adj= %4d \n", adj[i]);
	//}
}

int getIndex( int i, int j,  int rowSize) {
	return i*rowSize + j;
}

/*
* This method checks the frontier vector F for all 0's to indicate entire graph has traversed.
* If any un-traversed node has been found i.e if any value in F vector equals 1, then this function returns true.
* Else, function returns false.
*/
bool isClear(int F[], int rowNo, int columnNo, int size) {
	int t;
	for (t = 0; t<size; t++) {
		if (F[t] == 1) {
			return true;
		}
	}
	return false;
}
int main(int argc, char *argv[]) {
	int my_rank;
	int num_procs;
	MPI_Init(&argc, &argv);
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); //grab this process's rank
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	MPI_Comm row_comm, col_comm;        // row communication world and column communication world
	int i,m,n,k,l;
	//variables
	double t1, t2;                   // used to track execution time 
	t1 = MPI_Wtime();
	
	nor = sqrt(num_procs); //no of processor rows in 2-D partition, this is also equal to processor columns
	nov = nor; // number of input vertices
	nv = ceil(nov*1.0 / num_procs)* num_procs; // normalize the vertices count for equal adjacency matrix distribution
	d = 8; // density of the graph
	novpp = ceil(nv / pow(nor, 2));// number of vertices in each processor
	nvrp = novpp*nor;//number of vertices in row of processors
	rno = (my_rank) / nor + 1; //rowno
	cno = my_rank%nor + 1; //column no
	m = my_rank / nor;
	n = my_rank%nor;
	rowno = rno;
	colno = cno;
	//all variables for graph traversal
	GF = (int *)malloc(sizeof(int)*nv); //global frontier
	CF = (int *)malloc(sizeof(int)*novpp); //current local frontier
	CP = (int *)malloc(sizeof(int)*novpp); //local frontier
	LF = (int *)malloc(sizeof(int)*novpp); //local frontier
	T = (int *)malloc(sizeof(int)*novpp*nor); //next frontier 
	adj = (char*)malloc(sizeof(char)*nvrp*nvrp); //adjacency matrix
	//global frontier
	for (i = 0; i<nv; i++) {
		GF[i] = 0;
	}
	GF[nv] = 1;
	//local frontier
	for (i = 0; i<novpp; i++) {
		CP[i] = 0;
		CF[i] = 0;
		LF[i] = 0;
	}
	//next frontier of a row of processor
	for (i = 0; i<novpp*nor; i++) {
		T[i] = 0;
	}
	//adj adjacency matrix stored in each processor
	for (i = 0; i<nvrp*nvrp; i++) {
		adj[i] = 0;
	}

	if (my_rank == 0) {
		//printf("no of Vertices=%d \n nVertices=%d \n no of Vertices Per Processor=%d \n Ti_size=%d \n ",nov, nv, novpp, novpp*nor);
	}

	initialize();
	
	dens = 2 * o*(1.0 / (nvrp*nvrp));
	//printf("sum=%d, density=%d\n", dens, o);
	if (rowno != colno)
	{
		o = o * 2;
	}
	//printf("sum12=%d, density12=%d\n", dens, o);
	int recv_ddata;
	MPI_Reduce(&o, &recv_ddata, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	if (my_rank == 0) {
		//printf("sum=%d, density=%f\n", recv_ddata, recv_ddata*1.0 / (nov*nov));
	}
	/* Sending matrix data to processor below the diagonal, and also recieving data if a processor is below diagonal */
	//char temp;
	int chunkSize = (nvrp*nvrp)/num_procs;
	int  noofChunkstoSend = ceil(((nvrp*nvrp)*1.0) / (chunkSize*1.0));
	char* temp = adj;
	int r;
	//sending data 
	if (rowno < colno) {
		for (r = 0; r<noofChunkstoSend; r++) {
			int size = chunkSize;
			if (r == noofChunkstoSend - 1) {
				size = nvrp*nvrp - r*chunkSize;//size of sending data
			}
			MPI_Send(&temp, size, MPI_CHAR, (colno - 1)*nor + rowno - 1, 123, MPI_COMM_WORLD);
			temp += size;
		}
	}
	else if (rowno < colno) {
		for (r = 0; r<noofChunkstoSend; r++) {
			long long int size = chunkSize;
			if (r == noofChunkstoSend - 1) {
				size = nvrp*nvrp - r*chunkSize;
			}
			MPI_Recv(&temp, size, MPI_CHAR, (colno - 1)*nor + rowno - 1, 123, MPI_COMM_WORLD, &status);
			temp += size;
		}

		for (k = 0; k<nvrp; k++) {
			for (l = 0; l<nvrp; l++) {
				if (k<l) {
					t = adj[getIndex(k, l, nvrp)];
					adj[getIndex(k, l, nvrp)] = adj[getIndex(l, k, nvrp)];
					adj[getIndex(l, k, nvrp)] = t;
				}
			}
		}
	}

	// Padding zero to normalize, if input vertices count makes non-uniform distribution.
	if (n == nor - 1) {
		for (l = 0; l<nvrp; l++) {
			for (k = (nvrp - (nv - nov)); k<nvrp; k++) {
				adj[l*nvrp + k] = 0;
			}
		}
	}

	if (m == nor - 1) {
		for (l = 0;l<nvrp; l++) {
			for (k = (nvrp - (nv - nov)); k<nvrp; k++) {
				adj[k*nvrp + l] = 0;
			}
		}
	}

	//initializing the rowcomm and colcomm to communicate with the processors in 2-D distribution
	MPI_Comm_split(MPI_COMM_WORLD, colno, my_rank, &col_comm);
	MPI_Comm_split(MPI_COMM_WORLD, rowno, my_rank, &row_comm);
	//current local forntier
	for (k = 0; k< nvrp; k++) {
		CF[k] = GF[(rowno - 1)*nvrp*nor + (colno - 1)*nvrp + k];
		CP[k] = CF[k];
	}

	int* rec_buffer = (int*)malloc(sizeof(int)*nvrp);
	int* send_buffer = (int*)malloc(sizeof(int)*novpp);


	// Algorithm
	
	int roundNo = 1;
	while (isClear(GF, rowno, colno, nov)) {
		
		for (i = 0; i < novpp; i++) {
			send_buffer[i] = GF[i];
		}


		if (rowno != colno) {
			if (rowno > colno) {
				MPI_Send(send_buffer, novpp, MPI_INT, (colno - 1)*nor + rowno - 1, 123, MPI_COMM_WORLD);
				MPI_Recv(CF, novpp, MPI_INT, (colno - 1)*nor + rowno - 1, 123, MPI_COMM_WORLD, &status);
			}
			else {
				MPI_Recv(CF, novpp, MPI_INT, (colno - 1)*nor + rowno - 1, 123, MPI_COMM_WORLD, &status);
				MPI_Send(send_buffer, novpp, MPI_INT, (colno - 1)*nor + rowno - 1, 123, MPI_COMM_WORLD);
			}
		}
		MPI_Allgather(CF, novpp, MPI_INT, rec_buffer, novpp, MPI_INT, col_comm);


		// computing next frontier
		int val = 0;
		for (k = 0; k < nvrp; k++) {
			val = 0;
			for (l = 0; l < nvrp; l++) {
				val += adj[k*nvrp + l] * rec_buffer[l];
			}
			T[k] = val;
		}

		MPI_Alltoall(T, novpp, MPI_INT, rec_buffer, novpp, MPI_INT, row_comm);
		for (k = 0; k < nor; k++) {
			for (l = 0; l < nvrp; l++) {
				if (rec_buffer[k*nvrp + l] > 0) {
					T[k] = 1;
				}
			}
		}
		for (k = 0; k < novpp; k++) {
			if (CP[k] == 1 && T[k] == 1) {
				T[k] = 0;
			}
		}
		for (k = 0; k < novpp; k++) {
			if (CP[k] == o && T[k] == 1) {
				CP[k] = 1;
			}
		}

		for (k = 0; k < novpp; k++) {
			CF[k] = T[k];
		}

		MPI_Allgather(CF, novpp, MPI_INT, rec_buffer, novpp, MPI_INT, row_comm);
		MPI_Allgather(rec_buffer, nvrp, MPI_INT, CF, nvrp, MPI_INT, col_comm);

	}
	// End of an iteration.
	
	double t;
	MPI_Reduce(&t, &t2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // Computing maximum time elapsed among all processors
	t2 = MPI_Wtime();  // Marks end of the BFS.
	t1 = t2 - t1;
	if (my_rank == 0) {
		printf("time parallel bfs to traverse in %f\n", t1);
		fflush(stdout);
	}
	MPI_ERRORS_RETURN;
	MPI_Finalize();
	return EXIT_SUCCESS;
}
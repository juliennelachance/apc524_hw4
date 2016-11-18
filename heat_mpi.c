#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <time.h>

////////////////////////////////////////////////////////////////////////////////

int old_index, new_index;	// The indices which tell us which "generation" of solution we're on.

/* Dependencies. Listed in order of first usage in main. */

void do_one_iteration(int nx, int width, double dx, double dt, double heatgrid[2][nx][width], int numprocs, int myid, MPI_Status Stat){
	// Steps the heat grid forward by one time step iteration.
	int i, j;
	int source1, source2, dest1, dest2;
	double send_arr_first[nx], send_arr_last[nx];
	double recv_arr_first[nx], recv_arr_last[nx];

	if(numprocs>1){
		if(myid==0)
			source1 = dest2 = numprocs-1;
		else
			source1 = dest2 = myid-1;
		if(myid==(numprocs-1))
			dest1 = source2 = 0;
		else
			dest1 = source2 = myid+1;


		// Populate last column to send:
		for ( i = 0; i < nx; i++ ) {
			send_arr_first[i] = heatgrid[old_index][i][1];
			send_arr_last[i] = heatgrid[old_index][i][width-2];
		}
		// Send / Receive columns and update old matrix 
		MPI_Send(&send_arr_first[0], nx, MPI_DOUBLE, dest2, 1, MPI_COMM_WORLD);
		MPI_Recv(&recv_arr_last[0], nx, MPI_DOUBLE, source2, 1, MPI_COMM_WORLD, &Stat);
		MPI_Send(&send_arr_last[0], nx, MPI_DOUBLE, dest1, 2, MPI_COMM_WORLD);
		MPI_Recv(&recv_arr_first[0], nx, MPI_DOUBLE, source1, 2, MPI_COMM_WORLD, &Stat);

		// Use to replace old matrix values 
		for ( i = 0; i < nx; i++ ) {
			heatgrid[old_index][i][0] = recv_arr_first[i];
			heatgrid[old_index][i][width-1] = recv_arr_last[i];
		}
	}

		for ( i = 1; i < nx-1; i++ ) {
			// Inner points: 
			for ( j = 1; j < width-1; j++ ) {
				heatgrid[new_index][i][j] = heatgrid[old_index][i][j] + dt*( heatgrid[old_index][i-1][j] + 					heatgrid[old_index][i+1][j] + heatgrid[old_index][i][j-1] + heatgrid[old_index][i][j+1] - 
				4.0*heatgrid[old_index][i][j] )/(pow(dx,2));
			}
		}

	if(numprocs==1){
		for ( i = 1; i < nx-1; i++ ) {
		heatgrid[new_index][i][1] = heatgrid[old_index][i][1] + dt*( heatgrid[old_index][i-1][1] + heatgrid[old_index][i+1][1]
						+ heatgrid[old_index][i][width-3] + heatgrid[old_index][i][2] 
						- 4.0*heatgrid[old_index][i][1] )/(pow(dx,2));
		// Right-hand point: 
		heatgrid[new_index][i][width-2] = heatgrid[new_index][i][1];
		}

	}

}


////////////////////////////////////////////////////////////////////////////////
/* The main function. */
int main(int argc, char* argv[]) {
	time_t tstart, tend; 
	tstart = time(0);
	// This is the main file, which takes in the following terminal inputs:
	//   1) nx, to define grid of size nx^2

	int myid, numprocs;

	int i;
	int j;
	old_index=0;
	new_index=1;
	int temp;

	// Initialize MPI stuff 
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Status Stat;

	//printf("The number of processors is: %d \n ", numprocs);

	// Read input arguments from terminal command:
	int nx = atoi(argv[1]);
	
	// Constant choices:
	double kappa = 4;
	int timesteps = 100;

	// Compute mesh spacing and time step 
	double pi = M_PI;
	double dx = pi/(nx-1); // nx-1
	double t = 0.5*(pow(pi,2))/kappa;
	double dt = t/(timesteps);

	// Ensure the timestep is sufficiently small 
	while(dt > (dx*dx/(4*kappa))){
		timesteps = timesteps+10;
		dt = t/(timesteps);
	}

	int nxp = (nx/numprocs);
	int width = (nx/numprocs)+2;


	// Set up initial conditions, for each slice this time (including ghost cells)
	double heatgrid[2][nx][width];
	for (i = 0; i < nx; i++)
		for (j = 0; j < width; j++)
			heatgrid[0][i][j] = heatgrid[1][i][j] = 0.0;


	
	for (i = 1; i < (width-1); i++) {
		heatgrid[0][0][i] = heatgrid[1][0][i] = pow(cos( nxp*dx*myid + (i-1)*dx),2);
		heatgrid[0][nx-1][i] = heatgrid[1][nx-1][i] = pow(sin( nxp*dx*myid + (i-1)*dx),2);
	}


	int timecount = 0;
	for(i=0; i<timesteps; i++){
		// Do one iteration: 
		do_one_iteration(nx, width, dx, dt, heatgrid, numprocs, myid, Stat);
		// Change indices of new/old generations of solution 
		temp = old_index;
		old_index = new_index;
		new_index = temp;
	}

	/*for(i=0; i<nx; i++) {
		for(j=0; j<width; j++){
			printf("%f ", heatgrid[new_index][i][j]);
		}
		printf("\n");
	}*/

	// Calculate average temp
	double volAvg_sum_local = 0.0;
	double volAvg_sum_global = 0.0;

	for (i = 1; i < nx-1; i++) {
		for (j = 1; j < width-1; j++) {
			volAvg_sum_local += heatgrid[new_index][i][j];
		}
	}
	if(numprocs>1){
		MPI_Allreduce (&volAvg_sum_local, &volAvg_sum_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		if (myid==0) {
			printf("Volume avg temp = %f \n", volAvg_sum_global/pow((nx-2),2));
		}
	}
	else{
			printf("Volume avg temp = %f \n", volAvg_sum_local/pow((nx-2),2));
	}

  	char filename [80];
	sprintf (filename, "test_%d_%d.dat", nx, myid);

	FILE *f = fopen(filename, "wb");
	for(i=0; i<(nx-1); i++) {
		for(j=1; j<(width-1); j++){
			fprintf(f, "%f ", heatgrid[new_index][i][j]);
		}
		fprintf(f, "\n");
	}
	fclose(f);
/*
	// Externally we'll have to perform these operations on combined files 
	FILE *gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "set pm3d map \n");
	fprintf(gnuplot, "unset ytics \n");
	fprintf(gnuplot, "unset xtics \n");
	fprintf(gnuplot, "set xlabel 'x' \n");
	fprintf(gnuplot, "set ylabel 'y' \n");
	fprintf(gnuplot, "set title 'Heat Diffusion Equation; nx = %d' \n", nx);
	fprintf(gnuplot, "splot 'test_%d_0.dat' matrix w image\n", nx);
	fflush(gnuplot);
*/


	tend = time(0); 
	int runtime = difftime(tend, tstart);
	printf("Time to run = %d s \n", runtime);


	MPI_Finalize();

	return(0);
}
////////////////////////////////////////////////////////////////////////////////

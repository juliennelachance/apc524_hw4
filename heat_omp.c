#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <time.h>

////////////////////////////////////////////////////////////////////////////////

int nthreads;
int old_index, new_index;	// The indices which tell us which "generation" of solution we're on.

/* Dependencies. Listed in order of first usage in main. */
void usage(void) {
	// Offers the user some guidance on terminal inputs:
	printf("Usage: This function takes two arguments: \n"
		"1) The input nx for solution with grid size nx^2 (>type int). \n"
		"2) The input nthreads for number of parallel threads to use (>type int) \n"
		"Example terminal command: ./heat_omp 16 4 \n"
	);
}

void do_one_iteration(int nx, double dx, double dt, double heatgrid[2][nx][nx]){
	// Steps the heat grid forward by one time step iteration.
	int i, j;
	#pragma omp parallel for shared(heatgrid, nx, dx, dt) private(i,j)
	for ( i = 1; i < nx-1; i++ ) {
		// Inner points: 
		for ( j = 1; j < nx-1; j++ ) {
			heatgrid[new_index][i][j] = heatgrid[old_index][i][j] + dt*( heatgrid[old_index][i-1][j] + 						heatgrid[old_index][i+1][j] + heatgrid[old_index][i][j-1] + heatgrid[old_index][i][j+1] - 
				4.0*heatgrid[old_index][i][j] )/(pow(dx,2));
		}
		// Left-hand point:
		heatgrid[new_index][i][0] = heatgrid[old_index][i][0] + dt*( heatgrid[old_index][i-1][0] + heatgrid[old_index][i+1][0]
						+ heatgrid[old_index][i][nx-2] + heatgrid[old_index][i][1] 
						- 4.0*heatgrid[old_index][i][0] )/(pow(dx,2));
		// Right-hand point: 
		heatgrid[new_index][i][nx-1] = heatgrid[new_index][i][0];
	}
}


////////////////////////////////////////////////////////////////////////////////
/* The main function. */
int main(int argc, char* argv[]) {
	time_t tstart, tend; 
	tstart = time(0);
	// This is the main file, which takes in the following terminal inputs:
	//   1) nx, to define grid of size nx^2
	//   2) nthreads, the number of threads for parallel runs 
	int i;
	int j;
	old_index=0;
	new_index=1;
	int temp;

	// Display usage info if number of arguments incorrect:
	if (argc != 3) {
		usage(); 
		return 0;
	}
	// Read input arguments from terminal command:
	int nx = atoi(argv[1]);
	nthreads = atoi(argv[2]);

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

	// Set environment conditions
	omp_set_dynamic(0);     // Explicitly disable dynamic teams
	omp_set_num_threads(nthreads);

	// Set up initial conditions
	double heatgrid[2][nx][nx];

	#pragma omp parallel for shared(heatgrid) private(i,j)
	for (i = 0; i < nx; i++)
		for (j = 0; j < nx; j++)
			heatgrid[0][i][j] = heatgrid[1][i][j] = 0.0;
		
	#pragma omp parallel for shared(heatgrid) private(i)
	for (i = 0; i < nx; i++) {
		heatgrid[0][0][i] = heatgrid[1][0][i] = pow(cos(i*dx),2);
		heatgrid[0][nx-1][i] = heatgrid[1][nx-1][i] = pow(sin(i*dx),2);
	}

	/*for(i=0; i<nx; i++) {
		for(j=0; j<nx; j++){
			printf("%f ", heatgrid[new_index][i][j]);
		}
		printf("\n");
	}*/


	int timecount = 0;
	while (timecount < timesteps){
		// Do one iteration: 
		do_one_iteration(nx, dx, dt, heatgrid);
		// Change indices of new/old generations of solution 
		temp = old_index;
		old_index = new_index;
		new_index = temp;
		timecount = timecount+1;
	}

	temp = old_index;
	old_index = new_index;
	new_index = temp;

	

	
  	char filename1 [80];
	sprintf (filename1, "omp_heatgrid_nx_%d_nthreads_%d.dat", nx, nthreads);
	FILE *f = fopen(filename1, "wb");
	for(i=0; i<(nx); i++) {
		for(j=0; j<(nx); j++){
			fprintf(f, "%f ", heatgrid[new_index][i][j]);
		}
		fprintf(f, "\n");
	}
	fclose(f);

	// Plotting: 
/*
	FILE *gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "set pm3d map \n");
	fprintf(gnuplot, "unset ytics \n");
	fprintf(gnuplot, "unset xtics \n");
	fprintf(gnuplot, "set xlabel 'x' \n");
	fprintf(gnuplot, "set ylabel 'y' \n");
	fprintf(gnuplot, "set title 'Heat Diffusion Equation; nx = %d' \n", nx);
	fprintf(gnuplot, "splot 'test.dat' matrix w image\n");
	fflush(gnuplot);
*/

	double volAvg = 0.0;
	for(i=1; i<nx-1; i++) {
		for(j=1; j<nx-1; j++){
			volAvg = volAvg + heatgrid[new_index][i][j];
		}
	}
	volAvg = volAvg / pow(nx-2,2);
	//printf("Volume avg temp = %f \n", volAvg);

	tend = time(0); 
	int runtime = difftime(tend, tstart);
	//printf("Time to run = %d s \n", runtime);

  	char filename2 [80];
	sprintf (filename2, "volavg_time_omp_nx_%d_%d.dat", nx, nthreads);
	FILE *f2 = fopen(filename2, "wb");
		fprintf(f2, "Vol avg: %f \n", volAvg);
		fprintf(f2, "Runtime: %d s \n", runtime);
	fclose(f2);

	return(0);
}
////////////////////////////////////////////////////////////////////////////////




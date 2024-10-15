/**
 * Author: Aman Hogan-Bailey
 * Solves 2d heat equation using parralelization,
 * correct loop ordering, and different types
 * of scheduleing for loop iterations
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>

#define PI 3.14159265358979323846
#define GRID_SIZE 300

// Column major ordering
#define I2D(i, j, N) ((i) + (j) * (N + 1))

void initialize(double *u, int N, double dx, double dy);
void time_step(double *u, double *u_new, int N, double dt, double dx2, double dy2, char *schedule);
double compute_error(double *u, int N, double T, double dx, double dy);

int main() 
{
    double T = 1.0; // Final time
    double CFL = 0.25; // CFL confition
    double dx = PI / GRID_SIZE; // veclocity of space in x
    double dy = PI / GRID_SIZE; // veclocity of space in y
    double dx2 = dx * dx; // acceleration of space ein x
    double dy2 = dy * dy; // acceleration of space in y
    double dt = CFL * dx2; // change in time

    // Ensure time stepping does not overshoot T
    int nsteps = (int)ceil(T / dt);  
    dt = T / nsteps;

    double *u = (double *)malloc((GRID_SIZE + 1) * (GRID_SIZE + 1) * sizeof(double));
    double *u_new = (double *)malloc((GRID_SIZE + 1) * (GRID_SIZE + 1) * sizeof(double));

    // Initialize u with the initial condition
    initialize(u, GRID_SIZE, dx, dy);

    // Run with 24 threads and different schedules
    char *schedules[] = {"static", "dynamic", "guided"};
    int num_threads = 24;
    omp_set_num_threads(num_threads);  // Set 24 threads


    FILE *csv_file = fopen("heat_solver.csv", "a");
    if (csv_file == NULL) 
    {
        printf("Error opening file!\n");
    }

    fprintf(csv_file, "Scheduling,Time\n");

    for (int s = 0; s < 3; s++) 
    {
        printf("Running with %d threads and %s scheduling...\n", num_threads, schedules[s]);
        double start_time = omp_get_wtime();

        // For a given number of time steps, calcualute heat eq. update ...
        double time = 0.0;
        for (int t_step = 0; t_step < nsteps; t_step++)
        {
            time_step(u, u_new, GRID_SIZE, dt, dx2, dy2, schedules[s]);

            double *tmp = u;
            u = u_new;
            u_new = tmp;
            time += dt;
        }

        double end_time = omp_get_wtime();
        printf("Scheduling = %s, Time = %f seconds\n", schedules[s], end_time - start_time);
        fprintf(csv_file, "%s, %lf\n", schedules[s], end_time - start_time);
    }

    free(u);
    free(u_new);
    return 0;
}


/**
 * initializes grid using initial condition:
 * u(x,y,t=0) = sin(x) sin(y) using optimized ordering
 * @param u 2d grid array
 * @param N dim of grid array
 * @param dx velocity in space in x
 * @param dy velocity in space in y
 */
void initialize(double *u, int N, double dx, double dy)
{
    for (int j = 0; j <= N; j++) 
    {
        for (int i = 0; i <= N; i++) 
        {
            double x = i * dx;
            double y = j * dy;
            u[I2D(i, j, N)] = sin(x) * sin(y);
        }
    }
}


/**
 * applies heat equation update at each time step for 2d u using parallelization and
 * optimized column ordering and scheduling of loop iterations
 * @param u 2d grid array
 * @param N dim of grid array
 * @param u_new new u value
 * @param dt change in time
 * @param dx2 accelration in space in x
 * @param dy2 accelration in space in y
 * @param schedule type of schedule: satic, dynamic, and guided
 */
void time_step(double *u, double *u_new, int N, double dt, double dx2, double dy2, char *schedule)
{
    if (strcmp(schedule, "static") == 0)
    {
        #pragma omp parallel for schedule(static)
        for (int j = 1; j < N; j++) 
        {
            for (int i = 1; i < N; i++) 
            {
                u_new[I2D(i, j, N)] = u[I2D(i, j, N)] 
                    + dt/dx2 * (u[I2D(i+1, j, N)] - 2*u[I2D(i, j, N)] + u[I2D(i-1, j, N)])
                    + dt/dy2 * (u[I2D(i, j+1, N)] - 2*u[I2D(i, j, N)] + u[I2D(i, j-1, N)]);
            }
        }
    }

    else if (strcmp(schedule, "dynamic") == 0) 
    {
        #pragma omp parallel for schedule(dynamic)
        for (int j = 1; j < N; j++) 
        {
            for (int i = 1; i < N; i++) 
            {
                u_new[I2D(i, j, N)] = u[I2D(i, j, N)] 
                    + dt/dx2 * (u[I2D(i+1, j, N)] - 2*u[I2D(i, j, N)] + u[I2D(i-1, j, N)])
                    + dt/dy2 * (u[I2D(i, j+1, N)] - 2*u[I2D(i, j, N)] + u[I2D(i, j-1, N)]);
            }
        }
    } 
    
    else if (strcmp(schedule, "guided") == 0) 
    {
        #pragma omp parallel for schedule(guided)
        for (int j = 1; j < N; j++) 
        {
            for (int i = 1; i < N; i++) 
            {
                u_new[I2D(i, j, N)] = u[I2D(i, j, N)] 
                    + dt/dx2 * (u[I2D(i+1, j, N)] - 2*u[I2D(i, j, N)] + u[I2D(i-1, j, N)])
                    + dt/dy2 * (u[I2D(i, j+1, N)] - 2*u[I2D(i, j, N)] + u[I2D(i, j-1, N)]);
            }
        }
    }
}


/**
 * calcs error between numeric solution and exact solution 
 * using optimized loop ordering
 * @param u 2d grid array
 * @param T time interval
 * @param N dim of grid array
 * @param dx velocity in space in x
 * @param dy velocity in space in y
 */
double compute_error(double *u, int N, double T, double dx, double dy)
 {
    double diff = 0.0;
    for (int j = 0; j <= N; j++) 
    {
        for (int i = 0; i <= N; i++) 
        {
            double x = i * dx;
            double y = j * dy;
            double exact = exp(-2*T) * sin(x) * sin(y);
            double error = u[I2D(i, j, N)] - exact;
            diff += error * error * dx * dy;
        }
    }
    return sqrt(diff);
}
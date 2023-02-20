#include "common.h"
#include <mpi.h>
#include <iostream>
#include <cmath>

using namespace std;

// Put any static global variables here that you will use throughout the simulation.
int NUM_BLOCKS;
int tot_num_bins;
int NUM_PROC_X, NUM_PROC_Y, del_X, del_Y;



// Apply the force from neighbor to particle
void apply_force(particle_t& particle, particle_t& neighbor) {
    // Calculate Distance
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;

    // Check if the two particles should interact
    if (r2 > cutoff * cutoff)
        return;

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);

    // Very simple short-range repulsive force
    double coef = (1 - cutoff / r) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}

// Integrate the ODE
void move(particle_t& p, double size) {
    // Slightly simplified Velocity Verlet integration
    // Conserves energy better than explicit Euler method
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x += p.vx * dt;
    p.y += p.vy * dt;

    // Bounce from walls
    while (p.x < 0 || p.x > size) {
        p.x = p.x < 0 ? -p.x : 2 * size - p.x;
        p.vx = -p.vx;
    }

    while (p.y < 0 || p.y > size) {
        p.y = p.y < 0 ? -p.y : 2 * size - p.y;
        p.vy = -p.vy;
    }
}


void init_simulation(particle_t* parts, int num_parts, double size, int rank, int num_procs) {
	// You can use this space to initialize data objects that you may need
	// This function will be called once before the algorithm begins
	// Do not do any particle simulation here


    // The input is a copy of the entire set of the particles.
    // Each process fetches its local set in this function.

    NUM_BLOCKS =  size/cutoff;
    tot_num_bins = (NUM_BLOCKS+2)*(NUM_BLOCKS+2);

    // These variables set up the grid simulation
    /* 
    The simulation has NUM_BLOCKS*NUM_BLOCKS cells. These cells are divided into a grid of dimension
    NUM_PROC_X*NUM_PROC_Y = num_procs. 
    */
	
	NUM_PROC_X = sqrt(num_procs); // number of processor divisions in x
	NUM_PROC_Y = num_procs / NUM_PROC_X; // number of processor divisions in y
	del_X = (NUM_BLOCKS / NUM_PROC_X) + 1; // number of bins across each processor division in x
	del_Y = (NUM_BLOCKS / NUM_PROC_Y) + 1; // number of bins across each processor division in y
    int bins_per_proc = del_X*del_Y;

	// Now get the grid index for each processor
	int proc_X = rank % NUM_PROC_X; 
	int proc_Y = rank / NUM_PROC_X;
    /* 
    This means that the processor covers the chunks of bins defined as:
        bin_x \in [proc_X*del_X, (proc_X + 1)*del_X]
        bin_y \in [proc_Y*del_Y, (proc_Y + 1)*del_Y]
    */
    

    // Create the arrays for particle storage (also defined for each processor)
    // By creating arrays of size tot_num_bins (as opposed to bins_per_proc), we can use the same index scheme for particles across all processors
    // (But we might need to revisit this for memory efficiency)

    int* bins = (int*) malloc(tot_num_bins * sizeof(int)); //  represents the bins (w padding) that holds heads of linked lists
    int* part_links = (int*) malloc(tot_num_bins * sizeof(int)); // indexes the particles in the bins in the bins array
    // Instantiate the arrays with the correct starting vals
    for (int i = 0; i < tot_num_bins; i++) {
        bins[i] = -1;
    }
    for (int i = 0; i < num_parts; i++) {
        part_links[i] = -1;
    }

    // Now each processor goes through all the particles
    // If the particle is in its bin, then store it
    for (int i = 0; i < num_parts; ++i) {
    // Get what row and column the particle would be in, with padding
        int bins_x = (parts[i].x * NUM_BLOCKS / size) + 1;
        int bins_y = (parts[i].y * NUM_BLOCKS / size) + 1;
        
        // Store the particle if it's in the processor's array
        if (((bins_x/del_X) == proc_X) & ((bins_y/del_Y) == proc_Y)) {

            int bin_id = bins_x + (NUM_BLOCKS+2)*bins_y;
            part_links[i] = bins[bin_id];
            bins[bin_id] = i;
            parts[i].ax = parts[i].ay = 0;

        }
    }

	
}

void simulate_one_step(particle_t* parts, int num_parts, double size, int rank, int num_procs) {

    // Each process runs simulation for its local using the O(n) serial algorithm.
    // Redistribute the outgoing particles

    // CURRENT ISSUE: how to get part_links and bins defined in this scope??
    

    // loop over bins to compute forces, skipping the padding bins
    for (int i = 0; i < NUM_BLOCKS; i++) {
        for (int j = 0; j < NUM_BLOCKS; j++) {
            
            int bin_id_pad = (i+1) + (NUM_BLOCKS+2)*(j+1);
            int part_1_id = bins[bin_id_pad];
         
            while (part_1_id >= 0) {
                for (int m = -1; m <= 1; m++) {
                    for (int n = -1; n <=1; n++) {

                        int part_2_id = bins[bin_id_pad + n + (NUM_BLOCKS+2)*m];
                        
                        while (part_2_id >= 0) {
                            apply_force(parts[part_1_id], parts[part_2_id]);
                            part_2_id = part_links[part_2_id];
                        }
                    }
                }
                part_1_id = part_links[part_1_id];
            }
        }
    }
    
    // Mmve Particles
    for (int i = 0; i < num_parts; ++i) {
        move(parts[i], size);
    }

}

    

void gather_for_save(particle_t* parts, int num_parts, double size, int rank, int num_procs) {
    // Write this function such that at the end of it, the master (rank == 0)
    // processor has an in-order view of all particles. That is, the array
    // parts is complete and sorted by particle id.
}

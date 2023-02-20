#include "common.h"
#include <mpi.h>

// Put any static global variables here that you will use throughout the simulation.

void init_simulation(particle_t* parts, int num_parts, double size, int rank, int num_procs) {
	// You can use this space to initialize data objects that you may need
	// This function will be called once before the algorithm begins
	// Do not do any particle simulation here

    // BIN THE PARTICLES

    // Instantiate the bins
    NUM_BLOCKS = size/cutoff;
    tot_num_bins = (NUM_BLOCKS+2)*(NUM_BLOCKS+2);

    bins = (int*) malloc(tot_num_bins * sizeof(int)); // holds id of first particle in the linked list
    part_links = (int*) malloc(num_parts * sizeof(int)); // indexes the particles in the bins in the bins array
    
    for (int i = 0; i < tot_num_bins; i++) {
        bins[i] = -1;
    }

    for (int i = 0; i < num_parts; i++) {
        part_links[i] = -1;
        parts[i].ax = parts[i].ay = 0;
    }

    // Fill the bins

    for (int i = 0; i < num_parts; ++i) {
        // Get what row and column the particle would be in, with padding
        int dx = (parts[i].x * NUM_BLOCKS / size) + 1;
        int dy = (parts[i].y * NUM_BLOCKS / size) + 1;
        int bin_id = dx + (NUM_BLOCKS+2)*dy;

			// -1 if bin holds -1, the bin's particle id otherwise
			part_links[i] = bins[bin_id];
			bins[bin_id] = i;
    }


}

void simulate_one_step(particle_t* parts, int num_parts, double size, int rank, int num_procs) {
    // Write this function
}

void gather_for_save(particle_t* parts, int num_parts, double size, int rank, int num_procs) {
    // Write this function such that at the end of it, the master (rank == 0)
    // processor has an in-order view of all particles. That is, the array
    // parts is complete and sorted by particle id.
}

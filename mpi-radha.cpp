#include "common.h"
#include <mpi.h>
#include <iostream>
#include <cmath>

using namespace std;

// Put any static global variables here that you will use throughout the simulation.
int NUM_BLOCKS;
int tot_num_bins;

// These aren't static, but they need to be globa.///
int* bins;
int *part_links;


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

    NUM_BLOCKS =  size/cutoff;
    tot_num_bins = (NUM_BLOCKS+2)*(NUM_BLOCKS+2);
   
    bins = (int*) malloc(tot_num_bins * sizeof(int)); //  represents the bins (w padding) that holds heads of linked lists
    part_links = (int*) malloc(num_parts * sizeof(int)); // indexes the particles in the bins in the bins array

    for (int i = 0; i < tot_num_bins; i++) {
        bins[i] = -1;
    }

    // declare array that will 
    for (int i = 0; i < num_parts; i++) {
        part_links[i] = -1;
    }



}

void simulate_one_step(particle_t* parts, int num_parts, double size, int rank, int num_procs) {

    // assign particles to bins
    for (int i = 0; i < num_parts; ++i) {
        // Get what row and column the particle would be in, with padding
        int dx = (parts[i].x * NUM_BLOCKS / size) + 1;
        int dy = (parts[i].y * NUM_BLOCKS / size) + 1;
        int bin_id = dx + (NUM_BLOCKS+2)*dy;

        // eliminated if-statement checks here
        part_links[i] = bins[bin_id];
        bins[bin_id] = i;

        parts[i].ax = parts[i].ay = 0;
    }

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

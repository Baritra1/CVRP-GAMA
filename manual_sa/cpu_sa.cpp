// Copyright 2018 D-Wave Systems Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// ===========================================================================

#include <cstdint>
#include <math.h>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <omp.h>
#include "cpu_sa.h"

#define FASTRAND(rand) do {                       \
    uint64_t x = rng_state[0];                    \
    uint64_t const y = rng_state[1];              \
    rng_state[0] = y;                             \
    x ^= x << 23;                                 \
    rng_state[1] = x ^ y ^ (x >> 17) ^ (y >> 26); \
    rand = rng_state[1] + y;                      \
} while (0)

#define RANDMAX ((uint64_t)-1L)

using namespace std;

// thread-local RNG state
thread_local uint64_t rng_state[2];

double get_flip_energy(
    int var,
    std::int8_t *state,
    const vector<double>& h,
    const vector<int>& degrees,
    const vector<vector<int>>& neighbors,
    const vector<vector<double>>& neighbour_couplings
) {
    double energy = h[var];
    for (int n_i = 0; n_i < degrees[var]; n_i++) {
        energy += state[neighbors[var][n_i]] * neighbour_couplings[var][n_i];
    }
    return -2 * state[var] * energy;
}

template <VariableOrder varorder, Proposal proposal_acceptance_criteria>
void simulated_annealing_run(
    std::int8_t* state,
    const vector<double>& h,
    const vector<int>& degrees,
    const vector<vector<int>>& neighbors,
    const vector<vector<double>>& neighbour_couplings,
    const int sweeps_per_beta,
    const vector<double>& beta_schedule
) {
    const int num_vars = h.size();
    double *delta_energy = (double*)malloc(num_vars * sizeof(double));

    uint64_t rand;
    for (int var = 0; var < num_vars; var++) {
        delta_energy[var] = get_flip_energy(var, state, h, degrees,
                                            neighbors, neighbour_couplings);
    }

    bool flip_spin;
    for (int beta_idx = 0; beta_idx < (int)beta_schedule.size(); beta_idx++) {
        const double beta = beta_schedule[beta_idx];
        for (int sweep = 0; sweep < sweeps_per_beta; sweep++) {
            const double threshold = 44.36142 / beta;
            for (int varI = 0; varI < num_vars; varI++) {
                int var;
                if constexpr (varorder == Random) {
                    FASTRAND(rand);
                    var = rand % num_vars;
                } else {
                    var = varI;
                }
                if (delta_energy[var] >= threshold) continue;

                flip_spin = false;
                if constexpr (proposal_acceptance_criteria == Metropolis) {
                    if (delta_energy[var] <= 0.0) {
                        flip_spin = true;
                    } else {
                        FASTRAND(rand);
                        if (exp(-delta_energy[var]*beta) * RANDMAX > rand) {
                            flip_spin = true;
                        }
                    }
                } else {
                    FASTRAND(rand);
                    if (RANDMAX > rand * (1+exp(delta_energy[var]*beta))) {
                        flip_spin = true;
                    }
                }

                if (flip_spin) {
                    const std::int8_t multiplier = 4 * state[var];
                    for (int n_i = 0; n_i < degrees[var]; n_i++) {
                        int neighbor = neighbors[var][n_i];
                        delta_energy[neighbor] += multiplier *
                            neighbour_couplings[var][n_i] * state[neighbor];
                    }
                    state[var] *= -1;
                    delta_energy[var] *= -1;
                }
            }
        }
    }

    free(delta_energy);
}

double get_state_energy(
    std::int8_t* state,
    const vector<double>& h,
    const vector<int>& coupler_starts,
    const vector<int>& coupler_ends,
    const vector<double>& coupler_weights
) {
    double energy = 0.0;
    for (unsigned int var = 0; var < h.size(); var++) {
        energy += state[var] * h[var];
    }
    for (unsigned int c = 0; c < coupler_starts.size(); c++) {
        energy += state[coupler_starts[c]] * coupler_weights[c] *
                  state[coupler_ends[c]];
    }
    return energy;
}

int general_simulated_annealing(
    std::int8_t* states,
    double* energies,
    const int num_samples,
    const vector<double> h,
    const vector<int> coupler_starts,
    const vector<int> coupler_ends,
    const vector<double> coupler_weights,
    const int sweeps_per_beta,
    const vector<double> beta_schedule,
    const uint64_t seed,
    const VariableOrder varorder,
    const Proposal proposal_acceptance_criteria,
    callback interrupt_callback,
    void * const interrupt_function
) {
    const int num_vars = h.size();
    if (!((coupler_starts.size() == coupler_ends.size()) &&
                (coupler_starts.size() == coupler_weights.size()))) {
        throw runtime_error("coupler vectors have mismatched lengths");
    }

    // Build graph structure
    vector<int> degrees(num_vars, 0);
    vector<vector<int>> neighbors(num_vars);
    vector<vector<double>> neighbour_couplings(num_vars);

    for (unsigned int cplr = 0; cplr < coupler_starts.size(); cplr++) {
        int u = coupler_starts[cplr];
        int v = coupler_ends[cplr];
        if ((u < 0) || (v < 0) || (u >= num_vars) || (v >= num_vars)) {
            throw runtime_error("coupler indexes contain an invalid variable");
        }
        neighbors[u].push_back(v);
        neighbors[v].push_back(u);
        neighbour_couplings[u].push_back(coupler_weights[cplr]);
        neighbour_couplings[v].push_back(coupler_weights[cplr]);
        degrees[u]++;
        degrees[v]++;
    }

    cout << "Using threads " << omp_get_max_threads() << "\n";

    #pragma omp parallel for schedule(static) default(none) \
    shared(states, energies, h, degrees, neighbors, neighbour_couplings, \
           coupler_starts, coupler_ends, coupler_weights, sweeps_per_beta, \
           beta_schedule, seed, varorder, proposal_acceptance_criteria, \
           interrupt_function, interrupt_callback, num_samples, num_vars)
    for (int sample = 0; sample < num_samples; sample++) {
        // Initialize RNG for each thread/sample
        uint64_t thread_seed = seed + sample + 123456789ULL * omp_get_thread_num();
        rng_state[0] = thread_seed ? thread_seed : 1;
        rng_state[1] = thread_seed ^ 123456789ULL;
        
        std::int8_t *state = states + sample*num_vars;
        if (varorder == Random) {
            if (proposal_acceptance_criteria == Metropolis) {
                simulated_annealing_run<Random, Metropolis>(state, h, degrees,
                                                    neighbors, neighbour_couplings,
                                                    sweeps_per_beta, beta_schedule);
            } else {
                simulated_annealing_run<Random, Gibbs>(state, h, degrees,
                                                     neighbors, neighbour_couplings,
                                                     sweeps_per_beta, beta_schedule);
            }
        } else {
            if (proposal_acceptance_criteria == Metropolis) {
                simulated_annealing_run<Sequential, Metropolis>(state, h, degrees,
                                                     neighbors, neighbour_couplings,
                                                     sweeps_per_beta, beta_schedule);
            } else {
                simulated_annealing_run<Sequential, Gibbs>(state, h, degrees,
                                                      neighbors, neighbour_couplings,
                                                      sweeps_per_beta, beta_schedule);
            }
        }
        energies[sample] = get_state_energy(state, h, coupler_starts,
                                            coupler_ends, coupler_weights);
    }

    return num_samples;
}

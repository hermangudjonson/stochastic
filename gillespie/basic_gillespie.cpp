
/*
 * c++11 implementation of basic gillespie stochastic simulation algorithm
 *
 * General scheme for simulating reactions with 
 * arbitrary propensities and stoichiometries
 */

#include "basic_gillespie.h"

// main Gillespie loop
void BasicGillespie::run(double TMAX,unsigned int RMAX) {
  // initialization
  double current_time = 0.0;
  // initialize species
  init_species();
  // initialize reaction counts
  rxn_count.assign(reactions.size(),0);
  // vector to hold cumulative propensities for each iteration
  std::vector<double> propensities(reactions.size(),0.0);
  // pause times iterator
  auto pause_it = pause_times.begin();
  
  while (current_time < TMAX && num_rxns < RMAX) {
    double total_propensity = 0.0;
    // calculate all reaction propensities
    for (std::size_t i = 0; i < reactions.size(); ++i) {
      double rprop = reactions[i]->propensity(current_state,current_time);      
      propensities[i] = rprop + total_propensity;
      total_propensity += rprop;      
    }

    // get time of next reaction
    // rescaled unit exponential to total_propensity
    double next_rxn_time = exp(gen) / total_propensity;

    // for event based recording, tells us whether we need a new record
    if (tkeeper.new_entry(times,current_time)) {
      concentrations.push_back(std::vector<int>(species.size(),0));
    }
    
    // update time
    double next_time = current_time + next_rxn_time;
    if (!(next_time > current_time)) {
      std::cout << "Timestep too small! Time not progressing!!" << std::endl;
    }

    // calculating condition for pausing
    bool pause = false;
    if (pause_it != pause_times.end()) {
      if (*pause_it < next_time) {
	// need to pause time
	pause = true;
      }
    }

    if (pause) {
      current_time = nextafter(*pause_it,next_time); // time just after pause time
      ++pause_it; // update pause time iterator
    }
    else {
      // no pause, proceeding with reaction
      current_time = next_time;
    
      // step time and record until we pass event
      while (tkeeper.step_time(times,current_time)) {
	record_state(times[tkeeper.tindex-1],concentrations[tkeeper.tindex-1]);
      }

      // choose reaction to perform
      double target_prop = uniform(gen) * total_propensity;
      size_t chosen_rxn = 0;
      while(target_prop >= propensities[chosen_rxn]) {
	++chosen_rxn;
      }

      // only perform if a valid reaction is chosen
      // if propensities are all 0, this will be out of range
      if (chosen_rxn < reactions.size()) {
	// perform reaction
	reactions[chosen_rxn]->perform(current_state);
	++rxn_count[chosen_rxn];
	++num_rxns;
	// count is one step messed up for event recording
	// count also displays total performed not just recorded for set times
      }
    }
  }
}



/*
 * c++11 implementation of basic gillespie stochastic simulation algorithm
 *
 * General scheme for simulating reactions with 
 * arbitrary propensities and stoichiometries
 */

#ifndef BASIC_GILLESPIE_H
#define BASIC_GILLESPIE_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <fstream>
#include <random>
#include <unordered_map>
#include <vector>
#include "timekeeper.h"

// State represents current species composition during simulation
typedef std::unordered_map<std::string,int> State;
// Params represent parameters associated with a particular reaction
typedef std::unordered_map<std::string,double> Params;

/*
 * Class to contain reaction objects
 * (hopefully makes reaction specification easier and more concise)
 *
 * 4 things specify a reaction: propensity, result, name, parameters
 * Default reaction can be created by specifying these things directly
 */
class Reaction {
public:
  Reaction() {}
  Reaction(std::function<double(State&,double)> pr,std::function<void(State&)> pe,
	   std::string n,Params par) : prop(pr),perf(pe),name(n),parameters(par) {}
  virtual double propensity(State &current,double t) {return prop(current,t);}
  virtual void perform(State &current) {return perf(current);}
  std::string name;
  Params parameters;
private:
  // simply going to store these functions 
  std::function<double(State&,double)> prop;
  std::function<void(State&)> perf;
};

/*
 * Class to organize reactions and handle simulation
 *
 * Current simulation state is represented by an unordered map of
 * species:copy_number pairs
 */
class BasicGillespie {
public:
  // initialize process with initial conditions and times
  BasicGillespie(State IC,std::vector<double> &ts,double PRC=1e-15) : current_state(IC), tkeeper(times,ts,PRC) {tkeeper.init_times(times,0.0);}
  // initialize process without set times for recording
  BasicGillespie(State IC,double PRC=1e-15) : current_state(IC), tkeeper(times,PRC) {tkeeper.init_times(times,0.0);}

  void add_reaction(std::shared_ptr<Reaction> r) {reactions.push_back(r);}
  void run(double TMAX = std::numeric_limits<double>::max(),
	   unsigned int RMAX = std::numeric_limits<unsigned int>::max());
  // insert pause for when propensities depend on time
  void add_pause(double t)
  {
    pause_times.insert(std::lower_bound(pause_times.begin(),pause_times.end(),t),t);
  }
  
  void print()
  {
    for (auto t : times){std::cout << t << ' ';}
    std::cout << std::endl;
    for (auto sp : species){std::cout << sp << ' ';}
    std::cout << std::endl;
    concentrations_out(std::cout,' ');
    for (auto r : reactions){std::cout << r->name << ' ';}
    std::cout << std::endl;
    rxn_count_out(std::cout,' ');
    params_out(std::cout,' ');
  }
  void write(std::string filename,bool include_header=true,bool append=false)
  {
    std::ofstream to_file;
    if (append) {
      to_file.open(filename,std::ios::out | std::ios::app);
    }
    else {
      to_file.open(filename,std::ios::out);
    }
    if (to_file.is_open()) {
      if (include_header) {
	for (auto t : times){to_file << t << '\t';}
	to_file << std::endl;
	for (auto sp : species){to_file << sp << '\t';}
	to_file << std::endl;
	for (auto r : reactions){to_file << r->name << '\t';}
	to_file << std::endl;
	params_out(to_file,'\t');
      }
      concentrations_out(to_file,'\t');
      rxn_count_out(to_file,'\t');
    }
  }
  void concentrations_out(std::ostream &os,char sep='\t')
  {
    for (std::size_t i = 0; i < times.size(); ++i) {
      for (auto c : concentrations[i]){os << c << sep;}
      os << std::endl;
    }
  }
  void rxn_count_out(std::ostream &os,char sep='\t')
  {
    for (auto rc : rxn_count){os << rc << sep;}
    os << std::endl;
  }
  void params_out(std::ostream &os,char sep='\t')
  {
    for (auto r : reactions) {
      os << '{';
      for (auto p : r->parameters){os << p.first << ':' << p.second << ',';}
      os << '}';
      os << sep;
    }
    os << std::endl;
  }
  void set_seed(int s) {gen.seed(s);} // for setting seed manually
  
private:
  void init_species()
  {
    // initialize species list for record keeping
    for (auto sp : current_state) {
      species.push_back(sp.first);
    }
    std::sort(species.begin(),species.end()); // sorted order for consistency
    
    // initialize concentration records
    concentrations = std::vector< std::vector<int> >(times.size(),
						     std::vector<int>(species.size(),0));
  }
  void record_state(double time,std::vector<int> &dest)
  {
    // iterate through species and record current values
    for (std::size_t i = 0; i < species.size(); ++i) {
      dest[i] = current_state[species[i]]; // retrieve particular species value
    }
  }
  std::vector< std::shared_ptr<Reaction> > reactions;
  // recording results of simulation
  std::vector<std::string> species;
  std::vector< std::vector<int> > concentrations;
  
  // for current simulation run
  unsigned int num_rxns = 0;
  State current_state;
  std::vector<unsigned int> rxn_count;
  std::vector<double> pause_times;
  
  // for random number generation
  std::mt19937_64 gen{std::random_device{}()};
  std::exponential_distribution<double> exp{1.0}; // unit exp rv
  std::uniform_real_distribution<double> uniform{0.0,1.0};

  // for time keeping
  std::vector<double> times;
  TimeKeeper tkeeper;
};

#endif

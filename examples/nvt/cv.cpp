#include <cmath>
#include <iostream>
#include "mmd.hpp"

using namespace std;
using namespace mmd;

int main() {

  static const double dt = 0.002;
  static const double tau_T = 0.05;
  static const long long nsteps = 500000;
  static const double L = 7.4;
  static const double tstar = 100.0 / 121.0; /* Ar temperature scale is 121 K */
 
  const_well_params_LJ_cutoff_potential pot(1.0, 1.0, L, 2.5);
  nose_hoover_velocity_verlet_pbc time_int(tstar, tau_T); 

  cout << "Initializing simulation... to 100K, \n\n";
  simulation sim(molecular_id::Test, "liquid256_init.xyz", &pot, dt, tstar, 
                 time_int, L);

  cout << "Initial temperature... " << temperature(sim) * 121.0 << '\n';

  sim.add_callback(check_momentum(50*dt, 1e-6));

  // store mass positions and velocities
  sim.add_callback(save_vector_with_time_callback<7>("cv_data.csv", 50 * dt, 
                   euktpiv));

  sim.add_callback(print_vector_with_time_callback<7>(cout, 1000*dt, euktpiv));

  cout << "Running simulation for " << nsteps << " steps...\n\n";
  cout << "time E=U+K U K T P I V\n" << "======================\n";

  mprof::tic();
  sim.simulate(nsteps);
  mprof::toc();

  return 0;
}

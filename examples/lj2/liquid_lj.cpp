#include <cmath>
#include <iostream>
#include "mmd.hpp"

using namespace std;
using namespace mmd;

int main() {

  static const double dt = 0.002;
  static const long long nsteps = 100000;
  static const double L = 6.8;
  static const double tstar = 100.0 / 121.0; /* Ar temperature scale is 121 K */
 
  const_well_params_LJ_cutoff_potential pot(1.0, 1.0, L, 2.5);
  periodic_velocity_verlet time_int(L); 

  cout << "Initializing simulation... to 100K, \n\n";
  simulation sim(molecular_id::Test, "liquid256_init.xyz", &pot, dt, tstar, 
                 time_int, L);

  cout << "Initial temperature... " << temperature(sim) * 121.0 << '\n';

  sim.add_mutator(equilibrate_temperature(tstar, 1e-2, 1000));

  //sim.add_callback(check_energy(50*dt, 1e-3));
  sim.add_callback(check_momentum(50*dt, 1e-6));

  // store mass positions and velocities
  sim.add_callback(save_values_with_time_callback("liquid_mom.csv", 50*dt, 
        {
          [](const simulation& sim) -> double { 
            return momentum(sim)(0); // get momentum in the x-direction
          },
          [](const simulation& sim) -> double { 
            return momentum(sim)(1); // get momentum in the y-direction
          },
          [](const simulation& sim) -> double { 
            return momentum(sim)(2); // get momentum in the z-direction
          }
        }));

  // store mass positions and velocities
  sim.add_callback(save_values_with_time_callback("liquid_data.csv", 50 * dt, 
        {
          kinetic_energy, potential_energy, total_energy, temperature, pressure,
          ideal_pressure, virial_pressure
        }));

  sim.add_callback(print_energy_and_momentum_with_time_callback(1000*dt));
  sim.add_callback(save_xyz_callback("liquid_lj.xyz", 100*dt, positions));

  cout << "Running simulation for " << nsteps << " steps...\n\n";
  cout << "time U K E=U+K |p|\n" << "==================\n";

  mprof::tic();
  sim.simulate(nsteps);
  mprof::toc();

  return 0;
}

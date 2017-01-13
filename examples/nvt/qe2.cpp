#include <cmath>
#include <iostream>
#include "mmd.hpp"

using namespace std;
using namespace mmd;

int main() {

  static const double dt = 0.002;
  static const double tau_T = 0.05;
  static const long long nsteps = 100000;
  static const double L = 6.8;
  static const double tstar = 90.0 / 121.0; /* Ar temperature scale is 121 K */
  
  const auto ma = [](molecular_id id) -> double {
    if (id == molecular_id::Test1) return 1.0;
    return 2.0;
  };
 
  const_well_params_LJ_cutoff_potential pot(1.0, 1.0, L, 2.5);
  nose_hoover_velocity_verlet_pbc time_int(tstar, tau_T); 

  cout << "Initializing simulation... to 90K, \n\n";
  simulation sim({molecular_id::Test1, molecular_id::Test2}, ma, 
                 "liquid256_init.xyz", &pot, dt,
                 time_int, L);

  cout << "Testing mass accessor...";
  for (const auto& id : sim.get_molecular_ids()) 
    cout << "id = " << as_integer(id) << ", mass = " << ma(id) << '\n';
  cout << '\n';

  sim.simulate(nsteps / 10);
  sim.reset_clock();

  cout << "Initial temperature... " << temperature(sim) * 121.0 << '\n';

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
  sim.add_callback(save_vector_with_time_callback<7>("liquid_data.csv", 50 * dt, 
                   euktpiv));

 sim.add_callback(print_vector_with_time_callback<7>(cout, 1000*dt, euktpiv));
 sim.add_callback(save_xyz_callback("liquid_lj.xyz", 100*dt, positions));

  cout << "Running simulation for " << nsteps << " steps...\n\n";
  cout << "time E=U+K U K T P I V\n" << "======================\n";

  mprof::tic();
  sim.simulate(nsteps);
  mprof::toc();

  return 0;
}

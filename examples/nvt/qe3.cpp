#include <cmath>
#include <cassert>
#include <iostream>
#include <array>
#include <map>
#include <fstream>
#include "mmd.hpp"

using namespace std;
using namespace mmd;
using namespace arma;

int main() {

  static const size_t ntrials = 10;
  static const size_t N = 256;
  static const double dt = 0.002;
  static const double tau_T = 0.05;
  static const long long nsteps = 75000;
  static const double tstar = 90.0 / 121.0; /* Ar temperature scale is 121 K */
  static const double density_scale = 1686.85;

  // Half ot the molecules have a dimensionless mass of 1.0, the other half
  // have a dimensionless mass of 2.0
  const auto ma = [](molecular_id id) -> double {
    if (id == molecular_id::Test1) return 1.0;
    return 2.0;
  };

  const double density = 1760.0 / density_scale; // zero pressure density
  // average mass of a molecule is 1.5
  const double L = pow(1.5 * static_cast<double>(N) / density, 1.0 / 3.0);
  cout << "L: " << L << '\n';
 
  const_well_params_LJ_cutoff_potential pot(1.0, 1.0, L, 2.5);
  nose_hoover_velocity_verlet_pbc time_int_nvt(tstar, tau_T); 
  simulation sim({molecular_id::Test1, molecular_id::Test2}, ma, 
                 "liquid256_init.xyz", &pot, dt, 
                 time_int_nvt, L);

  for (size_t trial = 0; trial < ntrials; ++trial) {
    // initialize simulation with NVT time integration to equilibrate at 90K
    sim.set_time_int(time_int_nvt);
    sim.simulate(nsteps / 2);

    // now, switch time integrator to NVE ensemble
    sim.set_time_int(velocity_verlet_pbc);
    sim.simulate(nsteps / 20);
    sim.reset_clock();

    sim.add_callback(check_momentum(500*dt, 1e-3));
    sim.add_callback(save_xyz_callback(string("qe3_") + to_string(trial)
          + string(".xyz"), 10 * dt, positions));
    
    sim.add_callback(print_vector_with_time_callback<7>(cout, 1000*dt, euktpiv));
    cout << "Running simulation for " << nsteps << " steps...\n\n";
    cout << "time E=U+K U K T P I V\n" << "======================\n";

    mprof::tic();
    sim.simulate(nsteps);
    mprof::toc();
  }

  return 0;
}

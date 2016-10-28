#include <cmath>
#include <iostream>
#include <array>
#include <map>
#include <fstream>
#include "mmd.hpp"

using namespace std;
using namespace mmd;

int main() {

  static const double dt = 0.002;
  static const long long nsteps = 100000;
  static const double total_time = dt * nsteps;
  static const double L = 6.8;
  static const double tstar = 100.0 / 121.0; /* Ar temperature scale is 121 K */
  static const double density_scale = 1686.85;
  static const array<double, 2> density_bounds { { 950.0 / density_scale, 
                                                   1150.0 / density_scale } };
  static const size_t num_sims = 100;
  static const double drho = (density_bounds[1] - density_bounds[0]) / num_sims;

  map<double, double> pressures;

  #pragma omp parallel for schedule(dynamic)
  for (unsigned k = 0; k < num_sims; ++k) {

    const double density = k * drho + density_bounds[0];
   
    const_well_params_LJ_cutoff_potential pot(1.0, 1.0, L, 2.5);
    
    simulation sim(molecular_id::Test, "liquid256_init.xyz", density, &pot, dt, 
                   tstar, velocity_verlet_pbc, L);
 
    // equilibrate system at 100K
    sim.add_mutator(equilibrate_temperature(tstar, 1e-2, nsteps));
    sim.add_callback(check_momentum(100*dt, 1e-6));
    sim.simulate(nsteps / 10 - 100);
    sim.clear_all_callbacks();
    sim.add_mutator(raise_temperature(tstar, nsteps));
    sim.simulate(100);
    sim.clear_all_callbacks();
    sim.reset_clock();

    sim.add_callback([&](const simulation& sim) {
      if (fmod(sim.get_time(), 1000*dt) < dt) {
        cout << "Density: " << density << ", " 
             << sim.get_time() / total_time * 100 << "%\n";
      }
    });

    double pressure_sum = 0.0;
    size_t nsamples = 0;
    
    sim.add_callback([&](const simulation& sim) {
      if (fmod(sim.get_time(), 5*dt) < dt) {
        pressure_sum += pressure(sim);
        ++nsamples;
      }
    });
    sim.add_callback(check_energy(1000*dt, 1e-3));
    sim.add_callback(check_momentum(1000*dt, 1e-6));

    mprof::tic();
    sim.simulate(nsteps);
    cout << "Density: " << density << ", finished.\n";
    mprof::toc();

    pressures[density] = pressure_sum / nsamples;

  }

  ofstream outfile("pressure_data.csv", ofstream::out);

  for (const auto& pressure_data : pressures)
    outfile << pressure_data.first << ", " << pressure_data.second << '\n';

  outfile.close();

  return 0;
}

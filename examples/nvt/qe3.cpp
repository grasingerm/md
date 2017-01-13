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

  static const size_t N = 256;
  static const double dt = 0.002;
  static const double tau_T = 0.05;
  static const long long nsteps = 100000;
  static const double tstar = 90.0 / 121.0; /* Ar temperature scale is 121 K */
  static const double density_scale = 1686.85;

  // Half ot the molecules have a dimensionless mass of 1.0, the other half
  // have a dimensionless mass of 2.0
  const auto ma = [](molecular_id id) -> double {
    if (id == molecular_id::Test1) return 1.0;
    return 2.0;
  };

  const double density = 1.03; // zero pressure density
  // average mass of a molecule is 1.5
  const double L = pow(1.5 * static_cast<double>(N) / density, 1.0 / 3.0);
 
  const_well_params_LJ_cutoff_potential pot(1.0, 1.0, L, 2.5);

  // initialize simulation with NVT time integration to equilibrate at 90K
  nose_hoover_velocity_verlet_pbc time_int_nvt(tstar, tau_T); 
  simulation sim({molecular_id::Test1, molecular_id::Test2}, ma, 
                 "liquid256_init.xyz", &pot, dt, 
                 time_int_nvt, L);

  sim.simulate(nsteps / 10);
  sim.reset_clock();

  // now, switch time integrator to NVE ensemble
  sim.set_time_int(velocity_verlet_pbc);

  // initialize data structures for calculating diffusion coefficient
  ofstream outfile("diffusion_data.csv", ofstream::out);
  outfile << "t,dx1,dy1,dz1,dj1,dx2,dy2,dz2,dj2\n";

  // copy/store initial positions
  mat r_inits(sim.get_positions());

  // calculate diffusion coefficient
  sim.add_callback(check_momentum(500*dt, 1e-3));
  sim.add_callback([&](const simulation& sim) {

    if (fmod(sim.get_time(), 5*dt) < dt) {

      const auto& ids = sim.get_molecular_ids();
      const auto& positions = sim.get_positions();
      double dx1 = 0, dy1 = 0, dz1 = 0;
      double dx2 = 0, dy2 = 0, dz2 = 0;

      for (size_t i = 0; i < sim.get_N(); ++i) {
        vec r_rel = positions.col(i) - r_inits.col(i);
        
        // correct for periodic boundary conditions
        for (unsigned i = 0; i < 3; ++i) {
          if (r_rel(i) > L / 2.0)
            r_rel(i) -= L;
          else if (r_rel(i) < -L / 2.0)
            r_rel(i) += L;
        }

        if (ids[i] == molecular_id::Test1) {
          dx1 += r_rel(0)*r_rel(0);
          dy1 += r_rel(1)*r_rel(1);
          dz1 += r_rel(2)*r_rel(2);
        }
        else {
          dx2 += r_rel(0)*r_rel(0);
          dy2 += r_rel(1)*r_rel(1);
          dz2 += r_rel(2)*r_rel(2);
        }
      }

      dx1 /= 256;
      dy1 /= 256;
      dz1 /= 256;
      dx2 /= 256;
      dy2 /= 256;
      dz2 /= 256;
     
      // write row of data
      outfile << sim.get_time() 
              << ',' << dx1 << ',' << dy1 << ',' << dz1 << ',' << (dx1 + dy1 + dz1) / 3.0
              << ',' << dx2 << ',' << dy2 << ',' << dz2 << ',' << (dx2 + dy2 + dz2) / 3.0
              << '\n';
    }

  });

  sim.add_callback(print_vector_with_time_callback<7>(cout, 1000*dt, euktpiv));
  cout << "Running simulation for " << nsteps << " steps...\n\n";
  cout << "time E=U+K U K T P I V\n" << "======================\n";

  mprof::tic();
  sim.simulate(nsteps);
  mprof::toc();

  outfile.close();

  return 0;
}

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

  static const size_t ntrials = 3;
  static const size_t N = 256;
  static const double dt = 0.002;
  static const double tau_T = 0.05;
  static const long long nsteps = 75000;
  static const double tstar = 90.0 / 121.0; /* Ar temperature scale is 121 K */
  static const double density_scale = 1686.85;
  static const double pressure_scale = 42;

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
  
  double pressure_sum = 0;
  double temperature_sum = 0;
  size_t nsamples = 0;
  mat travelled(sim.get_positions());
  mat mean_sqs(3, 2);

  for (size_t trial = 0; trial < ntrials; ++trial) {
    sim.add_callback(print_vector_with_time_callback<7>(cout, 1000*dt, euktpiv));
    cout << "\n######################\n\tTrial:\t " << trial+1 
         << "\n######################\n\n";
    cout << "Initializing for " << nsteps / 2 + nsteps / 20 << " steps...\n\n";
    cout << "time E=U+K U K T P I V\n" << "======================\n";
    // initialize simulation with NVT time integration to equilibrate at 90K
    sim.set_time_int(time_int_nvt);
    sim.simulate(nsteps / 2);

    // now, switch time integrator to NVE ensemble
    sim.set_time_int(velocity_verlet_pbc);
    sim.simulate(nsteps / 20);
    sim.reset_clock();

    ofstream os(string("qe_") + to_string(trial) + string(".csv"));
    mat prev_positions(sim.get_positions());
    travelled.zeros();

    sim.add_callback([&](const simulation& sim) {
        if (fmod(sim.get_time(), 10*dt) < dt) {
          pressure_sum += pressure(sim);
          temperature_sum += temperature(sim);
          ++nsamples;

          const auto& positions = sim.get_positions();
          const auto& ids = sim.get_molecular_ids();

          os << sim.get_time();
          for (size_t k = 0; k < sim.get_N(); ++k) {
            vec dr = positions.col(k) - prev_positions.col(k);
            for (unsigned i = 0; i < 3; ++i) {
              if (dr(i) >= L / 2) {
                dr(i) -= L;
              }
              else if(dr(i) <= -L / 2) {
                dr(i) += L;
              }
              travelled(i, k) += dr(i);
              if (ids[k] == molecular_id::Test1) {
                mean_sqs(i, 0) += travelled(i, k) * travelled(i, k);
              }
              else {
                mean_sqs(i, 1) += travelled(i, k) * travelled(i, k);
              }
            }
          }
          for (unsigned j = 0; j < 2; ++j) {
            double sum = 0.0;
            for (unsigned i = 0; i < 3; ++i) {
              os << ',' <<  mean_sqs(i, j);
              sum += mean_sqs(i, j);
            }
            os << ',' << sum / 3.0;
          }
          os << '\n';
          prev_positions = sim.get_positions();
        }
    });
    sim.add_callback(check_momentum(500*dt, 1e-3));
    sim.add_callback(save_xyz_callback(string("qe3_") + to_string(trial)
          + string(".xyz"), 10 * dt, positions));
    
    cout << "Running simulation for " << nsteps << " steps...\n\n";
    cout << "time E=U+K U K T P I V\n" << "======================\n";

    mprof::tic();
    sim.simulate(nsteps);
    mprof::toc();
    
    cout << "\n<P> = " << pressure_sum / nsamples * pressure_scale << '\n';
    cout << "<T> = " << temperature_sum / nsamples * 121.0 << '\n';

    sim.clear_all_callbacks();
  }
   
  cout << "\n\nFinal ensemble averages:\n======================\n";
  cout << "\n<P> = " << pressure_sum / nsamples * pressure_scale << '\n';
  cout << "<T> = " << temperature_sum / nsamples * 121.0 << '\n';

  return 0;
}

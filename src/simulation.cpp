#include <iostream>
#include <fstream>
#include <random>
#include <string>
#include <stdexcept>
#include <cassert>
#include "simulation.hpp"
#include "distance.hpp"

namespace mmd {

using namespace std;
using namespace arma;

// default time integration scheme
const time_integrator simulation::default_time_int = velocity_verlet;

simulation::simulation(const molecular_id id, const char* fname, 
                       abstract_potential* pot, const double dt,
                       const double edge_length)
  : ma([](molecular_id) -> double { return 1.0; }), potentials(size_t{1}, pot), 
    time_int(simulation::default_time_int), dt(dt), t(0.0), volume(edge_length *
    edge_length * edge_length), kB(simulation::default_kB) {

  string line;
  ifstream infile(fname, ifstream::in);
  
  if (!infile.is_open()) 
    throw invalid_argument(string("Could not open data file: ") + string(fname));

  size_t n;
  infile >> n;

  // allocate memory for molecular degrees of freedom
  positions.resize(3, n); 
  velocities.zeros(3, n); 
  forces.zeros(3, n); 

  getline(infile, line); // read-in newline
  getline(infile, line); // comment line

  size_t i = 0;
  double x, y, z;

  // read positions from file and initialize random velocities
  while ((infile >> x >> y >> z) && i < n) {
    molecular_ids.push_back(id);

    positions(0, i) = x;
    positions(1, i) = y;
    positions(2, i) = z;

    ++i;
  }

  infile.close();

}

simulation::simulation(const molecular_id id, const char* fname, 
                       abstract_potential* pot, const double dt, 
                       const double tstar, const double edge_length) 
  : simulation(id, fname, pot, dt, edge_length) {

  const double vscale = sqrt(2.0 * tstar * 3.0 * (get_N() - 1) / 2.0);

  default_random_engine generator;
  normal_distribution<double> distx(0.0, vscale);
  normal_distribution<double> disty(0.0, vscale);
  normal_distribution<double> distz(0.0, vscale);

  for (size_t i = 0; i < velocities.n_cols; ++i) {
    velocities(0, i) = distx(generator);
    velocities(1, i) = disty(generator);
    velocities(2, i) = distz(generator);
  }

  const auto n = molecular_ids.size();

  // ensure initial momentum is zero
  // Note:: m is constant for this constructor so we only need to consider
  //        velocities
  for (unsigned i = 0; i < 3; ++i) {
    double vsum = 0.0;
    for (size_t j = 0; j < n; ++j) vsum += velocities(i, j);
    for (size_t j = 0; j < n; ++j) velocities(i, j) -= vsum / n;
  }

  // check momentum
#ifndef NDEBUG
  for (unsigned i = 0; i < 3; ++i) {
    double sum = 0.0;
    for (size_t j = 0; j < n; ++j) sum += velocities(i, j);
    assert(abs(sum) < 1.0e-9);
  }
#endif

  // scale temperature
  const double sqrt_alpha = sqrt(tstar / temperature(*this));
  for (size_t j = 0; j < n; ++j)
    velocities.col(j) *= sqrt_alpha;

  assert(abs(temperature(*this) - tstar) < 1e-9);
}

simulation::simulation(const molecular_id id, const char* fname_pos,
                       const char* fname_vel, abstract_potential* pot, 
                       const double dt, const double edge_length) 
  : simulation(id, fname_pos, pot, dt, edge_length) {

  string line;
  ifstream infile(fname_vel, ifstream::in);
  
  if (!infile.is_open()) 
    throw invalid_argument(string("Could not open data file: ") 
                           + string(fname_vel));

  size_t n;
  infile >> n;

  getline(infile, line); // read in newline
  getline(infile, line); // comment line

  size_t i = 0;
  double u, v, w;

  // read positions from file and initialize random velocities
  while ((infile >> u >> v >> w) && i < n) {

    velocities(0, i) = u;
    velocities(1, i) = v;
    velocities(2, i) = w;

    ++i;

  }

  infile.close();
}

void simulation::simulate(const unsigned long long nsteps) {

  try {
    for (unsigned k = 0; k < nsteps; ++k) {

      // Integrate in time
      time_int(potentials, molecular_ids, positions, velocities, forces, dt, ma);

      // Update clock 
      t += dt;

      // Call all of the callback functions
      for (auto& cb : callbacks) cb(*this);
      
    }
  }
  catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << '\n';
    throw e;
  }
}

double kinetic_energy(const simulation& sim) {
  
  const auto& molecular_ids = sim.get_molecular_ids();
  const auto& ma = sim.get_mass_accessor();
  const auto& velocities = sim.get_velocities();

  double sum = 0.0;
  for (size_t i = 0; i < molecular_ids.size(); ++i)
    sum += ma(molecular_ids[i]) * dot(velocities.col(i), velocities.col(i));

  return 0.5 * sum;

}

arma::vec momentum(const simulation& sim) {
  
  const auto& molecular_ids = sim.get_molecular_ids();
  const auto& ma = sim.get_mass_accessor();
  const auto& velocities = sim.get_velocities();

  arma::vec vsum(3); vsum.zeros();
  for (size_t i = 0; i < molecular_ids.size(); ++i) 
    vsum += (ma(molecular_ids[i]) * velocities.col(i));

  return vsum;

}

inline double _temperature_kernel(const double kinetic_energy, const double N,
                                  const double kB) {
  return 2.0 * kinetic_energy / (3.0 * (N - 1) * kB);
}

double temperature(const simulation& sim) {
  return _temperature_kernel(kinetic_energy(sim), sim.get_N(), sim.get_kB());
}

inline double _ideal_pressure_kernel(const double temperature, const double N,
                                     const double volume, const double kB){
  return N * temperature * kB / volume;
}

double ideal_pressure(const simulation& sim) {
  return _ideal_pressure_kernel(temperature(sim), sim.get_N(), sim.get_volume(),
                                sim.get_kB());
}

/* TODO: can't I just use the forces calculated from integration? sum Fi dot ri?
 * TODO: the periodic BC concept has become a mess..., feels like it would be
 *       much better to encapsulate the concept of a metric in the simulation
 *       object itself (that then forwards that information along), either
 *       using inheritance or composition. Gut says composition
 */
double virial_pressure(const simulation& sim) {
  const auto n = sim.get_N();
  const auto& positions = sim.get_positions();
  double vp = 0.0;

  for (const auto potential : sim.get_potentials())
    #pragma omp parallel for reduction(+:vp) schedule(dynamic)
    for (auto i = size_t{0}; i < n-1; ++i)
      for (auto j = i+1; j < n; ++j)
        // rij dot Fij
        vp += dot(rij_pbc(positions, i, j, sim.get_edge_length()), 
                  potential->force_ij(sim.get_molecular_ids(), 
                                      sim.get_positions(), i, j));

  return vp / (3.0 * sim.get_volume());
}

array<double, 3> ktp(const simulation& sim) {

  array<double, 3> result;

  result[0] = kinetic_energy(sim);
  result[1] = _temperature_kernel(result[0], sim.get_N(), sim.get_kB());
  result[2] = _ideal_pressure_kernel(result[1], sim.get_N(), sim.get_volume(),
                                     sim.get_kB()) + virial_pressure(sim);

  return result;

}

array<double, 5> euktp(const simulation& sim) {

  array<double, 5> result;

  result[2] = kinetic_energy(sim);
  result[3] = _temperature_kernel(result[2], sim.get_N(), sim.get_kB());
  result[4] = _ideal_pressure_kernel(result[3], sim.get_N(), sim.get_volume(),
                                     sim.get_kB()) + virial_pressure(sim);
  result[1] = potential_energy(sim);
  result[0] = result[1] + result[2];

  return result;

}

} // namespace mmd

#include <fstream>
#include <random>
#include <string>
#include <stdexcept>
#include <cassert>
#include "simulation.hpp"

namespace mmd {

using namespace std;
using namespace arma;

// default time integration scheme
const time_integrator simulation::default_time_int = velocity_verlet;

simulation::simulation(const molecular_id id, const char* fname, 
                       abstract_potential* pot, const double dt)
  : ma([](molecular_id) -> double { return 1.0; }), potentials(size_t{1}, pot), 
    time_int(simulation::default_time_int), dt(dt), t(0.0) {

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

    cout << "reading in values: " << x << ", " << y << ", " << z << '\n';

    positions(0, i) = x;
    positions(1, i) = y;
    positions(2, i) = z;

    ++i;
  }

  infile.close();

}

simulation::simulation(const molecular_id id, const char* fname, 
                       abstract_potential* pot, const double dt, 
                       const double vscale) : simulation(id, fname, pot, dt) {

  default_random_engine generator;
  normal_distribution<double> dist(0.0, vscale);

  for (size_t i = 0; i < velocities.n_cols; ++i) {
    velocities(0, i) = dist(generator);
    velocities(1, i) = dist(generator);
    velocities(2, i) = dist(generator);
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
}

simulation::simulation(const molecular_id id, const char* fname_pos,
                       const char* fname_vel, abstract_potential* pot, 
                       const double dt) 
  : simulation(id, fname_pos, pot, dt) {

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
  for (unsigned k = 0; k < nsteps; ++k) {

    // Integrate in time
    time_int(potentials, molecular_ids, positions, velocities, forces, dt, ma);

    // Update clock 
    t += dt;

    // Call all of the callback functions
    for (auto& cb : callbacks) cb(*this);
    
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

} // namespace mmd

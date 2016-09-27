#include <fstream>
#include <random>
#include <string>
#include <stdexcept>
#include <cassert>
#include "simulation.hpp"

namespace mmd {

using namespace std;
using namespace arma;

simulation::simulation(const molecular_id id, const char* fname, 
                       const abstract_potential* pot, const double dt,
                       const double vscale) 
  : ma([](molecular_id) { return 1.0; }), potentials(pot, 1), 
    time_int(simulation::default_time_int), dt(dt), t(0.0) {

  string line;
  ifstream infile(fname, ifstream::in);
  
  if (!infile.is_open) 
    throw invalid_argument(string("Could not open data file: ") + string(fname));

  int n;
  infile >> n;

  // allocate memory for molecular degrees of freedom
  positions.resize(3, n); 
  velocities.resize(3, n); 
  forces.zeros(3, n); 

  getline(infile, line); // comment line

  size_t i = 0;
  double x, y, z;
  default_random_engine generator;
  normal_distribution<double> dist(0.0, vscale);

  // read positions from file and initialize random velocities
  while ((infile >> x >> y >> z) && i < n) {
    molecular_ids.push_back(id);

    positions(i, 0) = x;
    positions(i, 1) = y;
    positions(i, 2) = z;

    velocities(i, 0) = dist(generator);
    velocities(i, 1) = dist(generator);
    velocities(i, 2) = dist(generator);
  }

  infile.close();

  // ensure initial momentum is zero
  // Note:: m is constant for this constructor so we only need to consider
  //        velocities
  for (unsigned j = 0; j < 3; ++j) {
    double vsum = 0.0;
    for (size_t i = 0; i < n; ++i) vsum += velocities(i, j);
    for (size_t i = 0; i < n; ++i) velocities(i, j) -= vsum / n;
  }

  assert(abs(sum(velocities)) < 1.0e-9);
}

void simulation::simulate(const unsigned nsteps) {
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

  auto sum = arma::zeros<vec>(3);
  for (size_t i = 0; i < molecular_ids.size(); ++i)
    sum += ma(molecular_ids[i]) * velocities.col(i);

  return sum;

}

} // namespace mmd

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
                       const double vscale = 1.0) 
  : potentials(pot, 1), time_int(simulation::default_time_int), dt(dt) {

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
  const n = molecules.size();
  double t = 0.0;
  vector<array<double, D>> forces(n, {0.0, 0.0});

  for (unsigned k = 0; k < nsteps; ++k) {

    for (auto& force : forces) {
      force[0] = 0.0;
      force[1] = 0.0;
    }

    // Add single particle potentials to instantaneous forces
    for (const auto ∇φ1 : ∇φ1s) {
      for (unsigned i = 0; i < n; ++i) {
        forces[i] += ∇φ1(molecules[i].q);
      }
    }

    // Add pairwise particle potentials to instantaneous forces
    for (const auto ∇φ2 : ∇φ2s) {
      for (unsigned i = 0; i < n-1; ++i) {
        for (unsigned j = i; j < n; ++j) {
          const unordered_pair<molecular_id> key { molecules[i].id, 
                                                   molecules[j].id };
          const auto& params = interactions[key];
          const auto ∇φij = ∇φ2(params[0], params[1], 
                                molecules[i].q - molecules[j].q);
          forces[i] += ∇φij;
          forces[j] -= ∇φij;
        }
      }
    }

    // Calculate instantaneous acceleration and integrate
    for (unsigned i = 0; i < n; ++i) {
      const mass = masses[molecules[i].id];
      const array<double, D> a { forces[i][0] / mass, forces[i][1] / mass };
      time_int(molecules[i], a, Δt); 
    }

    t += Δt;

    for (const auto callback : callbacks) callback(molecules, t);

  }
}

} // namespace mmd

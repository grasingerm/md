#include <cmath>
#include <iostream>
#include "callbacks.hh"

namespace mmd {

using namespace std;

save_xyz_callback::operator()(const simulation& sim) {
  
  outfile << sim.molecular_ids().size() << '\n';
  outfile << "time: " << sim.get_time() << '\n';
  
  const auto& data = da(sim);
  for (size_t j = 0; j < data.n_cols; ++j) {
    outfile << sim.molecular_ids()[j];
    for (unsigned i = 0; i < 3; ++i)
      outfile << ' ' << data(i, j);
    outfile << '\n';
  }

}

save_values_with_time_callback::operator()(const simulation& sim) {
  outfile << sim.get_time();
  for (const auto& va : vas) outfile << delim << va(sim);
  outfile << '\n';
}

callback check_energy(const double dt, const double eps) {
 
  double init_energy;
  bool init_call = true;
  
  return [=](const simulation& sim) mutable {
    if (fmod(sim.get_time(), dt) < sim.get_dt()) {
      if (init_call) {
        init_energy = total_energy(sim);
        init_call = false;
      }
      else {
        if (abs(init_energy - total_energy(sim)) / abs(init_energy) > eps) {
          cerr << "Energy is diverging outside the tolerance limit\n";
          cerr << "Initial energy: " << init_energy << ", current energy: "
               << total_energy(sim) << ", tolerance: " << eps << '\n';
          throw runtime_error("Energy is not being conserved");
        }
      }
    }
  };

}

callback check_momentum(const double dt, const double eps) {
 
  arma::vec init_momentum;
  bool init_call = true;
  
  return [=](const simulation& sim) mutable {
    if (fmod(sim.get_time(), dt) < sim.get_dt()) {
      if (init_call) {
        init_momentum = momentum(sim);
        init_call = false;
      }
      else {
        if (norm(init_momentum - momentum(sim)) / norm(init_momentum) > eps) {
          cerr << "System momentum is diverging outside the tolerance limit\n";
          cerr << "Initial momentum: " << init_momentum << ", current momentum: "
               << momentum(sim) << ", tolerance: " << eps << '\n';
          throw runtime_error("Momentum is not being conserved");
        }
      }
    }
  };

}

}

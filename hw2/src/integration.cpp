#include "integration.hpp"

namespace mmd {

using namespace arma;

void euler(std::vector<abstract_potential*>& potentials, 
           std::vector<molecular_id>& molecular_ids, 
           mat& positions, mat& velocities, mat& forces, const double dt) {
 
  _check_arg_sizes(molecular_ids, mat& positions, mat& velocities, mat& forces);

  forces.zeros();
  for (const auto& potential : potentials)
    potential->increment_forces(molecular_ids, positions, forces);

  for (auto i = size_t{0}; i < molecular_ids.size(); ++i) {
    velocities.col(i) += forces.col(i) / m * dt;
    positions.col(i) += velocities.col(i) * dt; 
  }

}

void velocity_verlet(std::vector<abstract_potential*>& potentials, 
                     std::vector<molecular_id>& molecular_ids, 
                     mat& positions, mat& velocities, mat& forces, 
                     const double dt) {
  
  throw "Not yet implemented"; 
  _check_arg_sizes(molecular_ids, mat& positions, mat& velocities, mat& forces);

  forces.zeros();
  for (const auto& potential : potentials)
    potential->increment_forces(molecular_ids, positions, forces);

  for (auto i = size_t{0}; i < molecular_ids.size(); ++i) {
    velocities.col(i) += forces.col(i) / m * dt;
    positions.col(i) += velocities.col(i) * dt; 
  }

}

} // namespace mmd

#include "integration.hpp"
#include "boundary.hpp"
#include "debug.hpp"

namespace mmd {

using namespace arma;

void euler(std::vector<abstract_potential *> &potentials,
           std::vector<molecular_id> &molecular_ids, mat &positions,
           mat &velocities, mat &forces, const double dt,
           const mass_accessor &ma) {

  _check_arg_sizes(molecular_ids, positions, velocities, forces);

  forces.zeros();
  for (const auto &potential : potentials)
    potential->increment_forces(molecular_ids, positions, forces);

  for (auto i = size_t{0}; i < molecular_ids.size(); ++i) {
    velocities.col(i) += forces.col(i) / ma(molecular_ids[i]) * dt;
    positions.col(i) += velocities.col(i) * dt;
  }
}

void euler_pbc(std::vector<abstract_potential *> &potentials,
               std::vector<molecular_id> &molecular_ids, mat &positions,
               mat &velocities, mat &forces, const double dt,
               const mass_accessor &ma, const double edge_length) {

  _check_arg_sizes(molecular_ids, positions, velocities, forces);

  forces.zeros();
  for (const auto &potential : potentials)
    potential->increment_forces(molecular_ids, positions, forces);

  for (auto i = size_t{0}; i < molecular_ids.size(); ++i) {
    velocities.col(i) += forces.col(i) / ma(molecular_ids[i]) * dt;
    positions.col(i) += velocities.col(i) * dt;
    apply_pbc(positions, edge_length, i);
  }
}

void velocity_verlet(std::vector<abstract_potential *> &potentials,
                     std::vector<molecular_id> &molecular_ids, mat &positions,
                     mat &velocities, mat &forces, const double dt,
                     const mass_accessor &ma) {

  _check_arg_sizes(molecular_ids, positions, velocities, forces);

  // calculate F(t)
  forces.zeros();
  for (const auto &potential : potentials)
    potential->increment_forces(molecular_ids, positions, forces);

  // calculate v(t+dt/2) and r(t+dt)
  for (auto i = size_t{0}; i < molecular_ids.size(); ++i) {
    const double m = ma(molecular_ids[i]);
    velocities.col(i) += forces.col(i) / (2 * m) * dt;
    positions.col(i) += velocities.col(i) * dt;
  }

  // recalculate forces, we want F(t+dt)
  forces.zeros();
  for (const auto &potential : potentials)
    potential->increment_forces(molecular_ids, positions, forces);

  // calculate v(t+dt)
  for (auto i = size_t{0}; i < molecular_ids.size(); ++i) {
    const double m = ma(molecular_ids[i]);
    velocities.col(i) += forces.col(i) / (2 * m) * dt;
  }
}

void velocity_verlet_pbc(std::vector<abstract_potential *> &potentials,
                         std::vector<molecular_id> &molecular_ids,
                         mat &positions, mat &velocities, mat &forces,
                         const double dt, const mass_accessor &ma,
                         const double edge_length) {

  _check_arg_sizes(molecular_ids, positions, velocities, forces);

  // calculate F(t)
  forces.zeros();
  for (const auto &potential : potentials)
    potential->increment_forces(molecular_ids, positions, forces);

  // calculate v(t+dt/2) and r(t+dt)
  for (auto i = size_t{0}; i < molecular_ids.size(); ++i) {
    const double m = ma(molecular_ids[i]);
    velocities.col(i) += forces.col(i) / (2 * m) * dt;
    positions.col(i) += velocities.col(i) * dt;
    apply_pbc(positions, edge_length, i);
  }

  // recalculate forces, we want F(t+dt)
  forces.zeros();
  for (const auto &potential : potentials)
    potential->increment_forces(molecular_ids, positions, forces);

  // calculate v(t+dt)
  for (auto i = size_t{0}; i < molecular_ids.size(); ++i) {
    const double m = ma(molecular_ids[i]);
    velocities.col(i) += forces.col(i) / (2 * m) * dt;
  }
}

void quenched_velocity_verlet::
operator()(std::vector<abstract_potential *> &potentials,
           std::vector<molecular_id> &molecular_ids, mat &positions,
           mat &velocities, mat &forces, const double dt,
           const mass_accessor &ma) const {

  _check_arg_sizes(molecular_ids, positions, velocities, forces);

  // calculate F(t)
  forces.zeros();
  for (const auto &potential : potentials)
    potential->increment_forces(molecular_ids, positions, forces);

  // calculate v(t+dt/2) and r(t+dt)
  for (auto i = size_t{0}; i < molecular_ids.size(); ++i) {
    const double m = ma(molecular_ids[i]);
    velocities.col(i) +=
        (forces.col(i) / m - eta * velocities.col(i)) * dt / 2.0;
    positions.col(i) += velocities.col(i) * dt;
  }

  // recalculate forces, we want F(t+dt)
  forces.zeros();
  for (const auto &potential : potentials)
    potential->increment_forces(molecular_ids, positions, forces);

  // calculate v(t+dt)
  for (auto i = size_t{0}; i < molecular_ids.size(); ++i) {
    const double m = ma(molecular_ids[i]);
    velocities.col(i) = (velocities.col(i) + forces.col(i) * dt / (2.0 * m)) /
                        (1.0 + eta * dt / 2.0);
  }
}

} // namespace mmd

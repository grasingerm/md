#include "integration.hpp"
#include "boundary.hpp"
#include "debug.hpp"
#include "simulation.hpp"

namespace mmd {

using namespace arma;

void euler(simulation &sim) {
  auto &forces = sim.get_forces();
  const auto &potentials = sim.get_potentials();
  const auto &molecular_ids = sim.get_molecular_ids();
  const auto &ma = sim.get_mass_accessor();
  auto &positions = sim.get_positions();
  auto &velocities = sim.get_velocities();
  const auto dt = sim.get_dt();

  forces.zeros();
  for (const auto &potential : potentials)
    potential->increment_forces(molecular_ids, positions, forces);

  for (auto i = size_t{0}; i < molecular_ids.size(); ++i) {
    velocities.col(i) += forces.col(i) / ma(molecular_ids[i]) * dt;
    positions.col(i) += velocities.col(i) * dt;
  }
}

void euler_pbc(simulation &sim) {
  auto &forces = sim.get_forces();
  const auto &potentials = sim.get_potentials();
  const auto &molecular_ids = sim.get_molecular_ids();
  const auto &ma = sim.get_mass_accessor();
  auto &positions = sim.get_positions();
  auto &velocities = sim.get_velocities();
  const auto dt = sim.get_dt();

  forces.zeros();
  for (const auto &potential : potentials)
    potential->increment_forces(molecular_ids, positions, forces);

  for (auto i = size_t{0}; i < molecular_ids.size(); ++i) {
    velocities.col(i) += forces.col(i) / ma(molecular_ids[i]) * dt;
    positions.col(i) += velocities.col(i) * dt;
    apply_pbc(positions, sim.get_edge_length(), i);
  }
}

void velocity_verlet(simulation &sim) {
  auto &forces = sim.get_forces();
  const auto &potentials = sim.get_potentials();
  const auto &molecular_ids = sim.get_molecular_ids();
  const auto &ma = sim.get_mass_accessor();
  auto &positions = sim.get_positions();
  auto &velocities = sim.get_velocities();
  const auto dt = sim.get_dt();

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

void velocity_verlet_pbc(simulation &sim) {
  auto &forces = sim.get_forces();
  const auto &potentials = sim.get_potentials();
  const auto &molecular_ids = sim.get_molecular_ids();
  const auto &ma = sim.get_mass_accessor();
  auto &positions = sim.get_positions();
  auto &velocities = sim.get_velocities();
  const auto dt = sim.get_dt();

  // calculate F(t)
  forces.zeros();
  for (const auto &potential : potentials)
    potential->increment_forces(molecular_ids, positions, forces);

  // calculate v(t+dt/2) and r(t+dt)
  for (auto i = size_t{0}; i < molecular_ids.size(); ++i) {
    const double m = ma(molecular_ids[i]);
    velocities.col(i) += forces.col(i) / (2 * m) * dt;
    positions.col(i) += velocities.col(i) * dt;
    apply_pbc(positions, sim.get_edge_length(), i);
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

void quenched_velocity_verlet::operator()(simulation &sim) const {
  auto &forces = sim.get_forces();
  const auto &potentials = sim.get_potentials();
  const auto &molecular_ids = sim.get_molecular_ids();
  const auto &ma = sim.get_mass_accessor();
  auto &positions = sim.get_positions();
  auto &velocities = sim.get_velocities();
  const auto dt = sim.get_dt();

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

void nose_hoover_velocity_verlet_pbc::operator()(simulation &sim) {
  auto &forces = sim.get_forces();
  const auto &potentials = sim.get_potentials();
  const auto &molecular_ids = sim.get_molecular_ids();
  const auto &ma = sim.get_mass_accessor();
  auto &positions = sim.get_positions();
  auto &velocities = sim.get_velocities();
  const auto dt = sim.get_dt();
  const auto T = temperature(sim);

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
    apply_pbc(positions, sim.get_edge_length(), i);
  }

  // update eta
  eta += dt / tau_Tsq * (T / T_set - 1);

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

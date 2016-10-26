#ifndef __INTEGRATION_HPP__
#define __INTEGRATION_HPP__

#include "mmd_types.hpp"
#include "potentials.hpp"
#include <armadillo>
#include <functional>
#include <vector>

namespace mmd {

/*! \brief Time integration using Euler integration
 *
 * Update the position and velocity of a molecule due to the instantaneous
 * forces acting on it for a given time step. The molecular positions
 * (velocities, and forces) are stored in a matrix such that the ith column
 * vector is the position (velocity, and force) vector of the ith molecule.
 * This function is a mutator that updates the positions, velocities, and
 * forces. This calculation is first order accurate in dt.
 * Requires: dt > 0
 *
 * \param   sim            Simulation object
 */
void euler(simulation &sim);

/*! \brief Time integration using Euler integration
 *
 * Update the position and velocity of a molecule due to the instantaneous
 * forces acting on it for a given time step. The molecular positions
 * (velocities, and forces) are stored in a matrix such that the ith column
 * vector is the position (velocity, and force) vector of the ith molecule.
 * This function is a mutator that updates the positions, velocities, and
 * forces. It applies periodic boundary conditions during integration.
 * This calculation is first order accurate in dt.
 * Requires: dt > 0
 *
 * \param   sim   Simulation object
 */
void euler_pbc(simulation &sim);

/*! \brief Time integration using Velocity Verlet algorithm
 *
 * Update the position and velocity of a molecule due to the instantaneous
 * forces acting on it for a given time step. The molecular positions
 * (velocities, and forces) are stored in a matrix such that the ith column
 * vector is the position (velocity, and force) vector of the ith molecule.
 * This function is a mutator that updates the positions, velocities, and
 * forces. This calculation is fourth order accurate in dt.
 * Requires: dt > 0
 *
 * \param   sim            Simulation object
 */
void velocity_verlet(simulation &sim);

/*! \brief Time integration using Velocity Verlet algorithm
 *
 * Update the position and velocity of a molecule due to the instantaneous
 * forces acting on it for a given time step. The molecular positions
 * (velocities, and forces) are stored in a matrix such that the ith column
 * vector is the position (velocity, and force) vector of the ith molecule.
 * This function is a mutator that updates the positions, velocities, and
 * forces. It applies periodic boundary conditions during integration.
 * This calculation is fourth order accurate in dt.
 * Requires: dt > 0
 *
 * \param   sim    Simulation object
 */
void velocity_verlet_pbc(simulation &sim);

/*! \brief Time integration using a "quenching", Velocity Verlet algorithm
 *
 * Update the position and velocity of a molecule due to the instantaneous
 * forces acting on it for a given time step. NOTE: this integration scheme
 * slowly removes kinetic energy from the system.
 * Requires: dt > 0
 */
class quenched_velocity_verlet {
public:
  /*! \brief Constructor for quenched velocity verlet
   *
   * \param   eta   Quenching parameter
   * \return        Quenched time integration scheme
   */
  quenched_velocity_verlet(const double eta) : eta(eta) {}

  /*! \brief Time integration using a "quenching", Velocity Verlet algorithm
   *
   * Update the position and velocity of a molecule due to the instantaneous
   * forces acting on it for a given time step. The molecular positions
   * (velocities, and forces) are stored in a matrix such that the ith column
   * vector is the position (velocity, and force) vector of the ith molecule.
   * This function is a mutator that updates the positions, velocities, and
   * forces. This calculation is fourth order accurate in dt. NOTE: this
   * integration scheme slowly removes kinetic energy from the system.
   * Requires: dt > 0
   *
   * \param   sim            Simulation object
   */
  void operator()(simulation &sim) const;

private:
  double eta;
};

/*! \brief Time integration for NVT ensemble using a Verlet with periodic BCs
 *
 * Update the position and velocity of a molecule due to the instantaneous
 * forces acting on it for a given time step. Uses the Nose-Hoover thermostat
 * equations of motion and velocity Verlet integration.
 * Requires: dt > 0
 */
class nose_hoover_velocity_verlet_pbc {
public:
  /*! \brief Constructor for Nose-Hoover velocity verlet with periodic BCs
   *
   * \param   T_set         Set temperature
   * \param   tau_T         Thermostat time constant
   * \return                Nose-Hoover, periodic velocity verlet
   */
  nose_hoover_velocity_verlet_pbc(const double T_set, const double tau_T)
      : eta(0.0), T_set(T_set), tau_Tsq(tau_T*tau_T) {}

  /*! \brief Time integration using Nose-Hoover Verlet with periodic BCs
   *
   * Update the position and velocity of a molecule due to the instantaneous
   * forces acting on it for a given time step. The molecular positions
   * (velocities, and forces) are stored in a matrix such that the ith column
   * vector is the position (velocity, and force) vector of the ith molecule.
   * This function is a mutator that updates the positions, velocities, and
   * forces. This calculation is fourth order accurate in dt.
   * Requires: dt > 0
   *
   * \param   sim           Simulation object
   */
  void operator()(simulation &sim);

private:
  mutable double eta;
  double T_set;
  double tau_Tsq;
};

} // namespace mmd

#endif

#ifndef __INTEGRATION_HPP__
#define __INTEGRATION_HPP__

#include <vector>
#include <functional>
#include <armadillo>
#include "mmd_types.hpp"
#include "potentials.hpp"

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
 * \param   potentials    Potential energy functions (e.g. spring, LJ, etc.)
 * \param   molecular_ids Collection of molecular identities
 * \param   positions     Matrix of molecular positions
 * \param   velocities    Matrix of molecular velocities
 * \param   forces        Matrix of molecular forces
 * \param   dt            Time step size
 * \param   ma            Function for looking up molecular mass
 */
void euler(std::vector<abstract_potential*>& potentials, 
           std::vector<molecular_id>& molecular_ids, 
           arma::mat& positions, arma::mat& velocities, 
           arma::mat& forces, const double dt, const mass_accessor& ma);

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
 * \param   potentials    Potential energy functions (e.g. spring, LJ, etc.)
 * \param   molecular_ids Collection of molecular identities
 * \param   positions     Matrix of molecular positions
 * \param   velocities    Matrix of molecular velocities
 * \param   forces        Matrix of molecular forces
 * \param   dt            Time step size
 * \param   ma            Function for looking up molecular mass
 * \param   edge_length   Edge length of the control volume cube
 */
void euler_pbc(std::vector<abstract_potential*>& potentials, 
               std::vector<molecular_id>& molecular_ids, 
               arma::mat& positions, arma::mat& velocities, 
               arma::mat& forces, const double dt, const mass_accessor& ma,
               const double edge_length);

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
 * \param   potentials    Potential energy functions (e.g. spring, LJ, etc.)
 * \param   molecular_ids Collection of molecular identities
 * \param   positions     Matrix of molecular positions
 * \param   velocities    Matrix of molecular velocities
 * \param   forces        Matrix of molecular forces
 * \param   dt            Time step size
 * \param   ma            Function for looking up molecular mass
 */
void velocity_verlet(std::vector<abstract_potential*>& potentials, 
                     std::vector<molecular_id>& molecular_ids, 
                     arma::mat& positions, arma::mat& velocities, 
                     arma::mat& forces, const double dt,
                     const mass_accessor& ma);

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
 * \param   potentials    Potential energy functions (e.g. spring, LJ, etc.)
 * \param   molecular_ids Collection of molecular identities
 * \param   positions     Matrix of molecular positions
 * \param   velocities    Matrix of molecular velocities
 * \param   forces        Matrix of molecular forces
 * \param   dt            Time step size
 * \param   ma            Function for looking up molecular mass
 * \param   edge_length   Edge length of the control volume cube
 */
void velocity_verlet_pbc(std::vector<abstract_potential*>& potentials, 
                         std::vector<molecular_id>& molecular_ids, 
                         arma::mat& positions, arma::mat& velocities, 
                         arma::mat& forces, const double dt,
                         const mass_accessor& ma, const double edge_length);

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
   * \param   potentials    Potential energy functions (e.g. spring, LJ, etc.)
   * \param   molecular_ids Collection of molecular identities
   * \param   positions     Matrix of molecular positions
   * \param   velocities    Matrix of molecular velocities
   * \param   forces        Matrix of molecular forces
   * \param   dt            Time step size
   * \param   ma            Function for looking up molecular mass
   */
  void operator()(std::vector<abstract_potential*>& potentials, 
                  std::vector<molecular_id>& molecular_ids, 
                  arma::mat& positions, arma::mat& velocities, 
                  arma::mat& forces, const double dt,
                  const mass_accessor& ma) const;

private:
  double eta;
};

/*! \brief Time integration using a Velocity Verlet algorithm with periodic BCs
 *
 * Update the position and velocity of a molecule due to the instantaneous
 * forces acting on it for a given time step.
 * Requires: dt > 0
 */
class periodic_velocity_verlet {
public:
  /*! \brief Constructor for velocity verlet with periodic BCs
   *
   * \param   edge_length   Edge length of control volume
   * \return                Periodic velocity verlet
   */
  periodic_velocity_verlet(const double edge_length) : edge_length(edge_length) {}

  /*! \brief Time integration using Velocity Verlet algorithm with periodic BCs
   *
   * Update the position and velocity of a molecule due to the instantaneous
   * forces acting on it for a given time step. The molecular positions 
   * (velocities, and forces) are stored in a matrix such that the ith column 
   * vector is the position (velocity, and force) vector of the ith molecule.
   * This function is a mutator that updates the positions, velocities, and 
   * forces. This calculation is fourth order accurate in dt.
   * Requires: dt > 0
   *
   * \param   potentials    Potential energy functions (e.g. spring, LJ, etc.)
   * \param   molecular_ids Collection of molecular identities
   * \param   positions     Matrix of molecular positions
   * \param   velocities    Matrix of molecular velocities
   * \param   forces        Matrix of molecular forces
   * \param   dt            Time step size
   * \param   ma            Function for looking up molecular mass
   */
  inline void operator()(std::vector<abstract_potential*>& potentials, 
                  std::vector<molecular_id>& molecular_ids, 
                  arma::mat& positions, arma::mat& velocities, 
                  arma::mat& forces, const double dt,
                  const mass_accessor& ma) const {
    velocity_verlet_pbc(potentials, molecular_ids, positions, velocities, forces,
                        dt, ma, edge_length);
  }

private:
  double edge_length;
};

} // namespace mmd

#endif

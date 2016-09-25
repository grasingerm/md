#ifndef __INTEGRATION_HPP__
#define __INTEGRATION_HPP__

#include <vector>
#include <functional>
#include <armadillo>
#include "potentials.hpp"

namespace mmd {

using time_integrator = 
  std::function<void(std::vector<abstract_potential*>&, 
                     std::vector<molecular_id>&, arma::mat&, arma::mat&,
                     arma::mat&, const double)>;

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
 */
void euler(std::vector<abstract_potential*>& potentials, 
           std::vector<molecular_id>& molecular_ids, 
           arma::mat& positions, arma::mat& velocities, 
           arma::mat& forces, const double dt);

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
 */
void velocity_verlet(std::vector<abstract_potential*>& potentials, 
                     std::vector<molecular_id>& molecular_ids, 
                     arma::mat& positions, arma::mat& velocities, 
                     arma::mat& forces, const double dt);

} // namespace mmd

#endif

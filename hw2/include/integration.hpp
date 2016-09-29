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

} // namespace mmd

#endif

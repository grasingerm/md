#ifndef __SIMULATION_HPP__
#define __SIMULATION_HPP__

#include <vector>
#include <armadillo>
#include "integration.hpp"
#include "potentials.hpp"
#include "callbacks.hpp"

namespace mmd {

/*! \brief Simulation object that tracks molecular degrees of freedom
 *
 * Simulation object that tracks positions, velcoities, and forces on molecules
 * throughut time.
 * Requires: dt > 0
 */
class simulation {

private: static const time_integrator default_time_int = velocity_verlet;

public:

  /*! \brief Simulation class with all molecules of the same type
   *
   * \param   id      Id of all molecules in the simulation
   * \param   fname   Filename that contains positions of molecules
   * \param   pot     Potential function that acts on molecules
   * \param   dt      Time step size
   * \param   vscale  Scale for random velocities that are generated 
   * \return        Molecular simulation object
   */
  simulation(const molecular_id id, const char* fname, 
             const abstract_potential* pot, const double dt, 
             const double vscale = 1.0);

  ~simulation() {}

  /*! \brief Simulate a number of time steps
   *
   * \param   nsteps  Number of time steps to simulate
   */
  simulate(const unsigned nsteps);

  /*! \brief Add a callback function
   *
   * \param   cb    Callback function to add
   */
  inline void add_callback(const callback cb) { callbacks.push_back(cb); }
  
  /*! \brief Get the current molecular positions
   *
   * \return  Current molecular positions
   */
  inline const std::vector<molecular_id>& get_molecular_ids() const { 
    return molecular_ids;
  }
  
  /*! \brief Get the simulation mass accessor
   *
   * \return  Mass accessor
   */
  inline const mass_accessor& get_mass_accessor() const { 
    return ma;
  }

  /*! \brief Get the current molecular positions
   *
   * \return  Current molecular positions
   */
  inline const arma::mat& get_positions() const { return positions; }

  /*! \brief Get the current molecular velocities
   *
   * \return  Current molecular velocities
   */
  inline const arma::mat& get_velocities() const { return velocities; }

  /*! \brief Get the current intermolecular forces
   *
   * \return  Current intermolecular forces
   */
  inline const arma::mat& get_forces() const { return forces; }
  
  /*! \brief Get the current time
   *
   * \return  Current time
   */
  inline double get_time() const { return t; }

  /*! \brief Get the dt
   *
   * \return  dt
   */
  inline double get_dt() const { return dt; }

private:

  std::vector<molecular_ids> molecular_ids;
  mass_accessor ma;
  arma::mat positions;
  arma::mat velocities;
  arma::mat forces;
  std::vector<abstract_potential*> potentials;
  time_integrator time_int;
  double dt;
  double t;
  std::vector<callback> callbacks;
};

/*! \brief Calculate the total potential energy for a simulation
 *
 * \param   sim   Simulation object
 * \return        Potential energy
 */
inline double potential_energy(const simulation& sim) {
  double sum = 0.0;
  for (const auto& potential : potentials)
    sum += potential(sim.get_molecular_ids(), sim.get_positions());
  return sum;
}

/*! \brief Calculate the total kinetic energy for a simulation
 *
 * \param   sim   Simulation object
 * \return        Kinetic energy
 */
double kinetic_energy(const simulation& sim);

/*! \brief Calculate the total energy for a simulation
 *
 * \param   sim   Simulation object
 * \return        Energy
 */
inline double total_energy(const simulation& sim) {
  return potential_energy(sim) + total_energy(sim);
}

/*! \brief Calculate the total momentum for a simulation
 *
 * \param   sim   Simulation object
 * \return        Total momentum
 */
arma::vec momentum(const simulation& sim);

} // namespace mmd

#endif

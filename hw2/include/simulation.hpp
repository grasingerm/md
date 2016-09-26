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

private:

  std::vector<molecular_ids> molecular_ids;
  arma::mat positions;
  arma::mat velocities;
  arma::mat forces;
  std::vector<abstract_potential*> potentials;
  time_integrator time_int;
  double dt;
  std::vector<callback> callbacks;
};

} // namespace mmd

#endif

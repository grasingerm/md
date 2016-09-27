#ifndef __CALLBACKS_HPP__
#define __CALLBACKS_HPP__

#include <initializer_list>
#include <vector>
#include <functional>
#include <fstream>
#include <armadillo>
#include "simulation.hpp"

namespace mmd {

using callback = std::function<void(const simulation&)>;
using data_accessor = std::function<const arma::mat&(const simulation&)>;
using value_accessor = std::function<const double(const simulation&)>;

/*! \brief Wrapper for accessing molecular positions
 *
 * \param   sim   Simulation object
 * \return        Molecular positions
 */
inline const arma::mat& positions(const simulation& sim) { 
  return sim.get_positions(); 
}

/*! \brief Wrapper for accessing molecular velocities
 *
 * \param   sim   Simulation object
 * \return        Molecular velocities
 */
inline const arma::mat& velocities(const simulation& sim) { 
  return sim.get_velocities(); 
}

/*! \brief Wrapper for accessing intermolecular forces
 *
 * \param   sim   Simulation object
 * \return        Molecular forces
 */
inline const arma::mat& forces(const simulation& sim) { 
  return sim.get_forces(); 
}

/*! \brief Callback for saving xyz data to file
 */
class save_xyz_callback {
public:
  /*! \brief Constructor for callback function that saves data in xyz format
   *
   * \param   fname   File name of the output file
   * \param   dt      Frequency with which to write data
   * \param   da      Function for accessing simulation data
   * \return          Callback function
   */
  save_xyz_callback(const char* fname, const double dt, data_accessor da) 
    : outfile(fname), dt(dt), da(da) {}
  
  ~save_xyz_callback() { outfile.close(); }
  
  operator()(const simulation&);
private:
  ofstream outfile;
  double dt;
  data_accessor da;
};

/*! \brief Callback for saving data values to a delimited file
 */
class save_values_with_time_callback {
public:
  /*! \brief Constructor for callback function that saves delimited data
   *
   * \param   fname   File name of the output file
   * \param   dt      Frequency with which to write data
   * \param   vas     Functions for accessing simulation values
   * \param   delim   Character delimiter
   * \return          Callback function
   */
  save_values_with_time_callback(const char* fname, const double dt,
      const initializer_list<value_accessor>& vas, const char delim=',') 
    : outfile(fname), dt(dt), vas(vas), delim(delim) {}
  
  ~save_values_with_time_callback() { outfile.close(); }
  
  operator()(const simulation&);

  /*! Static factory for constructing an energy and momentum saving callback
   *
   * \param   fname   File name of the output file
   * \param   dt      Frequency with which to write data
   * \param   delim   Character delimiter
   * \return          Callback function
   */
  inline static save_energy_and_momentum_with_time_callback
    (const char* fname, const double dt, const char delim=',') {

    return save_values_with_time_callback(fname, dt, 
        {potential_energy, kinetic_energy, total_energy, momentum}, delim);

  }

private:
  ofstream outfile;
  double dt;
  std::vector<value_accessor> vas;
  char delim;
};

/*! \brief Check conservation of energy
 *
 * \param   dt  Frequency with which to check energy
 * \param   eps Tolerance of check
 * \return      Callback
 */
callback check_energy(const double dt, const double eps=1e-7);

/*! \brief Check conservation of momentum
 *
 * \param   dt  Frequency with which to check momentum
 * \param   eps Tolerance of check
 * \return      Callback
 */
callback check_momentum(const double dt, const double eps=1e-7);

} // namespace mmd

#endif

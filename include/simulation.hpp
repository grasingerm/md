#ifndef __SIMULATION_HPP__
#define __SIMULATION_HPP__

#include "integration.hpp"
#include "mmd_types.hpp"
#include "potentials.hpp"
#include <armadillo>
#include <array>
#include <vector>

/* TODO: using default_time_int as a default param might be cleaner than all of
 *       the many overloads
 *
 * TODO: uhhh, there are 1000 overloads for the simulation constructor, can we
 *       simplify this nonsense?
 */

namespace mmd {

/*! \brief Simulation object that tracks molecular degrees of freedom
 *
 * Simulation object that tracks positions, velcoities, and forces on molecules
 * throughut time.
 * Requires: dt > 0
 */
class simulation {

private:
  static const time_integrator default_time_int;
  static constexpr double default_kB = 1.0;

public:
  /*! \brief Simulation class with all molecules of the same type
   *
   * \param   id            Id of all molecules in the simulation
   * \param   fname         Filename that contains positions of molecules
   * \param   pot           Potential function that acts on molecules
   * \param   dt            Time step size
   * \param   edge_length   Edge length of the control volume
   * \return                Molecular simulation object
   */
  simulation(const molecular_id id, const char *fname, abstract_potential *pot,
             const double dt, const double edge_length);

  /*! \brief Simulation class with all molecules of the same type
   *
   * \param   id            Id of all molecules in the simulation
   * \param   fname         Filename that contains positions of molecules
   * \param   pot           Potential function that acts on molecules
   * \param   dt            Time step size
   * \param   ti            Time integrator
   * \param   edge_length   Edge length of the control volume
   * \return                Molecular simulation object
   */
  simulation(const molecular_id id, const char *fname, abstract_potential *pot,
             const double dt, const time_integrator &ti,
             const double edge_length)
      : simulation(id, fname, pot, dt, edge_length) {
    time_int = ti;
  }

  /*! \brief Simulation class with all molecules of the same type
   *
   * \param   ids           Id of all molecules in the simulation
   * \param   ma            Mass accessor
   * \param   fname         Filename that contains positions of molecules
   * \param   pot           Potential function that acts on molecules
   * \param   dt            Time step size
   * \param   ti            Time integrator
   * \param   edge_length   Edge length of the control volume
   * \return                Molecular simulation object
   */
  simulation(const std::vector<molecular_id>& ids, const mass_accessor& ma,
             const char *fname, 
             abstract_potential *pot, const double dt, 
             const time_integrator &ti,
             const double edge_length);

  /*! \brief Simulation class with all molecules of the same type
   *
   * \param   ids           Id of all molecules in the simulation
   * \param   ma            Mass accessor
   * \param   fname         Filename that contains positions of molecules
   * \param   pot           Potential function that acts on molecules
   * \param   dt            Time step size
   * \param   ti            Time integrator
   * \param   tstar         Initial temperature
   * \param   edge_length   Edge length of the control volume
   * \return                Molecular simulation object
   */
  simulation(const std::vector<molecular_id>& ids, const mass_accessor& ma,
             const char *fname, 
             abstract_potential *pot, const double dt, 
             const time_integrator &ti, const double tstar,
             const double edge_length);

  /*! \brief Simulation class with all molecules of the same type
   *
   * \param   id            Id of all molecules in the simulation
   * \param   fname         Filename that contains positions of molecules
   * \param   pot           Potential function that acts on molecules
   * \param   dt            Time step size
   * \param   tstar         Non-dimensional temperature
   * \param   edge_length   Edge length of the control volume
   * \return                Molecular simulation object
   */
  simulation(const molecular_id id, const char *fname, abstract_potential *pot,
             const double dt, const double tstar, const double edge_length);

  /*! \brief Simulation class with all molecules of the same type
   *
   * \param   id            Id of all molecules in the simulation
   * \param   fname         Filename that contains positions of molecules
   * \param   pot           Potential function that acts on molecules
   * \param   dt            Time step size
   * \param   tstar         Non-dimensional temperature
   * \param   ti            Time integrator
   * \param   edge_length   Edge length of the control volume
   * \return                Molecular simulation object
   */
  simulation(const molecular_id id, const char *fname, abstract_potential *pot,
             const double dt, const double tstar, const time_integrator &ti,
             const double edge_length)
      : simulation(id, fname, pot, dt, tstar, edge_length) {
    time_int = ti;
  }

  /*! \brief Simulation class with all molecules of the same type
   *
   * \param   id            Id of all molecules in the simulation
   * \param   fname         Filename that contains positions of molecules
   * \param   density       Non-dimensional density
   * \param   pot           Potential function that acts on molecules
   * \param   dt            Time step size
   * \param   tstar         Non-dimensional temperature
   * \param   edge_length   Edge length of the control volume
   * \return                Molecular simulation object
   */
  simulation(const molecular_id id, const char *fname, const double density,
             abstract_potential *pot, const double dt, const double tstar,
             const double edge_length);

  /*! \brief Simulation class with all molecules of the same type
   *
   * \param   id            Id of all molecules in the simulation
   * \param   fname         Filename that contains positions of molecules
   * \param   density       Non-dimensional density
   * \param   pot           Potential function that acts on molecules
   * \param   dt            Time step size
   * \param   tstar         Non-dimensional temperature
   * \param   ti            Time integrator
   * \param   edge_length   Edge length of the control volume
   * \return                Molecular simulation object
   */
  simulation(const molecular_id id, const char *fname, const double density,
             abstract_potential *pot, const double dt, const double tstar,
             const time_integrator &ti, const double edge_length)
      : simulation(id, fname, density, pot, dt, tstar, edge_length) {
    time_int = ti;
  }

  /*! \brief Simulation class with all molecules of the same type
   *
   * \param   id            Id of all molecules in the simulation
   * \param   fname_pos     Filename that contains positions of molecules
   * \param   fname_vel     Filename that contains positions of molecules
   * \param   pot           Potential function that acts on molecules
   * \param   dt            Time step size
   * \param   vscale        Scale for random velocities that are generated
   * \param   edge_length   Edge length of the control volume
   * \return                Molecular simulation object
   */
  simulation(const molecular_id id, const char *fname_pos,
             const char *fname_vel, abstract_potential *pot, const double dt,
             const double edge_length);

  /*! \brief Simulation class with all molecules of the same type
   *
   * \param   id            Id of all molecules in the simulation
   * \param   fname_pos     Filename that contains positions of molecules
   * \param   fname_vel     Filename that contains positions of molecules
   * \param   pot           Potential function that acts on molecules
   * \param   dt            Time step size
   * \param   vscale        Scale for random velocities that are generated
   * \param   ti            Time integrator
   * \param   edge_length   Edge length of the control volume
   * \return                Molecular simulation object
   */
  simulation(const molecular_id id, const char *fname_pos,
             const char *fname_vel, abstract_potential *pot, const double dt,
             const time_integrator &ti, const double edge_length)
      : simulation(id, fname_pos, fname_vel, pot, dt, edge_length) {
    time_int = ti;
  }

  /*! \brief Simulation class with all molecules of the same type
   *
   * \param   id            Id of all molecules in the simulation
   * \param   density       Non-dimensional density
   * \param   pot           Potential function that acts on molecules
   * \param   dt            Time step size
   * \param   tstar         Non-dimensional temperature
   * \param   edge_length   Edge length of the control volume
   * \return                Molecular simulation object
   */
  simulation(const molecular_id id, const double density,
             abstract_potential *pot, const double dt, const double tstar,
             const double edge_length);

  /*! \brief Simulation class with all molecules of the same type
   *
   * \param   id            Id of all molecules in the simulation
   * \param   density       Non-dimensional density
   * \param   pot           Potential function that acts on molecules
   * \param   dt            Time step size
   * \param   tstar         Non-dimensional temperature
   * \param   ti            Time integrator
   * \param   edge_length   Edge length of the control volume
   * \return                Molecular simulation object
   */
  simulation(const molecular_id id, const double density,
             abstract_potential *pot, const double dt, const double tstar,
             const time_integrator &ti, const double edge_length)
      : simulation(id, density, pot, dt, tstar, edge_length) {
    time_int = ti;
  }

  ~simulation() {}

  /*! \brief Simulate a number of time steps
   *
   * \param   nsteps  Number of time steps to simulate
   */
  void simulate(const unsigned long long nsteps);

  /*! \brief Add a callback function
   *
   * \param   cb    Callback function to add
   */
  inline void add_callback(const callback cb) { callbacks.push_back(cb); }

  /*! \brief Add a mutator function
   *
   * \param   cb    Mutator function to add
   */
  inline void add_mutator(const mutator mt) { mutators.push_back(mt); }

  /*! \brief Get the current molecular positions
   *
   * \return  Current molecular positions
   */
  inline const std::vector<molecular_id> &get_molecular_ids() const {
    return molecular_ids;
  }

  /*! \brief Get the number of molecules
   *
   * \return  N
   */
  inline size_t get_N() const { return molecular_ids.size(); }

  /*! \brief Get the simulation mass accessor
   *
   * \return  Mass accessor
   */
  inline const mass_accessor &get_mass_accessor() const { return ma; }

  /*! \brief Get the current molecular positions
   *
   * \return  Current molecular positions
   */
  inline const arma::mat &get_positions() const { return positions; }

  /*! \brief Get the current molecular positions
   *
   * \return  Current molecular positions
   */
  inline arma::mat &get_positions() { return positions; }

  /*! \brief Get the current molecular velocities
   *
   * \return  Current molecular velocities
   */
  inline const arma::mat &get_velocities() const { return velocities; }

  /*! \brief Get the current molecular velocities
   *
   * \return  Current molecular velocities
   */
  inline arma::mat &get_velocities() { return velocities; }

  /*! \brief Get the current intermolecular forces
   *
   * \return  Current intermolecular forces
   */
  inline const arma::mat &get_forces() const { return forces; }

  /*! \brief Get the current intermolecular forces
   *
   * \return  Current intermolecular forces
   */
  inline arma::mat &get_forces() { return forces; }

  /*! \brief Get the current intermolecular forces
   *
   * \return  Current intermolecular forces
   */
  inline const std::vector<abstract_potential *> &get_potentials() const {
    return potentials;
  }

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

  /*! \brief Get the volume
   *
   * \return  Volume
   */
  inline double get_edge_length() const { return edge_length; }

  /*! \brief Get the volume
   *
   * \return  Volume
   */
  inline double get_volume() const { return volume; }

  /*! \brief Get the simulation Boltzmann's constant
   *
   * \return  Boltzmann's constant
   */
  inline double get_kB() const { return kB; }

  /*! \brief Reset the clock
   */
  inline void reset_clock() { t = 0.0; }

  /*! \brief Clear callbacks
   */
  inline void clear_callbacks() { callbacks.clear(); }

  /*! \brief Clear mutators
   */
  inline void clear_mutators() { mutators.clear(); }

  /*! \brief Clear mutators
   */
  inline void clear_all_callbacks() { clear_callbacks(); clear_mutators(); }

  /*! \brief Set time integrator
   */
  inline void set_time_int(const time_integrator &ti) { time_int = ti; }

private:
  std::vector<molecular_id> molecular_ids;
  mass_accessor ma;
  arma::mat positions;
  arma::mat velocities;
  arma::mat forces;
  std::vector<abstract_potential *> potentials;
  time_integrator time_int;
  double dt;
  double t;
  std::vector<callback> callbacks;
  std::vector<mutator> mutators;
  double edge_length;
  double volume;
  double kB;
};

/*! \brief Calculate the total potential energy for a simulation
 *
 * \param   sim   Simulation object
 * \return        Potential energy
 */
inline double potential_energy(const simulation &sim) {
  double sum = 0.0;
  for (const auto &potential : sim.get_potentials())
    sum += potential->potential_energy(sim.get_molecular_ids(),
                                       sim.get_positions());
  return sum;
}

/*! \brief Calculate the total kinetic energy for a simulation
 *
 * \param   sim   Simulation object
 * \return        Kinetic energy
 */
double kinetic_energy(const simulation &sim);

/*! \brief Calculate the total energy for a simulation
 *
 * \param   sim   Simulation object
 * \return        Energy
 */
inline double total_energy(const simulation &sim) {
  return potential_energy(sim) + kinetic_energy(sim);
}

/*! \brief Calculate the total momentum for a simulation
 *
 * \param   sim   Simulation object
 * \return        Total momentum
 */
arma::vec momentum(const simulation &sim);

/*! \brief Calculate the instantaneous temperature for a simulation
 *
 * \param   sim   Simulation object
 * \return        Temperature
 */
double temperature(const simulation &sim);

/*! \brief Calculate the instantaneous pressure due to ideal gas EOS
 *
 * \param   sim   Simulation object
 * \return        Ideal gas pressure
 */
double ideal_pressure(const simulation &sim);

/*! \brief Calculate the instantaneous virial pressure
 *
 * \param   sim   Simulation object
 * \return        Virial pressure
 */
double virial_pressure(const simulation &sim);

/*! \brief Calculate the instantaneous pressure for a simulation
 *
 * \param   sim   Simulation object
 * \return        Pressure
 */
inline double pressure(const simulation &sim) {
  return ideal_pressure(sim) + virial_pressure(sim);
}

/*! \brief Calculate kinetic energy, temperature and pressure for a simulation
 *
 * Because kinetic energy, temperature, and pressure are all related, this
 * function will be more efficient than calling each calculation separately
 *
 * \param   sim   Simulation object
 * \return        Kinetic energy, temperature, and pressure
 */
std::array<double, 3> ktp(const simulation &sim);

/*! \brief Calculate the energies, temperature and pressure for a simulation
 *
 * Because energies, temperature, and pressure are all related, this
 * function will be more efficient than calling each calculation separately
 *
 * \param   sim   Simulation object
 * \return        K+U, U, K, temperature, and pressure
 */
std::array<double, 5> euktp(const simulation &sim);

/*! \brief Calculate the energies, temperature and pressure for a simulation
 *
 * Because energies, temperature, and pressure are all related, this
 * function will be more efficient than calling each calculation separately
 *
 * \param   sim   Simulation object
 * \return        K+U, U, K, temperature, pressure, ideal gas term, virial term
 */
std::array<double, 7> euktpiv(const simulation &sim);

} // namespace mmd

#endif

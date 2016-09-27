#ifndef __POTENTIALS_HPP__
#define __POTENTIALS_HPP__

#include <armadillo>
#include <vector>
#include <initializer_list>
#include <stdexcept>
#include "debug.hpp"

namespace mmd {

/*! \brief Abstract base class for potentials
 *
 * Defines a public interface for all potentials (e.g. spring, LJ, bond-angle,
 * etc.)
 */
class abstract_potential {
public:
  virtual ~abstract_potential()=0;

  /*! \brief Calculates total potential energy of the current configuration
   *
   * Calculates the total potential energy of the current configuration given
   * each molecule type and molecule position. The molecular positions are
   * stored in a matrix such that the ith column vector is the position vector 
   * of the ith molecule.
   *
   * \param    molecular_ids   Collection of molecular identities
   * \param    positions       Matrix of molecular positions 
   * \return                   Potential energy of the current configuration
   */
  inline double potential_energy(const std::vector<molecular_id>& molecular_ids, 
                                 const arma::mat& positions) const {
    _check_arg_sizes(molecular_ids, positions);
    return _potential_energy(molecular_ids, positions);
  }
  
  /*! \brief Calculates total potential energy of the current configuration
   * 
   * Calculates the total potential energy of the current configuration given
   * each molecule type and molecule position. The molecular positions (and 
   * forces) are stored in a matrix such that the ith column vector is the 
   * position (force) vector of the ith molecule. This function is a mutator
   * that increments the matrix of forces.
   *
   * \param    molecular_ids   Collection of molecular identities
   * \param    positions       Collection of molecular positions
   * \param    forces          Collection of force vectors
   */
  inline void increment_forces(const std::vector<molecular_id>& molecular_ids, 
                               const arma::mat& positions, arma::mat& forces) 
    const {
    _check_arg_sizes(molecular_ids, positions, forces);
    _increment_forces(molecular_ids, positions, forces);
  }
private:
  virtual double _potential_energy(const std::vector<molecular_id>&, 
                                   const arma::mat&) const=0;
  virtual void _increment_forces(const std::vector<molecular_id>&, 
                                 const arma::mat&, arma::mat&) const=0;
};

/*! Public interface for a 6-12 Lennard-Jones pairwise potential
 */
class abstract_LJ_potential : public abstract_potential {
public:
  virtual ~abstract_LJ_potential()=0;

  /*! \brief Get the well depth of the potential (epsilon)
   *
   * \param    id1   Molecular id of the first molecule
   * \param    id2   Molecular id of the second molecule
   * \return         Depth of the potential well
   */
  inline double get_well_depth(const molecular_id id1, const molecular_id id2) 
    const { 
    _get_well_depth(id1, id2); 
  }
  
  /*! \brief Get the distance at which the potential is zero (sigma)
   *
   * \param    id1   Molecular id of the first molecule
   * \param    id2   Molecular id of the second molecule
   * \return         Zero of the potential
   */
  inline double get_rzero (const molecular_id id1, const molecular_id id2) 
    const { 
    _get_rzero(id1, id2); 
  }

private:
  virtual double _potential_energy(const std::vector<molecular_id>&, 
                                   const arma::mat&, arma::mat&) const;
  virtual void _increment_forces(const std::vector<molecular_id>&, 
                                 const arma::mat&, arma::mat&) const;
  virtual double _get_well_depth(const molecular_id, const molecular_id) const=0;
  virtual double _get_rzero(const molecular_id, const molecular_id) const=0;
};

/*! \brief 6-12 Lennard-Jones pairwise potential
 *
 * Lennard-Jones potential to be used in a simulation where the potential
 * function between every pair of molecules is the same (i.e. to be used in
 * a simulation where all of the molecules are of the same type).
 */
class const_well_params_LJ_potential : public abstract_potential {
public:
  /*! \brief Constructor for a LJ potential with constant well parameters
   *
   * \param    well_depth    Depth of the potential well
   * \param    rzero         Finite distance at which potential is zero
   * \return                 LJ potential
   */
  const_well_params_LJ_potential(const double well_depth, const double rzero)
    : _well_depth(well_depth), _rzero(rzero) {}

  ~const_well_params_LJ_potential() {}

private:
  virtual double _get_well_depth(molecular_id, molecular_id) const { 
    return _well_depth;
  }

  virtual double _get_rzero(molecular_id, molecular_id) const {
    return _rzero;
  }

  double _well_depth;
  double _rzero;
};

/*! \brief Public interface for a spring potential
 */
class abstract_spring_potential : public abstract_potential {
public:
  virtual ~abstract_spring_potential()=0;

  /*! \brief Get the spring constant
   *
   * \param    id    Molecular id of the molecule
   * \return         Spring constant
   */
  inline double get_k(molecular_id id) const { 
    _get_k(id); 
  }

private:
  virtual double _potential_energy(const std::vector<molecular_id>&, arma::mat&, 
                                   arma::mat&) const;
  virtual void _increment_forces(const std::vector<molecular_id>&, arma::mat&, 
                                 arma::mat&) const;
  virtual double _get_k(molecular_id) const=0;
};

/*! \brief Potential due to a spring with constant k
 *
 * Potential due to a spring in which the spring constant, k, is the same
 * for all molecules, regardless of molecular id.
 * k should always be greater than 0.
 */
class const_k_spring_potential : public abstract_spring_potential {
public:
  /*! \brief Constructor for potential with same spring constant for all molecules
   *
   * Constructor for potential with same spring constant for all molecules.
   * Requires: k > 0
   *
   * \param    k     Spring constant
   * \return         Spring potential
   */
  const_k_spring_potential(const double k) : _k(k) {}

  virtual ~const_k_spring_potential() {}

private:
  virtual double _get_k(molecular_id) const { return _k };

  double _k;
};

/*! \brief Potential due to a spring with a polynomial potential
 *
 * Potential due to a spring in which the potential, a x^n + b x^-1 + ..., 
 * is the same for all molecules, regardless of molecular id.
 */
class const_poly_spring_potential : public abstract_potential {
public:
  /*! \brief Constructor for potential with same spring potential for all molecules
   *
   * Constructor for potential with same spring constant for all molecules.
   *
   * \param    coeffs   Polynomial coefficients
   * \return            Spring potential
   */
  const_poly_spring_potential(const std::initializer_list<double>& coeffs); 

  virtual ~const_poly_spring_potential() {}

private:
  virtual double _potential_energy(const std::vector<molecular_id>&, arma::mat&, 
                                   arma::mat&) const;
  virtual void _increment_forces(const std::vector<molecular_id>&, arma::mat&, 
                                 arma::mat&) const;

  std::vector<double> pcoeffs;
  std::vector<double> fcoeffs;
};

} // namespace mmd

#endif

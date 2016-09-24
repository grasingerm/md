#ifndef __POTENTIALS_HPP__
#define __POTENTIALS_HPP__

#include <armadillo>
#include <vector>

namespace mmd {

//! Abstract base class for potentials
//!
//! Defines a public interface for all potentials (e.g. spring, LJ, bond-angle,
//! etc.)
class abstract_potential {
public:
  //! Calculates total potential energy of the current configuration
  //!
  //! 
  inline double potential_energy(const std::vector<molecular_id>& molecular_ids, 
                                 arma::mat& positions) {
    return _potential_energy(molecular_ids, positions);
  }
  
  //! Calculates total potential energy of the current configuration
  //!
  //! 
  inline void update_forces(const std::vector<molecular_id>& molecular_ids, 
                            arma::mat& positions, arma::mat& forces) {
    _update_forces(molecular_ids, positions, forces);
  }
private:
  virtual double _potential_energy(const std::vector&, arma::mat&, arma::mat&)=0;
  virtual void _update_forces(const std::vector&, arma::mat&, arma::mat&)=0;
};

}

#endif

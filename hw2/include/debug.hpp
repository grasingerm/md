#ifndef __DEBUG_HPP__
#define __DEBUG_HPP__

#include <vector>
#include <armadillo>
#include "molecular.hpp"

namespace mmd {

inline void _check_arg_sizes(const std::vector<molecular_id>& molecular_ids, 
                             const arma::mat& positions) {
  #ifndef MMD_DISABLE_RUNTIME_CHECKS
  if (molecular_ids.size() != positions.n_cols)
    throw std::invalid_argument("The size of molecular ids and the number of "
                                "columns in the position matrix must be equal");
  #endif
  return;
}

inline void _check_arg_sizes(const std::vector<molecular_id>& molecular_ids, 
                             const arma::mat& positions, 
                             const arma::mat& forces) { 
  #ifndef MMD_DISABLE_RUNTIME_CHECKS
  if (molecular_ids.size() != positions.n_cols || 
      molecular_ids.size() != forces.n_cols)
    throw std::invalid_argument("The size of molecular ids, the number of "
                                "columns in the position matrix, and the "
                                "number of columns in the forces matrix "
                                "must be equal");
  #endif
  return;
}

inline void _check_arg_sizes(const std::vector<molecular_id>& molecular_ids, 
                             const arma::mat& positions, 
                             const arma::mat& velocities,
                             const arma::mat& forces) { 
  #ifndef MMD_DISABLE_RUNTIME_CHECKS
  if (molecular_ids.size() != positions.n_cols || 
      molecular_ids.size() != velocities.n_cols ||
      molecular_ids.size() != forces.n_cols)
    throw std::invalid_argument("The size of molecular ids, the number of "
                                "columns in the position matrix, the "
                                "number of columns in the forces matrix, and "
                                "the number of columns in the velocities "
                                "must all be equal");
  #endif
  return;
}

} // namespace mmd

#endif

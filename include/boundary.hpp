#ifndef __BOUNDARY_HPP__
#define __BOUNDARY_HPP__

#include <armadillo>

namespace mmd {

void apply_pbc(arma::mat &positions, const double edge_length, const size_t i);

inline void apply_pbcs(arma::mat &positions, const double edge_length) {
  /* TODO: parallelize this... */
  for (auto i = size_t{0}; i < positions.n_cols; ++i)
    apply_pbc(positions, edge_length, i);
}

} // namespace mmd

#endif

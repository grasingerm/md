#include "distance.hpp"

namespace mmd {

using namespace arma;
  
double rij2(const mat& positions, const size_t i, const size_t j) {
  const double *ri = positions.colptr(i), *rj = positions.colptr(j);
  double rij2 = 0;
  for (auto k = size_t{0}; k < positions.n_rows; ++k) {
    double rijk = ri[k] - rj[k];
    rij2 += rijk * rijk;
  }
  return rij2;
}

double rij2_pbc(const arma::mat& positions, const size_t i, const size_t j,
                const double edge_length) {
  const double *ri = positions.colptr(i), *rj = positions.colptr(j);
  double rij2 = 0;
  for (auto k = size_t{0}; k < positions.n_rows; ++k) {
    double rijk = ri[k] - rj[k];

    /* correct for periodic boundary conditions using nearest image convention 
     */
    if (rijk > edge_length / 2.0) rijk -= edge_length;
    else if (rijk < -edge_length / 2.0) rijk += edge_length;

    rij2 += rijk * rijk;
  }
  return rij2;
}

arma::vec rij_pbc(const arma::mat& positions, const size_t i, const size_t j,
                  const double edge_length) {

  arma::vec rij = positions.col(i) - positions.col(j);

  /* correct for periodic boundary conditions using nearest image convention 
   */
  if (rij(0) > edge_length / 2.0) rij(0) -= edge_length;
  else if (rij(0) < -edge_length / 2.0) rij(0) += edge_length;

  if (rij(1) > edge_length / 2.0) rij(1) -= edge_length;
  else if (rij(1) < -edge_length / 2.0) rij(1) += edge_length;

  if (rij(2) > edge_length / 2.0) rij(2) -= edge_length;
  else if (rij(2) < -edge_length / 2.0) rij(2) += edge_length;

  return rij;
}

} // namespace mmd

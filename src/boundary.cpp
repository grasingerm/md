#include "boundary.hpp"

namespace mmd {

using namespace arma;

void apply_pbc(mat &positions, const double edge_length, const size_t i) {
  /* unrolled loop over each dimension */
  if (positions(0, i) > edge_length) {
    positions(0, i) -= edge_length;
  }
  else if (positions(0, i) < 0.0) {
    positions(0, i) += edge_length;
  }

  if (positions(1, i) > edge_length) {
    positions(1, i) -= edge_length;
  }
  else if (positions(1, i) < 0.0) {
    positions(1, i) += edge_length;
  }

  if (positions(2, i) > edge_length) {
    positions(2, i) -= edge_length;
  }
  else if (positions(2, i) < 0.0) {
    positions(2, i) += edge_length;
  }
}

} // namespace mmd

#ifndef __DISTANCE_HPP__
#define __DISTANCE_HPP__

#include <armadillo>

namespace mmd {

/*! \brief Get the square of the distance between molecules i and j
 *
 * Calculates the square of the distance between molecules i and j,
 * Requires: 0 <= i <= positions.n_cols, 0 <= j <= positions.n_cols
 * Guarantees: >= 0
 *
 * \param     positions     Molecular positions
 * \param     i             Index of molecule
 * \param     j             Index of molecule
 * \return                  Distance between the molecules squared
 */
double rij2(const arma::mat &positions, const size_t i, const size_t j);

/*! \brief Get the square of the distance between molecules i and j with PBC
 *
 * Get the square of the distance between molecules i and j using periodic
 * boundary conditions and the nearest image convention.
 * Requires:
 *  * 0 <= i <= positions.n_cols
 *  * 0 <= j <= positions.n_cols
 *  * edge_length > 0
 * Guarantees: >= 0
 *
 * \param     positions     Molecular positions
 * \param     i             Index of molecule
 * \param     j             Index of molecule
 * \param     edge_length   Edge length of control volume
 * \return                  Distance between the molecules squared
 */
double rij2_pbc(const arma::mat &positions, const size_t i, const size_t j,
                const double edge_length);

/*! \brief Get the relative position between molecules i and j with PBC
 *
 * Get the relative position between molecules i and j using periodic
 * boundary conditions and the nearest image convention.
 * Requires:
 *  * 0 <= i <= positions.n_cols
 *  * 0 <= j <= positions.n_cols
 *  * edge_length > 0
 *
 * \param     positions     Molecular positions
 * \param     i             Index of molecule
 * \param     j             Index of molecule
 * \param     edge_length   Edge length of control volume
 * \return                  Relative position between the molecules
 */
arma::vec rij_pbc(const arma::mat &positions, const size_t i, const size_t j,
                  const double edge_length);

} // namespace mmd

#endif

#include "potentials.hpp"

namespace mmd {
  
using namespace arma;

double abstract_LJ_potential::_potential_energy
  (const std::vector<molecular_id>& molecular_ids, const arma::mat& positions) 
  const {
  
  const auto n = molecular_ids.size();
  double potential = 0;
  
  for (auto i = size_t{0}; i < n-1; ++i) {
    for (auto i = j+1; j < n; ++j) {

      const double* ri = positions.colptr(i), rj = positions.colptr(j);
      double rij2 = 0;
      for (auto k = size_t{0}; k < positions.n_rows; ++k)
        rij2 += (ri[k] - rj[k]) * (ri[k] - rj[k]);

      const double rzero = get_rzero(molecular_ids[i], molecular_ids[j]);
      const double rat2 = (rzero * rzero) / rij2;
      const double rat6 = rat2 * rat2 * rat2;
      
      potential += 4.0 * get_well_depth(molecular_ids[i], molecular_ids[j]) *
                   (rat6*rat6 - rat6);

    }
  }

  return potential;
}

void abstract_LJ_potential::_increment_forces
  (const std::vector<molecular_id>& molecular_ids, const arma::mat& positions, 
   arma::mat& forces) const {
  
  const auto n = molecular_ids.size();
  
  for (auto i = size_t{0}; i < n-1; ++i) {
    for (auto i = j+1; j < n; ++j) {

      const double* ri = positions.colptr(i), rj = positions.colptr(j);
      double rij2 = 0;
      for (auto k = size_t{0}; k < positions.n_rows; ++k)
        rij2 += (ri[k] - rj[k]) * (ri[k] - rj[k]);

      const double rzero = get_rzero(molecular_ids[i], molecular_ids[j]);
      const double rat2 = (rzero * rzero) / rij2;
      const double rat4 = rat2 * rat2;
      const double rat8 = rat4 * rat4;
      const double well_depth = get_well_depth(molecular_ids[i], 
                                               molecular_ids[j]);
     
      const auto fij = positions.col(i) * well_depth * 
                       (48 * rat8*rat6 - 24 * rat8);
      forces.col(i) += fij;
      forces.col(j) -= fij; // fji = -fij

    }
  }

}

double abstract_spring_potential::_potential_energy
  (const std::vector<molecular_id>& molecular_ids, arma::mat& positions) const {
  
  const auto n = molecular_ids.size();
  double potential = 0;
  
  for (auto i = size_t{0}; i < n; ++i) {

    const double r2 = dot(positions.col(i), positions.col(i));
    potential += 0.5 * get_k() * r2;

  }

  return potential;
}

void abstract_spring_potential::_increment_forces
  (const std::vector<molecular_id>& molecular_ids, const arma::mat& positions, 
   arma::mat& forces) const {
  
  const auto n = molecular_ids.size();
  
  for (auto i = size_t{0}; i < n; ++i)
    forces.col(i) -= get_k() * positions.col(i);

}

const_poly_spring_potential(const initializer_list<double>& coeffs) 
  : pcoeffs(coeffs) {
  
  for (size_t i = 1; i < pcoeffs.size(); ++i)
    fcoeffs.push_back(pcoeffs[i]*i);

}

double const_poly_spring_potential::_potential_energy
  (const std::vector<molecular_id>& molecular_ids, arma::mat& positions) const {
  
  const auto n = molecular_ids.size();
  double potential = 0;
  
  for (auto i = size_t{0}; i < n; ++i) {

    const double r2 = dot(positions.col(i), positions.col(i));
    for (auto p = size_t{0}; p < pcoeffs.size(); ++p) {
      double x = 1;
      for (auto k = size_t{0}; k < p; ++k) x *= r2;
      potential += pcoeffs[p] * x;
    }

  }

  return potential;
}

void abstract_spring_potential::_increment_forces
  (const std::vector<molecular_id>& molecular_ids, const arma::mat& positions, 
   arma::mat& forces) const {
  
  const auto n = molecular_ids.size();
  
  for (auto i = size_t{0}; i < n; ++i) {
    const double r2 = dot(positions.col(i), positions.col(i));
    for (auto p = size_t{0}; p < pcoeffs.size(); ++p) {
      double x = 1;
      for (auto k = size_t{0}; k < p; ++k) x *= r2;
      forces.col(i) += fcoeffs[p] * x * positions.col(i);
    }

  }

}

} // namespace mmd

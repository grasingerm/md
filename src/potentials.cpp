#include "potentials.hpp"

namespace mmd {
  
using namespace arma;

double abstract_LJ_potential::_potential_energy
  (const std::vector<molecular_id>& molecular_ids, const arma::mat& positions) 
  const {
  
  const auto n = molecular_ids.size();
  double potential = 0;
  
  #pragma omp parallel for private(i, j, ri, rj, rij2, rijk, rzero, rat2, rat6) \
    reduction(+:potential) schedule(dynamic)
  for (auto i = size_t{0}; i < n-1; ++i) {
    for (auto j = i+1; j < n; ++j) {

      const double *ri = positions.colptr(i), *rj = positions.colptr(j);
      double rij2 = 0;
      for (auto k = size_t{0}; k < positions.n_rows; ++k) {
        double rijk = ri[k] - rj[k];
        rij2 += rijk * rijk;
      }
      
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
  
  for (auto i = size_t{0}; i < n; ++i) {
    for (auto j = i+1; j < n; ++j) {

      const double *ri = positions.colptr(i), *rj = positions.colptr(j);
      double rij2 = 0;
      for (auto k = size_t{0}; k < positions.n_rows; ++k) {
        double rijk = ri[k] - rj[k];
        rij2 += rijk * rijk;
      }

      const double rzero = get_rzero(molecular_ids[i], molecular_ids[j]);
      const double rat2 = (rzero * rzero) / rij2;
      const double rat6 = rat2 * rat2 * rat2;
      const double rat8 = rat6 * rat2;
      const double well_depth = get_well_depth(molecular_ids[i], 
                                               molecular_ids[j]);
     
      const auto fij = (positions.col(i) - positions.col(j)) * well_depth * 
                       (48.0 * rat8*rat6 - 24.0 * rat8);
      forces.col(i) += fij;
      forces.col(j) -= fij; // fji = -fij

    }
  }

}

double abstract_spring_potential::_potential_energy
  (const std::vector<molecular_id>& molecular_ids, const arma::mat& positions) const {
  
  const auto n = molecular_ids.size();
  double potential = 0;
  
  #pragma omp parallel for private(i, r2) reduction(+:potential) schedule(dynamic)
  for (auto i = size_t{0}; i < n; ++i) {

    const double r2 = dot(positions.col(i), positions.col(i));
    potential += 0.5 * get_k(molecular_ids[i]) * r2;

  }

  return potential;
}

void abstract_spring_potential::_increment_forces
  (const std::vector<molecular_id>& molecular_ids, const arma::mat& positions, 
   arma::mat& forces) const {
  
  const auto n = molecular_ids.size();
  
  for (auto i = size_t{0}; i < n; ++i)
    forces.col(i) -= get_k(molecular_ids[i]) * positions.col(i);

}

const_poly_spring_potential::const_poly_spring_potential
  (const std::initializer_list<double>& coeffs) : pcoeffs(coeffs) {
 
  for (size_t i = 1; i < pcoeffs.size(); ++i)
    fcoeffs.push_back(pcoeffs[i]*2*i);

}

double const_poly_spring_potential::_potential_energy
  (const std::vector<molecular_id>& molecular_ids, const arma::mat& positions) 
  const {
  
  const auto n = molecular_ids.size();
  double potential = 0;
  
  #pragma omp parallel for private(i, p, r2, x, k) reduction(+:potential)
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

void const_poly_spring_potential::_increment_forces
  (const std::vector<molecular_id>& molecular_ids, const arma::mat& positions, 
   arma::mat& forces) const {
  
  const auto n = molecular_ids.size();
  
  for (auto i = size_t{0}; i < n; ++i) {
    const double r2 = dot(positions.col(i), positions.col(i));
    for (auto p = size_t{0}; p < pcoeffs.size(); ++p) {
      double x = 1;
      for (auto k = size_t{0}; k < p; ++k) x *= r2;
      forces.col(i) -= fcoeffs[p] * x * positions.col(i);
    }

  }

}

double const_quad_spring_potential::_potential_energy
  (const std::vector<molecular_id>& molecular_ids, const arma::mat& positions) 
  const {
  
  const auto n = molecular_ids.size();
  double potential = 0;
  
  #pragma omp parallel for private(i, r2) reduction(+:potential)
  for (auto i = size_t{0}; i < n; ++i) {

    const double r2 = dot(positions.col(i), positions.col(i));
    potential += a * r2 * r2 + b * r2 + c;

  }

  return potential;
}

void const_quad_spring_potential::_increment_forces
  (const std::vector<molecular_id>& molecular_ids, const arma::mat& positions, 
   arma::mat& forces) const {
  
  const auto n = molecular_ids.size();
  
  for (auto i = size_t{0}; i < n; ++i) {

    const double r2 = dot(positions.col(i), positions.col(i));
    forces.col(i) -= (4 * a * r2 + 2 * b) * positions.col(i);

  }

}

// definitions for pure virtual destructors
abstract_potential::~abstract_potential() {}
abstract_LJ_potential::~abstract_LJ_potential() {}
abstract_spring_potential::~abstract_spring_potential() {}

} // namespace mmd

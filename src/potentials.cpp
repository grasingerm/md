#include <cmath>
#include "potentials.hpp"
#include "distance.hpp"

namespace mmd {
  
using namespace arma;

double abstract_LJ_potential::_potential_energy
  (const std::vector<molecular_id>& molecular_ids, const arma::mat& positions) 
  const {
  
  const auto n = molecular_ids.size();
  double potential = 0;
  
  for (auto i = size_t{0}; i < n-1; ++i) {
    for (auto j = i+1; j < n; ++j) {

      const double _rij2 = rij2(positions, i, j);
      const double rzero = get_rzero(molecular_ids[i], molecular_ids[j]);
      const double rat2 = (rzero * rzero) / _rij2;
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
    for (auto j = i+1; j < n; ++j) {
     
      const auto fij = _force_ij(molecular_ids, positions, i, j);
      forces.col(i) += fij;
      forces.col(j) -= fij; // fji = -fij

    }
  }

}

arma::vec abstract_LJ_potential::_force_ij
  (const std::vector<molecular_id>& molecular_ids, const arma::mat& positions, 
   const size_t i, const size_t j) const {
  
  const double _rij2 = rij2(positions, i, j);
  const double rzero = get_rzero(molecular_ids[i], molecular_ids[j]);
  const double rat2 = (rzero * rzero) / _rij2;
  const double rat6 = rat2 * rat2 * rat2;
  const double rat8 = rat6 * rat2;
  const double well_depth = get_well_depth(molecular_ids[i], 
                                           molecular_ids[j]);
 
  return (positions.col(i) - positions.col(j)) * well_depth * 
         (48.0 * rat8*rat6 - 24.0 * rat8);

}

double abstract_LJ_cutoff_potential::_potential_energy
  (const std::vector<molecular_id>& molecular_ids, const arma::mat& positions) 
  const {
  
  const auto n = molecular_ids.size();
  double potential = 0;
  
  for (auto i = size_t{0}; i < n-1; ++i) {
    for (auto j = i+1; j < n; ++j) {

      const double _rij2 = rij2_pbc(positions, i, j, _edge_length);

      if (_rij2 <= _rc2) {
      
        const double rzero = get_rzero(molecular_ids[i], molecular_ids[j]);
        const double rz2 = rzero * rzero;
        const double rz4 = rz2 * rz2;
        const double rz8 = rz4 * rz4;

        const double dudr_rc = 6.0 * (rz4*rz2*rzero) / (_rc7) - 
                               12.0 * (rz8*rz4*rzero) / (_rc13); 
        const double u_rc = ((rz8*rz4) / (_rc12) - (rz4*rz2) / (_rc6)); 

        const double rat2 = rz2 / _rij2;
        const double rat6 = rat2 * rat2 * rat2;
       
        potential += 4.0 * get_well_depth(molecular_ids[i], molecular_ids[j]) *
                     ((rat6*rat6 - rat6) - u_rc - 
                      (std::sqrt(_rij2) - _cutoff) * dudr_rc);

      }

    }
  }

  return potential;
}

void abstract_LJ_cutoff_potential::_increment_forces
  (const std::vector<molecular_id>& molecular_ids, const arma::mat& positions, 
   arma::mat& forces) const {
  
  const auto n = molecular_ids.size();
  
  for (auto i = size_t{0}; i < n-1; ++i) {
    for (auto j = i+1; j < n; ++j) {

      const auto fij = _force_ij(molecular_ids, positions, i, j);
      forces.col(i) += fij;
      forces.col(j) -= fij; // fji = -fij

    }
  }

}

arma::vec abstract_LJ_cutoff_potential::_force_ij
  (const std::vector<molecular_id>& molecular_ids, const arma::mat& positions, 
   const size_t i, const size_t j) const {
  
  const auto _rij = rij_pbc(positions, i, j, _edge_length);
  const double _rij2 = dot(_rij, _rij);

  if (_rij2 <= _rc2) {

    const double rzero = get_rzero(molecular_ids[i], molecular_ids[j]);
    const double rz2 = rzero * rzero;
    const double rz4 = rz2 * rz2;
    const double rz8 = rz4 * rz4;

    const double dudr_rc = 24.0 * (rz4*rz2*rzero) / (_rc7) - 
                           48.0 * (rz8*rz4*rzero) / (_rc13);

    const double rat2 = rz2 / _rij2;
    const double rat6 = rat2 * rat2 * rat2;
    const double rat8 = rat6 * rat2;
    const double well_depth = get_well_depth(molecular_ids[i], 
                                             molecular_ids[j]);
  
    /* TODO: this calculation isn't correct if rzero != 1.0 */
    return _rij * well_depth * 
           (48.0 * rat8*rat6 - 24.0 * rat8 + dudr_rc / std::sqrt(_rij2));
  
  }

  return arma::zeros(3);

}

double abstract_spring_potential::_potential_energy
  (const std::vector<molecular_id>& molecular_ids, const arma::mat& positions) const {
  
  const auto n = molecular_ids.size();
  double potential = 0;
  
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
abstract_LJ_cutoff_potential::~abstract_LJ_cutoff_potential() {}
abstract_spring_potential::~abstract_spring_potential() {}

} // namespace mmd

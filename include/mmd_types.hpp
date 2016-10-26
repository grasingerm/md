#ifndef __MMD_TYPES_HPP__
#define __MMD_TYPES_HPP__

#include "molecular.hpp"
#include <armadillo>

namespace mmd {

class simulation;
class abstract_potential;

using callback = std::function<void(const simulation &)>;
using mutator = std::function<void(simulation &)>;
using data_accessor = std::function<const arma::mat &(const simulation &)>;
using value_accessor = std::function<double(const simulation &)>;

using mass_accessor = std::function<double(molecular_id)>;
using time_integrator = std::function<void(simulation &)>;

} // namespace mmd

/*! \brief Helper function for printing the integer value of an enum class
 *
 * \param   value   Enumeration value
 * \return          Enumeration value converted to integer
 */
template <typename Enumeration>
auto as_integer(Enumeration const value) ->
    typename std::underlying_type<Enumeration>::type {

  return static_cast<typename std::underlying_type<Enumeration>::type>(value);
}

#endif

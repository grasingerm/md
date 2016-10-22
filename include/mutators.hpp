#ifndef __MUTATORS_HPP__
#define __MUTATORS_HPP__

namespace mmd {

/*! \brief Mutating callback that equilibrates initial temperature
 *
 * Mutating callback that equilibrates initial temperature.
 * Requires: eps > 0
 *
 * \param   tstar      Non-dimensional temperature
 * \param   eps        Tolerance for which equilibration is defined
 * \param   nsteps_eq  Number of time steps that temperature should be steady
 * \return             Mutator
 */
mutator equilibrate_temperature(const double tstar, const double eps,
                                const unsigned nsteps_eq = 1);

} // namespace mmd

#endif

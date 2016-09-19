#ifndef __INTEGRATION__
#define __INTEGRATION__

#include "molecule.hh"
#include "potentials.hh"

namespace mmd {

template <size_t D>
using time_integrator = void (*)(molecule&, const array<double, D>&,

template <size_t D>
void velocity_verlet(vector<molecule<D>>&, const vector<∇φs>&, const double Δt);

}

#endif

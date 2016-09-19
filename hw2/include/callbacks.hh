#ifndef __CALLBACKS__
#define __CALLBACKS__

#include "molecule.hh"

namespace mmd {

template <size_t D>
using callback = function<void(const vector<molecule<D>>&, const double)>;

callback save_pos_callback(const char*, const double);

}

#endif

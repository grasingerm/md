#ifndef __CALLBACKS_HPP__
#define __CALLBACKS_HPP__

#include <functional>

namespace mmd {

template <size_t D>
using callback = std::function<void(const vector<molecule<D>>&, const double)>;

callback save_pos_callback(const char*, const double);

}

#endif

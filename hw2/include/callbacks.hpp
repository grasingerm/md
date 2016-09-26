#ifndef __CALLBACKS_HPP__
#define __CALLBACKS_HPP__

#include <functional>

namespace mmd {

using callback = std::function<void(const simulation&)>;

callback save_pos_callback(const char*, const double);

}

#endif

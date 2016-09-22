#ifndef __POTENTIALS_HPP__
#define __POTENTIALS_HPP__

#include <array>
#include <functional>

namespace mmd {

template <size_t D> double pot_spring(const std::array& r);
template <size_t D> double pot_LJ(const std::array& r1, const std::array& r2);

template <size_t D>
class spring_potential {
};

}

#endif

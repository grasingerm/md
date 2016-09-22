#ifndef __SIMULATION_HPP__
#define __SIMULATION_HPP__

#include "potentials.hh"
#include "callbacks.hh"

namespace mmd {

template <size_t D>
class simulation {
private:
  vector<molecule<D>> molecules;
  mass_map masses;
  time_integrator time_int;
  double Δt;
  vector<∇φ_fp1> ∇φ1s;
  vector<∇φ_fp2> ∇φ2s;
  vector<callback> callbacks;

public:
  simulation<D>(vector<molecule<D>>& molecules, time_integrator time_int, 
                const double Δt) : molecules(molecules), time_int(time_int), 
                                   Δt(Δt) {}
  simulation<D>(vector<molecule<D>>&& molecules, time_integrator time_int, 
                const double Δt) : molecules(molecules), time_int(time_int), 
                                   Δt(Δt) {}
  simulation<D>(initializer_list<molecule<D>>& molecules, 
                time_integrator time_int, const double Δt) 
                : molecules(molecules), time_int(time_int), Δt(Δt) {}

  inline void set_Δt(const double Δt) { this->Δt = Δt; }
  inline void add_∇φ_fp1(∇φ_fp1 φ) { ∇φ1s.push(φ); };
  inline void add_∇φ_fp2(∇φ_fp2 φ) { ∇φ2s.push(φ); };
  inline void add_callback(callback cb) { callbacks.push(cb); }

  void simulate(const unsigned); 
};

}

#endif

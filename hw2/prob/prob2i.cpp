#include <iostream>
#include "mmd.hpp"

using namespace mmd;

int main() {
  static const double dt = 1e-2;

  // spring potential, k = 1.0
  const_k_spring_potential spr(1.0);

  simulation sim(molecular_id::Test, "2i_init.xyz", "2i_init.uvw", &spr, dt);

  // check energy is conserved
  sim.add_callback(check_energy(dt, 1e-7));
  
  // check momentum is conserved
  sim.add_callback(check_momentum(dt, 1e-7));
  
  // store energy and momentum
  sim.add_callback(
      save_values_with_time_callback::
      save_energy_and_momentum_with_time_callback("2i_data.csv", dt, ','));
  
  // store mass positions and velocities
  sim.add_callback(save_values_with_time_callback("2i_dofs.csv", dt, 
        {
          [](const simulation& sim) -> double { 
            return sim.get_positions()(0, 0) // get x position of mass
          },
          [](const simulation& sim) -> double { 
            return sim.get_velocities()(0, 0) // get x component of velocity
          }
        }));

  // run simulation
  sim.simulate(static_cast<unsigned>(20.0 / dt));

  return 0;
}

#include <cmath>
#include <iostream>
#include "mmd.hpp"

using namespace std;
using namespace mmd;

int main() {
  static const double dt = 1e-4;

  // spring potential, k = 1.0
  const_k_spring_potential spr(1.0);

  cout << "Initializing simulation...\n\n";
  simulation sim(molecular_id::Test, "2i_init.xyz", "2i_init.uvw", &spr, 100*dt);

  // Print out initial conditions for diagnositc purposes
  cout << "Initial positions...\n";
  cout << "x y z\n";
  for (size_t i = 0; i < sim.get_molecular_ids().size(); ++i)
    for (unsigned j = 0; j < 3; ++j) cout << sim.get_positions()(j, i) << ' ';
  cout << "\n\n";

  cout << "Initial velocities...\n";
  cout << "u v w\n";
  for (size_t i = 0; i < sim.get_molecular_ids().size(); ++i)
    for (unsigned j = 0; j < 3; ++j) cout << sim.get_velocities()(j, i) << ' ';
  cout << "\n\n";

  cout << "Initial forces...\n";
  cout << "Fx Fy Fz\n";
  for (size_t i = 0; i < sim.get_molecular_ids().size(); ++i)
    for (unsigned j = 0; j < 3; ++j) cout << sim.get_forces()(j, i) << ' ';
  cout << "\n\n";

  // check energy is conserved
  sim.add_callback(check_energy(dt, 1e-3));
  
  // check momentum is conserved
  //sim.add_callback(check_momentum(dt, 1e-7));
  
  // store energy and momentum
  sim.add_callback(
      save_energy_and_momentum_with_time_callback("2i_data.csv", 100*dt));
  
  // store mass positions and velocities
  sim.add_callback(save_values_with_time_callback("2i_dofs.csv", 100*dt, 
        {
          [](const simulation& sim) -> double { 
            return sim.get_positions()(0, 0); // get x position of mass
          },
          [](const simulation& sim) -> double { 
            return sim.get_velocities()(0, 0); // get x component of velocity
          }
        }));

  // print energy and momentum
  sim.add_callback(print_energy_and_momentum_with_time_callback(100*dt));

  // run simulation
  cout << "Running simulation...\n";
  sim.simulate(static_cast<unsigned long long>(ceil(20.0 / dt) + 1));

  return 0;
}

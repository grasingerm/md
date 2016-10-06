#include <cmath>
#include <iostream>
#include "mmd.hpp"

using namespace std;
using namespace mmd;

int main() {
  static const double dt = 0.002;

  // dimensionless LJ
  const_well_params_LJ_potential pot(1.0, 1.0);

  cout << "Initializing simulation...\n\n";
  simulation sim(molecular_id::Test, "10.txt", &pot, dt);

  // Print out initial conditions for diagnositc purposes
  cout << "Initial positions...\n";
  cout << "x y z\n";
  for (size_t i = 0; i < sim.get_molecular_ids().size(); ++i) {
    for (unsigned j = 0; j < 3; ++j) cout << sim.get_positions()(j, i) << ' ';
    cout << '\n';
  }
  cout << "\n\n";

  cout << "Initial velocities...\n";
  cout << "u v w\n";
  for (size_t i = 0; i < sim.get_molecular_ids().size(); ++i) {
    for (unsigned j = 0; j < 3; ++j) cout << sim.get_velocities()(j, i) << ' ';
    cout << '\n';
  }
  cout << "\n\n";

  cout << "Initial forces...\n";
  cout << "Fx Fy Fz\n";
  for (size_t i = 0; i < sim.get_molecular_ids().size(); ++i) {
    for (unsigned j = 0; j < 3; ++j) cout << sim.get_forces()(j, i) << ' ';
    cout << '\n';
  }
  cout << "\n\n";

  // check energy is conserved (every time steps)
  sim.add_callback(check_energy(100*dt, 1e-3));
  
  // check momentum is conserved
  // momentum is NOT conserved for this problem...
  sim.add_callback(check_momentum(100*dt, 1e-6));
  
  // store energy and momentum
  sim.add_callback(
      save_energy_and_momentum_with_time_callback("3_data.csv", 2.0 / 100.0));
  
  // store mass positions and velocities
  sim.add_callback(save_values_with_time_callback("3_mom.csv", 2.0 / 100.0, 
        {
          [](const simulation& sim) -> double { 
            return momentum(sim)(0); // get momentum in the x-direction
          },
          [](const simulation& sim) -> double { 
            return momentum(sim)(1); // get momentum in the y-direction
          },
          [](const simulation& sim) -> double { 
            return momentum(sim)(2); // get momentum in the z-direction
          }
        }));

  // save mass position for visualization
  sim.add_callback(save_xyz_callback("3.xyz", 2.0 / 100.0, positions));

  // print energy and momentum
  sim.add_callback(print_energy_and_momentum_with_time_callback(2.0 / 100.0));

  // run simulation
  const auto nsteps = static_cast<unsigned long long>(ceil(2.0 / dt) + 1);
  cout << "Time step size: " << sim.get_dt() << '\n';
  cout << "Total simulation time: " << nsteps*sim.get_dt() << '\n';
  cout << "Running simulation for " << nsteps << " steps...\n\n";
  cout << "time U K E=U+K |p|\n" << "==================\n";

  mprof::tic();
  sim.simulate(nsteps);
  mprof::toc();

  return 0;
}

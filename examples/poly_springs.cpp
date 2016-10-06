#include <cmath>
#include <iostream>
#include <vector>
#include "mmd.hpp"

using namespace std;
using namespace mmd;

int main() {
  static const double dt = 1e-4;

  // spring potential, k = 1.0
  const_quad_spring_potential spr(1.0, -2.0, 1.0);

  vector<string> names({"2ii_025_p", "2ii_100_p", "2ii_025_m", "2ii_100_m", "2ii_200"});

  vector<simulation> sims({
      { molecular_id::Test, "2ii_025_p_init.xyz", "2ii_025_p_init.uvw", &spr, dt },
      { molecular_id::Test, "2ii_100_p_init.xyz", "2ii_100_p_init.uvw", &spr, dt },
      { molecular_id::Test, "2ii_025_m_init.xyz", "2ii_025_m_init.uvw", &spr, dt },
      { molecular_id::Test, "2ii_100_m_init.xyz", "2ii_100_m_init.uvw", &spr, dt },
      { molecular_id::Test, "2ii_200_init.xyz", "2ii_200_init.uvw", &spr, dt }
    });

  // check energy is conserved (every 100 time steps)
  for (size_t i = 0; i < sims.size(); ++i) {
    sims[i].add_callback(check_energy(10*dt, 1e-6));
    
    // check momentum is conserved
    // momentum is NOT conserved for this problem...
    //sim.add_callback(check_momentum(dt, 1e-7));
    
    // store energy and momentum
    sims[i].add_callback(
        save_energy_and_momentum_with_time_callback(names[i] + string("_data.csv"), 20.0 / 1000.0));
    
    // store mass positions and velocities
    sims[i].add_callback(save_values_with_time_callback(names[i] + string("_dofs.csv"), 20.0 / 40.0, 
          {
            [](const simulation& sim) -> double { 
              return sim.get_positions()(0, 0); // get x position of mass
            },
            [](const simulation& sim) -> double { 
              return sim.get_velocities()(0, 0); // get x component of velocity
            }
          }));

    // save mass position for visualization
    sims[i].add_callback(save_xyz_callback(names[i] + string(".xyz"), 20.0 / 40.0, positions));

    // print energy and momentum
    sims[i].add_callback(print_energy_and_momentum_with_time_callback(20.0 / 100.0));
  }

  for (auto& sim : sims) {
    // run simulation
    const auto nsteps = static_cast<unsigned long long>(ceil(20.0 / dt) + 1);
    cout << "Time step size: " << sim.get_dt() << '\n';
    cout << "Total simulation time: " << nsteps*sim.get_dt() << '\n';
    cout << "Running simulation for " << nsteps << " steps...\n\n";
    cout << "time U K E=U+K |p|\n" << "==================\n";

    mprof::tic();
    sim.simulate(nsteps);
    mprof::toc();
  }

  return 0;
}

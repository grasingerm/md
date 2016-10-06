#include <cmath>
#include <iostream>
#include <fstream>
#include <utility>
#include "mmd.hpp"

using namespace std;
using namespace mmd;

int main() {
  ofstream outfile1("bonus_eta10000.csv", ofstream::out);
  ofstream outfile2("bonus_eta00100.csv", ofstream::out);
  ofstream outfile3("bonus_eta00001.csv", ofstream::out);

  vector<ofstream*> outfiles;
  outfiles.push_back(&outfile1); 
  outfiles.push_back(&outfile2); 
  outfiles.push_back(&outfile3); 

  static const double dt = 1e-4;

  const_well_params_LJ_potential pot(1.0, 1.0);
  std::vector<double> etas({1.0, 0.01, 0.0001});

  std::vector<std::string> files({"2.txt", "3.txt", "3_2.txt", "3_3.txt", 
                                  "4.txt", "4_2.txt", "4_3.txt", "5.txt", 
                                  "6.txt", "7.txt", "8.txt", "9.txt"});

  for (const auto& file : files) {
    for (unsigned i = 0; i < 3; ++i) {
      quenched_velocity_verlet qv(etas[i]); 
      simulation sim(molecular_id::Test, file.c_str(), &pot, dt, qv);
 
      // print energy and momentum
      sim.add_callback(print_energy_and_momentum_with_time_callback(2.0 / 4));

      // run simulation
      const auto nsteps = static_cast<unsigned long long>(ceil(2.0 / dt) + 1);
      cout << "Eta = " << etas[i] << ", data file = " << file << '\n';
      cout << "Time step size: " << sim.get_dt() << '\n';
      cout << "Total simulation time: " << nsteps*sim.get_dt() << '\n';
      cout << "Running simulation for " << nsteps << " steps...\n\n";
      cout << "time U K E=U+K |p|\n" << "==================\n";

      mprof::tic();
      sim.simulate(nsteps);
      mprof::toc();

      cout << "\n\n";

      const int N = file[0] - '0';
      *outfiles[i] << N << ',' 
                   << potential_energy(sim) / static_cast<double>(N) << '\n';

    }
  }

  return 0;
}

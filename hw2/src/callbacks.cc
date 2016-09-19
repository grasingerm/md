#include "callbacks.hh"

namespace mmd {

callback save_pos_callback(const char* fname, const double time_out) {
  ofstream outfile(fname, ofstream::out);
  return [&](const vector<molecule<1>>& mols, const double t) mutable {
    auto threshold = 1;
    if (t / time_out > threshold) {
      outfile << t;
      for (const auto& mol : mols)
        outfile << ',' << mol.q[0];
      outfile << '\n';

      threshold++; // next time to process
    }
  };
}

}

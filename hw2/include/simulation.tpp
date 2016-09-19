namespace mmd {

void simulation<D>::simulate(const unsigned nsteps) {
  const n = molecules.size();
  double t = 0.0;
  vector<array<double, D>> forces(n, {0.0, 0.0});

  for (unsigned k = 0; k < nsteps; ++k) {

    for (auto& force : forces) {
      force[0] = 0.0;
      force[1] = 0.0;
    }

    // Add single particle potentials to instantaneous forces
    for (const auto ∇φ1 : ∇φ1s) {
      for (unsigned i = 0; i < n; ++i) {
        forces[i] += ∇φ1(molecules[i].q);
      }
    }

    // Add pairwise particle potentials to instantaneous forces
    for (const auto ∇φ2 : ∇φ2s) {
      for (unsigned i = 0; i < n-1; ++i) {
        for (unsigned j = i; j < n; ++j) {
          const unordered_pair<molecular_id> key { molecules[i].id, 
                                                   molecules[j].id };
          const auto& params = interactions[key];
          const auto ∇φij = ∇φ2(params[0], params[1], 
                                molecules[i].q - molecules[j].q);
          forces[i] += ∇φij;
          forces[j] -= ∇φij;
        }
      }
    }

    // Calculate instantaneous acceleration and integrate
    for (unsigned i = 0; i < n; ++i) {
      const mass = masses[molecules[i].id];
      const array<double, D> a { forces[i][0] / mass, forces[i][1] / mass };
      time_int(molecules[i], a, Δt); 
    }

    t += Δt;

    for (const auto callback : callbacks) callback(molecules, t);

  }
}

}

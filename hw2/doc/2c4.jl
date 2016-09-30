using Winston;

data = readdlm("../prob/2i_data.csv", ',');
dofs = readdlm("../prob/2i_dofs.csv", ',');

plot(data[:, 1], data[:, 2], "k--", dofs[:, 1], dofs[:, 2], "k-", dofs[:, 1], dofs[:, 3], "k:");
legend(["U", "x", "v"], [0.05, 0.2]);
xlabel("time (sec)");
savefig("2c_xvE.png");


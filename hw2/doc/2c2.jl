using PyPlot;

data = readdlm("../prob/2i_dofs.csv", ',');

scatter(data[:, 2], data[:, 3]);
xlabel("x");
ylabel("v");
savefig("2c_manifold.png");


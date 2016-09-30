using PyPlot;

dofs = readdlm("../prob/2i_dofs.csv", ',');
data = readdlm("../prob/2i_data.csv", ',');

U(x) = x^2 / 2;
xs = linspace(-1.8, 1.8, 500);

plot(xs, map(U, xs));
scatter(dofs[:, 2], map(U, dofs[:, 2]));
xlabel("x");
ylabel("U");
savefig("2c_pot.png");


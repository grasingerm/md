using PyPlot;

names = ["2ii_025_p", "2ii_100_p", "2ii_025_m", "2ii_100_m", "2ii_200"];

U(x) = x^4 - 2*x^2 + 1;
xs = linspace(-1.8, 1.8, 500);


clf();
data1 = readdlm("../prob/2ii_025_p_dofs.csv", ',');
plot(data1[:, 1], data1[:, 2], "b-", data1[:, 1], data1[:, 3], "b--");
data2 = readdlm("../prob/2ii_025_m_dofs.csv", ',');
plot(data2[:, 1], data2[:, 2], "r-", data2[:, 1], data2[:, 3], "r--");
xlabel("time (sec)");
legend(["x, x0 = 1","v, v0 = 0.707", "x, x0 = -1","v, v0 = -0.707"]);
savefig("2d_2ii_025.png");

clf();
plot(xs, map(U, xs));
scatter(data1[:, 2], map(U, data1[:, 2]); marker="o");
scatter(data2[:, 2], map(U, data2[:, 2]); marker="x");
xlabel("x");
ylabel("U");
savefig("2d_2ii_025_pot.png");

clf();
data1 = readdlm("../prob/2ii_100_p_dofs.csv", ',');
plot(data1[:, 1], data1[:, 2], "b-", data1[:, 1], data1[:, 3], "b--");
data2 = readdlm("../prob/2ii_100_m_dofs.csv", ',');
plot(data2[:, 1], data2[:, 2], "r-", data2[:, 1], data2[:, 3], "r--");
xlabel("time (sec)");
legend(["x, x0 = 1","v, v0 = 1.414", "x, x0 = -1","v, v0 = -1.414"]);
savefig("2d_2ii_100.png");

clf();
plot(xs, map(U, xs));
scatter(data1[:, 2], map(U, data1[:, 2]); marker="o");
scatter(data2[:, 2], map(U, data2[:, 2]); marker="x");
xlabel("x");
ylabel("U");
savefig("2d_2ii_100_pot.png");

clf();
data = readdlm("../prob/2ii_200_dofs.csv", ',');
plot(data[:, 1], data[:, 2], "k-", data[:, 1], data[:, 3], "k--");
xlabel("time (sec)");
legend(["x","v"]);
savefig("2d_2ii_200.png");

clf();
plot(xs, map(U, xs));
scatter(data[:, 2], map(U, data[:, 2]));
xlabel("x");
ylabel("U");
savefig("2d_2ii_200_pot.png");

clf();
markers = ["o", "x", "o", "x", "^"];
labels = ["E = 0.25", "E = 1.0", "E = 0.25", "E = 1.0", "E = 2.0"];
for (name, marker, lab) in zip(names, markers, labels)
  data = readdlm("../prob/" * name * "_dofs.csv", ',');
  scatter(data[:, 2], data[:, 3]; marker=marker, label=lab);
end

xlabel("x");
ylabel("v");
legend();
xlim([-2.5; 2.5]);
ylim([-3.5; 3.5]);
savefig("2d_manifold.png");

using PyPlot;

data = readdlm("../prob/3_data.csv", ',');

plot(data[:, 1], data[:, 2], "k-"; label="U");
plot(data[:, 1], data[:, 3], "k--"; label="K");
plot(data[:, 1], data[:, 4], "k:"; label="H=K+U");
legend(; loc=5);
xlabel("time (sec)");
savefig("3_energy.png");

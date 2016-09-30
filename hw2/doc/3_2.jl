using PyPlot;

data = readdlm("../prob/3_mom.csv", ',');

plot(data[:, 1], data[:, 2], "k-"; label="p_x");
plot(data[:, 1], data[:, 3], "k--"; label="p_y");
plot(data[:, 1], data[:, 4], "k:"; label="p_z");
legend(; loc=2);
xlabel("time (sec)");
ylabel("momentum");
savefig("3_mom.png");

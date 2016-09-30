data = readdlm("../prob/2i_data.csv", ',');

plot(data[:, 1], data[:, 2], "k-"; label="U");
plot(data[:, 1], data[:, 3], "k--"; label="K");
plot(data[:, 1], data[:, 4], "k:"; label="H=K+U");
legend(; bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.);
xlabel("time (sec)");
savefig("2c_energy.png");


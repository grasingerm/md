using Winston;

const x0 = 0.0;
const v0 = sqrt(2);

x(t) = x0 * cos(t) + v0 * sin(t)
v(t) = -x0 * sin(t) + v0 * cos(t)

ts = linspace(0, 20, 1000);
xs = map(x, ts);
vs = map(v, ts);

data = readdlm("../2i_dof.csv", ',');

plot(ts, xs, "k-", data[:, 1], data[:, 2], "ko", ts, vs, "k-.", data[:, 1], data[:, 3], "k^");
legend(["x, analytical", "x, MD", "v, analytical", "v, MD"]);
xlabel("time (sec)");
ylim([-2.0, 2.5]);
savefig("2c_comparison.png");

xerr = map(i -> abs(data[i,2]-x(data[i,1])), 1:size(data, 1)); 
verr = map(i -> abs(data[i,3]-v(data[i,1])), 1:size(data, 1)); 
plot(data[:, 1], xerr, "ko", data[:, 1], verr, "k^");
legend(["|x - x'|", "|v - v'|"]);
xlabel("time (sec)");
savefig("2c_error.png");


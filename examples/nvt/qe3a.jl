using DataFrames
using PyPlot
using GLM

D1s = zeros(10);
D2s = zeros(10);
for i = 1:10
d, hds = readdlm("qe3_$i.csv", ','; header=true);
times = vec(d[:, 1]);
plot(times, vec(d[:, 10]), "-",
     times, vec(d[:, 11]), "--",
     times, vec(d[:, 12]), "-.");
title("Atomic species with mass, \$m = 1\$");
legend(["\$\\alpha = x\$", "\$\\alpha = y\$", "\$\\alpha = z\$"]; loc=4);
xlabel("Dimensionless time");
ylabel("\$\\frac{1}{2N_1} \\sum_{i=1}^{N_1} \\left[r_{i,\\alpha}(t) - r_{i,\\alpha}(0)\\right]^2\$");
savefig("diffusion_1_$i.png");
clf();

plot(times, vec(d[:, 14]), "-",
     times, vec(d[:, 15]), "--",
     times, vec(d[:, 16]), "-.");
title("Atomic species with mass, \$m = 2\$");
legend(["\$\\alpha = x\$", "\$\\alpha = y\$", "\$\\alpha = z\$"]; loc=4);
xlabel("Dimensionless time");
ylabel("\$\\frac{1}{2N_2} \\sum_{i=1}^{N_2} \\left[r_{i,\\alpha}(t) - r_{i,\\alpha}(0)\\right]^2\$");
savefig("diffusion_2_$i.png");
clf();

plot(times, vec(d[:, 13]), "-", times, vec(d[:, 17]), "--");
legend(["\$m = 1\$", "\$m = 2\$"]; loc=4);
xlabel("Dimensionless time");
ylabel("\$\\frac{1}{6N_j} \\sum_{i=1}^{N_j} \\sum_{\\alpha=x,y,z} \\left[r_{i,\\alpha}(t) - r_{i,\\alpha}(0)\\right]^2\$");
savefig("diffusion_comparison.png");
clf();

d = readtable("qe3_$i.csv");
lm1 = lm(D1 ~ time, d);
println("coeffs for D1: $(coef(lm1))");
D1s[i] = coef(lm1);
lm2 = lm(D2 ~ time, d);
println("coeffs for D2: $(coef(lm2))");
D2s[i] = coef(lm1);
end

println("avgs");
println("D1: $(sum(D1s) / length(D1s))");
println("D2: $(sum(D2s) / length(D2s))");

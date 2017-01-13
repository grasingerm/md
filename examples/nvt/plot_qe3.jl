using PyPlot;

data = readdlm("diffusion_data.csv", ','; skipstart=1);
ts = vec(data[:, 1]);
dx1s = vec(data[:, 2]);
dy1s = vec(data[:, 3]);
dz1s = vec(data[:, 4]);
dj1s = vec(data[:, 5]);
dx2s = vec(data[:, 6]);
dy2s = vec(data[:, 7]);
dz2s = vec(data[:, 8]);
dj2s = vec(data[:, 9]);

clf();
plot(ts, dx1s, "-", ts, dy1s, "--", ts, dz1s, ":")
legend(["Dx", "Dy", "Dz"]; loc=4);
title("Diffusion for species with \$m=1\$");
savefig("dx1s.png");
xlim([0, 100]);
savefig("dx1s-zoom.png");

clf();
plot(ts, dx2s, "-", ts, dy2s, "--", ts, dz2s, ":")
legend(["Dx", "Dy", "Dz"]; loc=4);
title("Diffusion for species with \$m=2\$");
savefig("dx2s.png");
xlim([0, 100]);
savefig("dx2s-zoom.png");

clf();
plot(ts, dx1s, "-", ts, dx2s, "--")
legend(["Dx1", "Dx2"]; loc=4);
xlim([0, 100]);
savefig("dxs-compared.png");
xlim([0, 15]);
savefig("dxs-compared-zoomed.png");

sum = 0.0;
for i=1:length(dj1s)-1
  sum += dj1s[i+1]-dj1s[i] / ts[i];
end
@show sum / (length(dj1s)-1);

sum = 0.0;
for i=1:length(dj2s)-1
  sum += dj2s[i+1]-dj2s[i] / ts[i];
end
@show sum / (length(dj2s)-1);

using PyPlot;

for (eta, etastr) in zip([1.0, 0.01, 0.0001], ["10000", "00100", "00001"])

  data = readdlm("../prob/bonus_eta" * etastr * ".csv", ',');
  plot(data[:,1], data[:,2]);

end
xlabel("N");
ylabel("U / N");
legend(["eta = 1.0","eta = 0.01","eta = 0.0001"]);
savefig("bonus.png");

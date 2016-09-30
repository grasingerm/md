using PyPlot;

markers = ["o", "x", "s", "p", "*", "+", "<", "v", ">", "^"];
for (eta, etastr) in zip([1.0, 0.01, 0.0001], ["10000", "00100", "00001"])

  clf();
  data = readdlm("../prob/bonus_eta" * etastr * ".csv", ',');
  for i = 1:size(data, 1)
    scatter(data[i,1], data[i,2]; marker=markers[convert(Int, data[i,1])]);
  end

  xlabel("N");
  ylabel("U / N");
  savefig("bonus_eta" * etastr * ".png");
end

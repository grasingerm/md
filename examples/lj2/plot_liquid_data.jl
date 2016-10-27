using PyPlot;

function rolling_avg(u)
  sum = 0.0;
  result = zeros(length(u));
  for i=1:length(u)
    sum += u[i];
    result[i] = sum / i;
  end
  return result;
end

tscale = 121.0;
σ = 3.4e-10;
ϵ = 1.65e-21;
pscale = ϵ / σ^3;

data = readdlm("liquid_data.csv",',');

ts = vec(data[:,1]);
Ts = vec(data[:,5])*tscale;
Ps = vec(data[:,6])*pscale*1e-6;
Is = vec(data[:,7])*pscale*1e-6;
Vs = vec(data[:,8])*pscale*1e-6;

plot(ts, rolling_avg(Ps), "k-", ts, rolling_avg(Is), "k--", ts, rolling_avg(Vs), "k-.");
legend(["P (MPa)", "I (MPa)", "V (MPa)"]);
figure();
plot(ts, rolling_avg(Ts));

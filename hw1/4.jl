using Winston;

const A12 = 12.13;
const A6 = 14.45;

f1(x) = 4 * ((1/x)^12 - (1/x)^6);
f2(x) = 2 * (A12*(1/x)^12 - A6*(1/x)^6);

zero1 = 1.0;
zero2 = (A6/A12)^(-1/6);
min1x = (1/2)^(-1/6);
min1y = f1(min1x);
min2x = (1/2 * A6/A12)^(-1/6);
min2y = f2(min2x);

r = linspace(0.95, 2.5, 100);

p = plot(r, map(f1, r), "k-", r, map(f2, r), "k--");
xlabel("r / σ");
ylabel("u / ϵ");
legend(["LJ"; "LJ nn FCC"], [0.4; 0.2]);
add(p, Points([zero1; zero2], [0.0; 0.0], color="blue", symbolkind="circle"));
add(p, Points([min1x; min2x], [min1y; min2y], color="green", symbolkind="cross"));
savefig("potentials.png");

println("\n\nDONE!\n");

const L = 7.16636

ccount(s, c) = length(split(s, c)) - 1;

open("qe3.csv", "w") do w
  write(w, "time,d1x,d1y,d1z,d1,d2x,d2y,d2z,d2,D1X,D1Y,D1Z,D1,D2X,D2Y,D2Z,D2\n");
  open("qe3.xyz") do f

    idx = 1
    r_curr = zeros(3, 256);
    r_prev = zeros(3, 256);
    d_total = zeros(4, 2);
    ids = Array{Int}(256);
    curr_time = 0.0
    init_flag = true;
    prev_flag = false;

    for line in eachline(f)
      line = strip(line);
      ns = ccount(line, ' ');
      if ns == 0
        @assert(parse(line) == 256);
        if init_flag
          init_flag = false
        else
          if !prev_flag
            r_prev = copy(r_curr);
            idx = 1
            prev_flag = true
          else
            idx = 1
            # record data
            sums = zeros(4, 2);
            for k = 1:256, alpha = 1:3
              j = ids[k];
              dr = abs(r_curr[alpha, k] - r_prev[alpha, k]);
              sums[alpha, j] += (dr < L / 2) ? dr^2 : (dr - L)^2;
            end
            map!(s -> s / 256, sums);
            sums[4, :] = map(j -> sum(sums[:, j]) / 3, 1:2);
            d_total += sums;

            # output data
            write(w, "$curr_time");
            for j=1:2, alpha=1:4
              write(w, ",$(sums[alpha, j])");
            end
            for j=1:2, alpha=1:4
              write(w, ",$(d_total[alpha, j])");
            end
            write(w, '\n');

            r_prev = copy(r_curr);
          end
        end
      elseif ns == 1
        curr_time = parse(split(line, ' ')[2]);
      else
#        println("Parsing $line...");
        columns = split(line, ' ');
        ids[idx] = parse(columns[1]) - 2;
        r_curr[1, idx] = parse(columns[2]);
        r_curr[2, idx] = parse(columns[3]);
        r_curr[3, idx] = parse(columns[4]);
        idx += 1
      end
    end
  end
end

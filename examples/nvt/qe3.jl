const L = 7.16636

ccount(s, c) = length(split(s, c)) - 1;

#open("qe_extra.csv", "w") do w_extra
open("qe3.csv", "w") do w
  write(w, "time,d1x,d1y,d1z,d1,d2x,d2y,d2z,d2\n");
  open("qe3.xyz") do f

    idx = 1
    r_curr = zeros(3, 256); # current position
    r_prev = zeros(3, 256); # previous position
    dr = zeros(3, 256); # total distance travelled
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
              dr_k = r_curr[alpha, k] - r_prev[alpha, k];
              if dr_k < -(L / 2)
                dr_k += L;
              elseif dr_k > (L / 2)
                dr_k -= L;
              end
              dr[alpha, k] += dr_k;
              sums[alpha, j] += (dr[alpha, k])^2;
#=              
              println("id = $j");
              println("r_c = $(r_curr[alpha, k])");
              println("r_p = $(r_prev[alpha, k])");
              println("dr_k = $(dr_k)");
              println("Σdr = $(dr[alpha, k])");
              println("Σ = $(sums[alpha, j])");
              println("L = $L");
=#
            end
            map!(s -> s / 256, sums);
            for j = 1:2, alpha = 1:3
              sums[4, j] += sums[alpha, j] / 3.0;
            end

            # output data
            write(w, "$curr_time");
            for j=1:2, alpha=1:4
              write(w, ",$(sums[alpha, j])");
            end
            write(w, '\n');
#=
            write(w_extra, "$curr_time");
            for k=1:256, alpha=1:3
              write(w_extra, ",$(dr[alpha, k])");
            end
            write(w_extra, '\n');
=#
            r_prev = copy(r_curr);
          end
        end
      elseif ns == 1
        curr_time = parse(split(line, ' ')[2]);
      else
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
#end

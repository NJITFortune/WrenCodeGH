function tims = zca(sig,Fs);

  intvl = 1/Fs;

  tim = 1/Fs:1/Fs:length(sig)/Fs;
  
  bint = zeros(length(sig),1);
  
  pos = find(sig > 0);
  
  bint(pos) = ones(length(pos),1);

  ups = find(diff(bint) == 1);

  zc = 
  
  
  crossing_intervals = (tim(1:end-1).*tim(2:end) <= 0);
  left_ends = (x(1:end-1))(crossing_intervals);
  right_ends = (x(2:end))(crossing_intervals);
  left_vals = (y(1:end-1))(crossing_intervals);
  right_vals = (y(2:end))(crossing_intervals);
  mid_points = (left_ends+right_ends)/2;

  zero_intervals = find(left_vals==right_vals);
 
  retval1 = mid_points(zero_intervals);
  left_ends(zero_intervals) = [];
  right_ends(zero_intervals) = [];
  left_vals(zero_intervals) = [];
  right_vals(zero_intervals) = [];
  retval2=left_ends-(right_ends-left_ends).*left_vals./(right_vals-left_vals);
  retval = union(retval1,retval2);

endfunction

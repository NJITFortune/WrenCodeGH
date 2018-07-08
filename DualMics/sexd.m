function [ sex ] = sexd( datam, dataf, sylnum, dat )
%[ out ] = sexd( md, fd, dat );
%   This takes the output of wrenxc.  Need to adjust for amplitude
%   differences??

%% Sex determination

% Is it male or female??

% figure(1); 
hold on;

for j = length(sylnum):-1:1;

    strt = min([dataf(sylnum(j)).sylind(1) datam(sylnum(j)).sylind(1)]);
    stp = max([dataf(sylnum(j)).sylind(2) datam(sylnum(j)).sylind(2)]);

    cnt(j) = length(find(dat(strt:stp) > 0)) - length(find(dat(strt:stp) < 0));
    mag(j) = sum(dat(strt:stp));

    plot(dat(strt:stp));
    
end;
    if mean(mag) < 0; sex = 0; end;
    if mean(mag) > 0; sex = 1; end;



%subplot(2,2,1); oscson(fem(strt_ind:stop_ind),Fs);
%subplot(2,2,2); 


end


function out = bs_corrhist(f, m, win, binwidth, overlap, pt)
%OUT = BS_CORRHIST(F, M, WIN, BINWIDTH, OVERLAP, PT);
% F and M are the spike cell arrays from bs_converter
% WIN is start and end times in seconds relative to onset of stimulus
% BINWIDTH is the width of the histogram window in milliseconds
% OVERLAP is the percent (max 0.9) of the overlap of windows.
% If PT is one, plot the data.

if nargin < 6; pt = 0; end;

binwidth = binwidth/1000; % convert to seconds
overlap = 1 - overlap;

%% The Histograms

% How many bins?

    len = win(2) - win(1);
    stepsize = binwidth*overlap;
    binnum = round( len / stepsize);
    fhist = zeros(binnum,1);
    mhist = fhist;
    out.tim = fhist;

% Fill the bins, cycle by bin

for i = 1:binnum;
    binstart = (i-1)*stepsize + win(1);
    binend   =     binstart + binwidth;   
    out.tim(i) = binstart + binwidth/2;
    
    % and cycle here by reps of the stimulus
    
    for j = 1:length(f);       
        frepbincnt = sum (f{j} > binstart & f{j} <= binend);
        fhist(i) = fhist(i) + frepbincnt;
    end;
    for j = 1:length(m);       
        mrepbincnt = sum (m{j} > binstart & m{j} <= binend);
        mhist(i) = mhist(i) + mrepbincnt;
    end;    
end;

out.mhist = mhist/max(mhist);
out.fhist = fhist/max(fhist);
out.reps = length(f);


%% Make the cross-correlation

[out.xcorr, lag] = xcorr(out.mhist, out.fhist);
out.xcorr(end+1) = 0;
[out.cmax, i] = max(abs(out.xcorr));
out.xctim = stepsize:stepsize:len*2;
out.xctim = out.xctim - len;
out.offset = lag(i)*stepsize;

%% Plot

if pt == 1;
   subplot(211);
    plot(out.tim, out.mhist, 'b', out.tim, out.fhist, 'm');
   subplot(212);
    plot(out.xctim, out.xcorr, 'k', out.xctim(i), out.xcorr(i), 'r*');
end;
   
end

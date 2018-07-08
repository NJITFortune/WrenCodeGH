function out = bs_xcorr(sig1, Fs1, tim1, sig2, Fs2, tim2, win, typ)
% out = bs_xcorr(sig1, Fs1, tim1, sig2, Fs2, tim2, win, typ)
% sig1 and sig2 can be either stimulus or swPSTH data
% Fs1 and Fs2 are the respective sample rates
% tim1 and tim2 are time bases for each
% win is the time window for analysis in seconds - must be close to
% stimulus length for accurate results, e.g. [0 7]
% typ is the classs of stimulus - 1 for stimulus and 0 for swPSTH, e.g. [1 0]
%
% output is the cross-correlation (xcorr), centered on zero.
% structure includes the time-matched inputs.


%% Prepare stimuli
xFs = min([Fs1 Fs2]);

tt1 = find(tim1 > win(1) & tim1 < win(2));
tt2 = find(tim2 > win(1) & tim2 < win(2));

if typ(1) == 1; % A stimulus
[b,a] = butter(2,(20/Fs1),'low');
s = filtfilt(b,a,abs(sig1(tt1) - mean(sig1(tt1))));
s = s / max(s);
ts = tim1(tt1);
elseif typ(1) == 0; % A swPSTH
    s = sig1(tt1);
    s = s / max(s);
    ts = tim1(tt1);
end;
if typ(2) == 1; % A stimulus
[b,a] = butter(2,(20/Fs2),'low');
r = filtfilt(b,a,abs(sig2(tt2) - mean(sig2(tt2))));
r = r / max(r);
tr = tim2(tt2);
elseif typ(2) == 0; % A swPSTH
    r = sig2(tt2);
    r = r / max(r);
    tr = tim2(tt2);
end;

%% Fix sample rates
if Fs1 ~= Fs2; % If different sample rates, need to downsample faster
    if Fs1 > Fs2;
        s = decimate(s,Fs1/Fs2);
        ts = decimate(ts,Fs1/Fs2);
    elseif Fs2 > Fs1;
        r = decimate(r,Fs2/Fs1);
        tr = decimate(tr,Fs2/Fs1);
    end;
end;

%% Do the correlation and filter it
xs = xcorr(s,r);
[b,a] = butter(2,(5/Fs2),'high');
xsf = filtfilt(b,a,xs);

%% Do the cross-spectral density

[Cxy,F] = mscohere(s,r,[],[],256,xFs);
figure(1); plot(F,Cxy);


%% Plot the results
figure(2);
subplot(212); plot(ts, s, 'b');
hold on;
subplot(212); plot(tr, r, 'r');

subplot(211);
midx = ceil(length(xs)/2);
out.xtim(1:midx-1) = -(midx-1)/xFs:1/xFs:-1/xFs;
out.xtim(midx) = 0;
out.xtim(midx+1:midx+midx-1) = 1/xFs:1/xFs:(midx-1)/xFs;
plot(out.xtim,[xs xsf]);

%% Output the results
out.s = s;
out.ts = ts;
out.r = r;
out.tr = tr;
out.xs = xs;
out.xsf = xsf;


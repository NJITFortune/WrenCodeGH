function out = wrenv(sig,Fs,zz);
% Usage: out = wrenv(signal, Fs, [start stop]);
% start and stop are the times (seconds) of a good silent section

tim = 1/Fs:1/Fs:length(sig)/Fs;

out = lpf(2*abs(hpf(sig,Fs) - mean(hpf(sig,Fs))),Fs,[100 7]);

out = out - mean(out(find(tim > zz(1) & tim < zz(2))));


function [freq tim] = TOE2Samples(toe,Fs,smoothlength);
% computes instantaneous frequency of input signal y based on sampling rate
% fs using threshold value "zerocross" (usually zero V)

maxtim = toe(end); 


tim = 1/Fs:1/Fs:maxtim;
instfreq = 1 ./ diff(toe);

freq = interp1(toe(2:end),instfreq,tim);

if nargin == 3;
smoothlength = smoothlength*Fs;
freq = smooth(freq, smoothlength);
end;


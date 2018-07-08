function out = toneharm(fun, dur, Fs);
% this makes a harmonic pair of sinewaves

tim = 1/Fs:1/Fs:dur;

nop = zeros(2*dur*Fs,1);

ramp = ones(length(tim),1);
rampdur = dur*Fs*0.1;

ramp(1:rampdur) = 0.0001:1/rampdur:1;
ramp(end-(rampdur-1):end) = 1:-1/rampdur:0.001;

tone1 = 0.5 * sin(2*pi*fun*tim); tone1 = tone1.*ramp';
tone2 = 0.5 * sin(2*pi*2*fun*tim); tone2 = tone2.*ramp';
tone3 = tone1 .+ tone2;

out = [tone1 nop' tone2 nop' tone3 nop'];
out = out';


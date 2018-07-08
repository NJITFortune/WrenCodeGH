function out = wrenalign(sig1, sig2, Fs);

hil1 = abs(hilbert(sig1));
hil2 = abs(hilbert(sig2));

xx = xcorr(hil1, hil2);
[maxamp xoff] = max(xx);

mid = round((length(xx)+1) / 2);

offset = mid - xoff;

tim = offset / Fs;

%printf ("Start offset of  %f seconds. \n", tim);

%%%%%% pad the signals - start with the front for alignment

startpad = zeros(abs(offset), 1);

if ( offset > 1 ); sig1a = [startpad' sig1']; sig2a = sig2; end;
if ( offset < 1 ); sig2a = [startpad' sig2']; sig1a = sig1; end;

%%%%%% pad the signals -  make them the same length

maxlen = max([ length(sig2a) length(sig1a) ]);
longzeros = zeros(maxlen,1);

sig1b = longzeros; sig1b(1:length(sig1a)) = sig1a;
sig2b = longzeros; sig2b(1:length(sig2a)) = sig2a;

%% Adding the Hilbert...
hil1 = medfilt1( abs ( hilbert ( sig1b ) ), 100);
hil2 = medfilt1( abs ( hilbert ( sig2b ) ), 100);

out.sig1 = sig1b;
out.sig2 = sig2b;
out.Fs = Fs;
out.offsetsamples = offset;
out.offsettime = tim;
out.hil1 = hil1;
out.hil2 = hil2;
out.tim = 1/Fs:1/Fs:length(sig1b)/Fs;


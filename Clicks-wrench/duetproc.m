function foo = duetproc(in)
% Usage out = duetproc(in)
% Where in is the structure produced by distantwren.m
% Depends on "hagaclicks" which depends on "clickplotter" "ginputc" and "oscson"
%

Fs = 1 / (in(1).ftim(2) - in(1).ftim(1)); % Get the samplerate
length(in)

%% Cycle through each duet in the file
for jj = 1:length (in) % length(in);
    figure(27); clf;
    ax(1) = subplot(211); specgram(in(jj).fc.fem, 2048, Fs, [], 2000); ylim([0 6000]); 
    ax(2) = subplot(212); specgram(in(jj).mc.mal, 2048, Fs, [], 2000); ylim([0 6000]); 
    linkaxes(ax,'x');

    foo(jj).m = hagaclics(in(jj).mc.mal, Fs, 0);
    foo(jj).f = hagaclics(in(jj).fc.fem, Fs, 0);

end;
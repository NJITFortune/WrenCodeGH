function out = wrenstartend(sig, Fs);
% function out = wrenstartend(signal, samplerate);
% signal = the data, 1 channel
% samplerate = samples per second
% out = the chopped sample
% Currently you click the immediate start of the song and
% then the immediate end of the song. This script then pads
% both ends, current 1 second each side. Edit the script to
% change the values of startbuf and endbuf if necessary.
% The script does not currently add zeros if there is less than 
% the pad on either side of the signal.

	startbuf = 1; endbuf = 1;

close all;

i=0;
tim = 1/Fs:1/Fs:length(sig)/Fs;

figure(1);
	subplot(211); specgram(sig,1024,Fs,[],1000);
	colormap('HOT'); ylim([400 5400]); caxis([-20 20]);
	subplot(212); plot(tim,sig); xlim([0 max(tim)]);

while i ~=1;

	k = ginput(2);

	subplot(212); hold on; 
	ax = axis; sx = [400 5400]; 
	plot([k(1) k(1)], [ax(3) ax(4)], 'g');
	plot([k(2) k(2)], [ax(3) ax(4)], 'r');
	hold off;

	subplot(211); hold on;
	plot([k(1) k(1)], [sx(1) sx(2)], 'g');
	plot([k(2) k(2)], [sx(1) sx(2)], 'r');
	hold off;

	i = input('OK? 1=yes, 0=redo');
end

	tt = find(tim > k(1)-startbuf & tim < k(2)+endbuf);
	out = sig(tt);

	close(1);
	ntim = 1/Fs:1/Fs:length(out)/Fs;

	figure(1);
        subplot(211); specgram(out,1024,Fs,[],1000);
        colormap('HOT'); ylim([400 5400]); caxis([-20 20]);
        subplot(212); plot(ntim,out); xlim([0 max(ntim)]);


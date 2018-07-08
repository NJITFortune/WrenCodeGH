function out = wrenpickpeaks(sig, Fs);
% function out = wrenpickpeaks(signal, samplerate);
% signal = the data
% samplerate = samples per seconds
% out is a structure with the times and frequencies of where the
% user clicked on the spectrogram.
% out.peaktime(x) = time in seconds from the beginning of the signal
% out.peakfreq(x) = frequency (Hz) of the peak that was clicked

close all;

i=0; p=0;

tim = 1/Fs:1/Fs:length(sig)/Fs;

% Get the smoothed envelope of the signal 
%% env = medfilt1(abs(hilbert(sig)),150);

figure(1);
	specgram(sig,1024,Fs,[],1000);
	colormap('HOT'); ylim([400 5400]); caxis([-20 20]);
figure(2);
	plot(tim,sig); xlim([0 max(tim)]); ax = axis;

while i ~= 1;

	figure(1);
	k = ginput(1);
		p=p+1; 
			out.peaktime(p)=k(1);
			out.peakfreq(p)=k(2);
	hold on; 
	plot([k(1) k(1)], [400 5400], 'r');
	plot(k(1), k(2), 'hg'); 
	hold off;

	figure(2); 
	y = find(tim >= k(1));
	hold on; 
	plot([k(1) k(1)], [ax(3) ax(4)], 'r');
	%%plot(k(1), env(y(1)), 'hg'); 
	hold off;
	xlim([ax(1) ax(2)]);

	a = input('Done? 1=yes, 0 or [return]=continue: ');
	if a == 1; i = 1; end;
end


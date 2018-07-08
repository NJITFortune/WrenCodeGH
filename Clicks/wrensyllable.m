function out = wrensyllable(struct);
% function out = wrensyllable(struct);
% struct needs to have struct.signal, struct.fs, struct.peaktime, struct.peakfreq.
% out gets out.syl which includes out.syl.start and out.syl.end.
% A bunch of analysis should be added to this script... currently it is only clicking.

close all;
Fs = struct.fs;
fdata = wrenfilter(struct.signal, Fs);

i=0; 
p=0;

sig = struct.signal;
tim = 1/Fs:1/Fs:length(struct.signal)/Fs;

env = medfilt1(abs(hilbert(fdata)),200);

figure(1); subplot(211);
	specgram(sig,1024,Fs,[],1000);
	colormap('HOT'); ylim([400 5400]); caxis([-20 20]);
figure(1); subplot(212);
	plot(tim,env); xlim([0 max(tim)]); ax = axis;

while i ~= 1;
	p=p+1;

	tt = find(tim >= struct.peaktime(p)-1 & tim < struct.peaktime(p)+1);

	figure(2);
	specgram(sig(tt),512,Fs,[],505); ylim([400 5400]);
	colormap('HOT'); caxis([-20 20]);
	hold on;
	foo = fdata(tt) - mean(fdata(tt));
	mx = max( abs(foo) );
	foo = 4400 + (foo * 1000/mx);
	xx = tim(tt) - min(tim(tt));
	plot(xx, foo, 'w');

	figure(2);
	hold on; plot(struct.peaktime(p)-min(tim(tt)), struct.peakfreq(p),'g*'); hold off;

	k = ginput(2);

	figure(2);
	hold on; 
	plot([k(1) k(1)], [400 5400], 'c');
	plot([k(2) k(2)], [400 5400], 'm');
	hold off;

	k(1) = k(1) + min(tim(tt));
	k(2) = k(2) + min(tim(tt));

	figure(1); subplot(211);
	hold on; 
	plot([k(1) k(1)], [400 5400], 'c');
	plot([k(2) k(2)], [400 5400], 'm');
	hold off;

        figure(1); subplot(212);
        hold on;
        plot([k(1) k(1)], [ax(3) ax(4)], 'g');
        plot([k(2) k(2)], [ax(3) ax(4)], 'r');
        hold off;

%	xlim([ax(1) ax(2)]);

	b = input('Ok? 1 if no: ');
		if b == 1


		else;

		out.syl(p).start = k(1);
		out.syl(p).end = k(2);
			
	a = input('Done? 1=yes, 0 or [return]=continue: ');
	if a == 1; i = 1; end;
		end;


end


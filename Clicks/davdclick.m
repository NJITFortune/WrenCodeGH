function struct = davdclick(sng, Fs);

tim = 1/Fs:1/Fs:length(sng)/Fs;
maxtim = tim(end);

bsx = 0; st = 0;

windwid = 3;
mvfd = windwid -0.5;

figure(1);

while bsx < tim(end);

%	subplot(211); 
%	specgram(sng,1024,Fs); axis([bsx bsx+5 500 5000]);
%	subplot(212);
%	plot(tim,sng); xlim([bsx bsx+5]);

	tt = find(tim >= bsx & tim <= bsx+windwid);
	oscson(sng(tt),Fs); 

	[x y] = ginput;

	xclik(st+1:st+length(x)) = x+bsx;
	st = length(xclik);

	bsx = bsx + mvfd;

end

struct.song = sng;
struct.Fs = Fs;
struct.time = tim;

for i = 1:length(xclik)-1;
	struct.syl(i,:) = [xclik(i) xclik(i+1)];
	struct.sylen(i) = diff(struct.syl(i,:));
	tmpa = find(tim > xclik(i) & tim < xclik(i)+0.1);
	tmpb = find(tim > xclik(i+1) & tim < xclik(i+1)+0.1);
	struct.ind(i,:) = [tmpa(1) tmpb(1)];

	% get loudness
	sy = struct.song(struct.ind(i,1):struct.ind(i,2));
	sy = sy - mean(sy);
	struct.loud(i) = sum(abs(sy));
	struct.meanloud(i) = struct.loud(i) / struct.sylen(i);

	% get freq
	fftdata = fft(sy);
	struct.spectrum{i} = abs(real(fftdata(1:round(end/2))));
	stepsize = Fs/round(length(sy));
	struct.fftfreqs{i} = stepsize:stepsize:Fs/2;
	mn = min([length(struct.spectrum{i}) length(struct.fftfreqs{i})]);
	struct.spectrum{i} = struct.spectrum{i}(1:mn);
	struct.specfilt{i} = medfilt1(struct.spectrum{i},20);
	struct.fftfreqs{i} = struct.fftfreqs{i}(1:mn);

	[gg hh] = max(struct.specfilt{i});
	struct.peakfreq(i) = struct.fftfreqs{i}(hh);

	thresh = find( struct.specfilt{i} > gg/4 );
	struct.minfreq(i) = struct.fftfreqs{i}(thresh(1));
	struct.maxfreq(i) = struct.fftfreqs{i}(thresh(end));
	struct.band(i) = struct.maxfreq(i) - struct.minfreq(i);

end


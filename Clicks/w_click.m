function out = w_click(sng, Fs);
% Usage: out = w_click(sng, Fs);

firstsample = 1;

% Make a time vector so that we know when is when
tim = 1/Fs:1/Fs:length(sng)/Fs;
maxtim = tim(end);

% bsx is the real start time, and we start at the beginning
bsx = 0; % st = 0;

% The width of the clicking window is 3 seconds
windwid = 3;

eos = 0;

% Here we scroll, roughly 3 seconds at a time, through the whole sample
while eos ~= 1;

	tt = find(tim >= bsx & tim <= bsx+windwid);

	[starts ends] = clicker(sng(tt),Fs);
	
if firstsample == 1; 
	out.syl(:,1) = starts; 
	out.syl(:,2) = ends; 

	for i = 1:length(starts);
		inds = find (tim > starts(i) & tim < ends(i));
		out.ind(i,:) = [min(inds) max(inds)];
	end;

end;

starts = starts + bsx;
ends = ends + bsx;

if firstsample > 1;
	aa = length(out.syl);
	bb = length(starts)+aa;
	out.syl(aa+1:bb,1) = starts;
	out.syl(aa+1:bb,2) = ends;

        for i = 1:length(starts);
                inds = find (tim > starts(i) & tim < ends(i));
                out.ind(aa+i,:) = [min(inds) max(inds)];
        end;

end;

	bsx = ends(end);

if bsx > maxtim - 2;
	eos = input('Enter 1 if you are finished clicking, 0 if not: ');
end;

	firstsample = firstsample + 1;

ns = length(starts); ts = length(out.syl(:,1));
fprintf('You added %i syllables for a total of %i syllables.\n', ns, ts);

end

out.song = sng;
out.Fs = Fs;
out.time = tim;

%for i = 1:length(xclik)-1;
%	struct.syl(i,:) = [xclik(i) xclik(i+1)];
%	struct.sylen(i) = diff(struct.syl(i,:));
%	tmpa = find(tim > xclik(i) & tim < xclik(i)+0.1);
%	tmpb = find(tim > xclik(i+1) & tim < xclik(i+1)+0.1);
%	struct.ind(i,:) = [tmpa(1) tmpb(1)];
%
%	% get loudness
%	sy = struct.song(struct.ind(i,1):struct.ind(i,2));
%	sy = sy - mean(sy);
%	struct.loud(i) = sum(abs(sy));
%	struct.meanloud(i) = struct.loud(i) / struct.sylen(i);
%
%	% get freq
%	fftdata = fft(sy);
%	struct.spectrum{i} = abs(real(fftdata(1:round(end/2))));
%	stepsize = Fs/round(length(sy));
%	struct.fftfreqs{i} = stepsize:stepsize:Fs/2;
%	mn = min([length(struct.spectrum{i}) length(struct.fftfreqs{i})]);
%	struct.spectrum{i} = struct.spectrum{i}(1:mn);
%	struct.specfilt{i} = medfilt1(struct.spectrum{i},20);
%	struct.fftfreqs{i} = struct.fftfreqs{i}(1:mn);
%
%	[gg hh] = max(struct.specfilt{i});
%	struct.peakfreq(i) = struct.fftfreqs{i}(hh);
%
%	thresh = find( struct.specfilt{i} > gg/4 );
%	struct.minfreq(i) = struct.fftfreqs{i}(thresh(1));
%	struct.maxfreq(i) = struct.fftfreqs{i}(thresh(end));
%	struct.band(i) = struct.maxfreq(i) - struct.minfreq(i);
%
%end


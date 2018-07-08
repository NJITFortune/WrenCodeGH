function struct = clicksyllables(sng, Fs);
% This function depends on oscson for plotting.
%
%
% Close all of the open figures.  That can be annoying! But it MUST be done...
close all;

% Make and time sequence and get the length of the signal.
tim = 1/Fs:1/Fs:length(sng)/Fs;
maxtim = tim(end);

% Initialize some variables.  yvals is for plotting verticle lines.
adv = 1; bsx = 0; st = 1; yvals = [0 5000];

% This is the width of the clicking window. 
windwid = 3;

% This is the overlap between the previous and the next window.
mvfd = windwid - 0.5;

%%%%%%%%%%%%%%%%%%% START OF CLICKING

% The main loop of the program - we go until we've done the whole thing.

while bsx < tim(end);

	% get our indices for the time window to display
	tt = find(tim >= bsx & tim <= bsx+windwid);

	if adv == 1; 
		close all; figure(1); oscson(sng(tt),Fs); 
		figprop = get(gcf,'Position'); set(gcf,'Position',[figprop(1) figprop(2) 800 400]);

		if st > 1; figure(1); 
		hold on; plot([xclick(st-1,2)-bsx xclick(st-1,2)-bsx], yvals, 'm'); hold off; 
		end;
	end;

	adv = 0;

	% get two clicks for the syllable
	clicksOK = 0;
	[x y] = ginput(2);

	% Plot what the user just clicked in a separate window with
	% 0.5 seconds on either side.

	figure(2);

		if x(1)+bsx > 0.5 & x(2)+bsx < tim(end)-0.5; 
			samp = find (tim >= bsx+x(1)-0.5 & tim < bsx+x(2)+0.5);
			oscson(sng(samp), Fs);
			magentaline = length(samp)/Fs - 0.5; 
			hold on;
                        plot([0.5 0.5], yvals, 'g', 'LineWidth',2);
                        plot([magentaline magentaline], yvals, 'm', 'LineWidth',2);
			hold off;
		end;
		if x(1)+bsx <= 0.5;
			samp = find (tim < bsx+x(2)+0.5);
			oscson(sng(samp), Fs);
			magentaline = length(samp)/Fs - 0.5; 
			hold on;
                        plot([x(1) x(1)], yvals, 'g', 'LineWidth',2);
                        plot([magentaline magentaline], yvals, 'm', 'LineWidth',2);
			hold off;
		end;
                if x(1)+bsx > 0.5 & x(2)+bsx > tim(end)-0.5;
                        samp = find (tim >= bsx+x(1)-0.5);
                        oscson(sng(samp), Fs);
                        hold on;
                        plot([0.5 0.5], yvals, 'g', 'LineWidth',2);
			matentaline = 0.5 + (x(2) - x(1));
                        plot([magentaline magentaline], yvals, 'm', 'LineWidth',2);
                        hold off;
                end;
                
        % The user gets the final say - ask whether the clicks are OK.
        	clicksOK = input('If good type 0 or hit return; bad type 1; if OK and done with sample 99: ');
	close(2);

	% If OK, then we add it to the final list
		if clicksOK ~= 1;

			xclick(st,1) = x(1) + bsx; xclick(st,2) = x(2) + bsx;
		
			if x(2) >= 2.5; 
				bsx = bsx + mvfd; 
				adv = 1;
			end;
			if x(2) < 2.5;
				figure(1); hold on; 
				plot([x(1) x(1)], yvals, 'g');
				plot([x(2) x(2)], yvals, 'm');
				hold off;
			end;

			clear x y; 
			st = st + 1;
		end
	if clicksOK == 99; bsx = bsx + 99; end;

end

%%%%%%%%%%%%%%%%%%% END OF CLICKING


struct.song = sng;
struct.Fs = Fs;
struct.time = tim;



%%%%%%%%%%%%%%%%%%% BEGIN Syllable-by-Syllable Analyses

for i = 1:length(xclick);
	struct.syl(i,1) = xclick(i,1);
	struct.syl(i,2) = xclick(i,2);

% Perhaps the most useful bits - the durations and the indices

	struct.sylen(i) = diff(struct.syl(i,:));
	tmpa = find(tim >= xclick(i,1) & tim <= xclick(i,2));
	struct.ind(i,:) = [tmpa(1) tmpa(end)];

% get loudness - which is kinda ridiculous, but kinda not

	sy = struct.song(struct.ind(i,1):struct.ind(i,2));
	sy = sy - mean(sy);
	struct.loud(i) = sum(abs(sy));
	struct.meanloud(i) = struct.loud(i) / struct.sylen(i);

% get frequency content of the syllable

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

% Trace the frequency of the syllable

	siglen = length(sy);
	tinytim = 1/Fs:1/Fs:siglen/Fs;
	stp = 0.005;
	startwin = 0;
	p = 1;
	nfft = 2048;

	while startwin < max(tinytim) + stp;

        	data = sy(find(tinytim > startwin & tinytim < startwin + stp));
        	L = length(data);

        	t_fftdata = fft(data,nfft)/L;
        	t_fftdata = 2 * abs(t_fftdata(1:nfft/2+1));
        	t_freqs = Fs/2*linspace(0,1,nfft/2+1);

        	[maxval idx] = max(t_fftdata);
        	peakamp(p) = maxval;
        	peakfreq(p) = t_freqs(idx);
        	stw(p) = startwin + stp/2;

        	p = p + 1; startwin = startwin + stp/2;

	end

	struct.freqtrace{i} = medfilt1(peakfreq(1:end-2),5);
	struct.freqtim{i} = stw(1:end-2);
	struct.slopemean(i) = mean(diff(peakfreq));
	struct.slopestd(i) = std(diff(peakfreq));
	struct.slopevar(i) = var(diff(peakfreq));
	[struct.freqtracepeak(i) indf] = max(peakfreq);
	struct.freqtracepeaktim(i) = stw(indf);
	struct.freqtracepeakpercent(i) = struct.freqtracepeaktim(i) / stw(end);

	clear peakfreq tinytim stw;

% Do Ofer's analysis

[q,w,e,r,t,y,u,k,o,p]=deriva(sy,Fs);

        struct.m_spec_deriv{i} = q;
        struct.m_AM{i} = w;
        struct.m_FM{i} = e;
        struct.m_Entropy{i} = r;
        struct.m_amplitude{i} = t;
        struct.gravity_center{i} = y;
        struct.m_PitchGoodness{i} = u;
        struct.m_Pitch{i} = k;
        struct.Pitch_chose{i} = o;
        struct.Pitch_weight{i} = p;

        struct.PW(i) = mean(p);
        struct.GC(i) = mean(y);

end

%%%%%%%%%%%%%%%%%%% END Syllable-by-Syllable Analyses

for j = 1:length(struct.slopemean)-1;
	struct.ISI(j) = struct.syl(j+1,1) - struct.syl(j,2);
end


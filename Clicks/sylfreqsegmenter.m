function out = sylfreqsegmenter(sig, Fs)
% This uses the frequency power to segment syllables
% Usage: out = syllabletracer(signal, samplerate);
% This is for use with gelisclics, but can be used standalone.
% It is designed specifically for duets of plain-tailed wrens.
% You will be prompted with what to do next!

% Make a time base, tinytim
	tinytim = 1/Fs:1/Fs:length(sig)/Fs;

% We will step through the sample with this precision (in seconds)
	stp = 0.005; % Step size
	wid = 0.020; % Width of the FFT

	nfft = 512;	% The size of the FFT
    
	p = 1;		% this is a counter - leave as is
    
% The main loop of the program
	startwin = 0;
	while startwin < tinytim(end) - wid;

	% Take a 'stp' window of the sample
        	data = sig(tinytim >= startwin & tinytim < startwin + wid);
        	L = length(data);

	% Perform an FFT on that data
        	t_fftdata = fft(data,nfft)/L;
        	t_fftdata = 2 * abs(t_fftdata(1:nfft/2+1));
        	t_freqs = Fs/2*linspace(0,1,nfft/2+1);

    		struct.fft{p} = t_fftdata; % Power
        	struct.freqs{p} = t_freqs; % Frequencies
	
	% Get the frequency at the peak amplitude
        	[maxval, idx] = max(t_fftdata); % Max power
        	peakamp(p) = maxval;
        	freq(p) = t_freqs(idx); % This is the frequency that has the most power

% 	% Get the slope to and from peak
% 		qwer = find(t_freqs <= 6000); qwer = qwer(end);
% 		struct.up(p) = mean(diff(t_fftdata(1:idx)));
% 		struct.dwn(p) = mean(diff(t_fftdata(idx:qwer)));
% 		struct.slopiness(p) = struct.up(p) - struct.dwn(p);

%         % Get the amplitude of the first harmonic
% 		[~, harmonic_idx] = find(t_freqs >= 2*freq(p));
% 		harmonic_amp(p) = t_fftdata(harmonic_idx(1));

	% This is the time base for the trace
        	stw(p) = startwin + stp/2;

	% Advance one time step forward
        	p = p + 1; startwin = startwin + stp/2;

end

	struct.freqtrace = medfilt1(freq(1:end-2),5);
	struct.freqtim = stw(1:end-2);
	struct.amp = peakamp(1:end-2);
%	struct.harmonic_amp = harmonic_amp(1:end-2);

tt = find(tinytim > 0 & tinytim < 4);
% figure(1); wrengram(sig(tt),Fs);
% figure(1); hold on; plot(struct.freqtim, struct.freqtrace, 'b', 'LineWidth', 2); hold off;
% figure(1); hold on; plot(struct.freqtim, struct.amp*10000, 'g', 'LineWidth', 2); hold off;

%[xx yy] = ginput(1);
%	xxx = find(struct.freqtim >= xx);
%	figure(2); plot(struct.freqs{xxx(1)}, struct.fft{xxx(1)});

doloop = 0.008;
while doloop ~= 99;
			
	% Make the slope steeper so the threshold is easier to set
	struct.ampex = 100000*medfilt1(struct.amp,7).^5;
	struct.sylidx = find(struct.ampex > doloop);
	len = length(struct.freqtim);

	% zz extracts the traces during syllables
	% The "base" is -10000 so that there are clear syllable boundaries in the diff
	zz = ones(len,1)*-10000;
	zz(struct.sylidx) = struct.freqtrace(struct.sylidx);
	
	% The diff yields syllable starts (>+10000) and ends (< -8000)
	dizzy = diff(zz);

	% Get start and end markers
		start_idx = find(dizzy > 9000);
		end_idx = find(dizzy < -6000);

	figure(2); wrengram(sig(tt),Fs); hold on; plot(struct.freqtim, zz, 'g', 'LineWidth', 2); 

for i = 1:length(start_idx); plot([stw(start_idx(i)) stw(start_idx(i))], [4000 5500], 'y'); end;
for i = 1:length(end_idx); plot([stw(end_idx(i)) stw(end_idx(i))], [4000 5500], 'r'); end;
	
	hold off;
    
    fprintf('Threshold is %f. \n', doloop);
	doloop = input('Enter revised threshold (when done, enter 99): ');

end;

	struct.onlysyls = zz;	
	struct.dizzy = dizzy;
	out.startims = stw(start_idx);
	out.endtims = stw(end_idx);
    out.tims = stw;

	clear peakfreq tinytim stw;




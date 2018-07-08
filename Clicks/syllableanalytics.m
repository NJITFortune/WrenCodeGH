function struct = syllableanalytics(in);

for i = 1:length(in.sylen);

	sy = in.song(in.ind(i,1):in.ind(i,2));
	sy = sy - mean(sy);

% get frequency content of the syllable

	fftdata = fft(sy);
	struct.spectrum{i} = abs(real(fftdata(1:round(end/2))));
	stepsize = in.Fs/round(length(sy));
	struct.fftfreqs{i} = stepsize:stepsize:in.Fs/2;

	mn = min([length(struct.spectrum{i}) length(struct.fftfreqs{i})]);
	struct.spectrum{i} = struct.spectrum{i}(1:mn);
	struct.specfilt{i} = medfilt1(struct.spectrum{i},20);
	struct.fftfreqs{i} = struct.fftfreqs{i}(1:mn);

	[gg hh] = max(struct.specfilt{i});
	struct.peakfreq(i) = struct.fftfreqs{i}(hh);

	thresh = find(struct.specfilt{i} > gg/4);
	struct.minfreq(i) = struct.fftfreqs{i}(thresh(1));
	struct.maxfreq(i) = struct.fftfreqs{i}(thresh(end));
	struct.band(i) = struct.maxfreq(i) - struct.minfreq(i);

% Trace the frequency of the syllable

	siglen = length(sy);
	tinytim = 1/in.Fs:1/in.Fs:siglen/in.Fs;
	stp = 0.005;
	startwin = 0;
	p = 1;
	nfft = 2048;

	while startwin < max(tinytim) + stp;

        	data = sy(find(tinytim > startwin & tinytim < startwin + stp));
        	L = length(data);

        	t_fftdata = fft(data,nfft)/L;
        	t_fftdata = 2 * abs(t_fftdata(1:nfft/2+1));
        	t_freqs = in.Fs/2*linspace(0,1,nfft/2+1);

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

[q,w,e,r,t,y,u,k,o,p]=deriva(sy,in.Fs);

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


function out = syldata(data, Fs)
% out = syldata(data, Fs);
% Give this script 1 syllable worth of data - no more, no less.
%
% For internal use only! This script supports hagaclics and dualclics.
%
% PROVIDES:
% sylen         Length of syllable in seconds
% loud          Total magnitude
% meanloud      Time normalized magnitude
% spec_p        Power from FFT of the whole syllable
% spec_f        Frequencies for the FFT
% spec_pfilt    A smoothed version of the spec_p
% spec_peakf    Max frequency of FFT
% spec_minf     Minimum freq above amplitude threshold
% spec_maxf     Maximum freq above amplitude threshold
% spec_band     Just the differences between maxf and minf
% trace_freq    Trace of the loudest frequency in syllable
% trace_tim     Time base for trace_freq
% trace_slopemean   mean(diff( of the trace ))
% trace_slopestd    std(diff( of the trace ));
% trace_slopevar    var(diff( of the trace ));
% trace_peakf       max( of the trace );
% trace_peakt       Time in seconds of the peak frequency
% trace_peakpt      Time as a percentage of the lenght of the syllable

% Perhaps the most useful bit - the duration in seconds

    out.sylen = length(data)/Fs;

% Get loudness - which is kinda ridiculous, but kinda not

    nodcdata = data - mean(data);

    out.loud = sum(abs(nodcdata));
	out.meanloud = out.loud / out.sylen;

% Get frequency content of the syllable - FFT, no temporal information

	fftdata = fft(nodcdata);
	out.spec_p = abs(real(fftdata(1:round(end/2))));
	stepsize = Fs/round(length(nodcdata));
	out.spec_f = stepsize:stepsize:Fs/2;

	mn = min([length(out.spec_p) length(out.spec_f)]);
	out.spec_p = out.spec_p(1:mn);
	out.spec_f = out.spec_f(1:mn);
	out.spec_pfilt = medfilt1(out.spec_p,20);

	[gg, hh] = max(out.spec_pfilt);
	out.spec_peakf = out.spec_f(hh);

	thresh = find(out.spec_pfilt > gg/4);
	out.spec_minf = out.spec_f(thresh(1));
	out.spec_maxf = out.spec_f(thresh(end));
	out.spec_band = out.spec_maxf - out.spec_minf;

% Trace the frequency of the syllable

	siglen = length(nodcdata);
	tinytim = 1/Fs:1/Fs:siglen/Fs;
%	stp = 0.002; % old
    stp = 0.005; % Testing new
	startwin = 0;
	p = 1;
%	nfft = 2048; % old
    nfft = 512; % Testing new
    
	while startwin < max(tinytim) + stp;
            % Matlab removed the find command here
        	ttdata = nodcdata(tinytim > startwin & tinytim < startwin + stp);
        	L = length(ttdata);

        	t_fftdata = fft(ttdata,nfft)/L;
        	t_fftdata = 2 * abs(t_fftdata(1:nfft/2+1));
        	t_freqs = Fs/2*linspace(0,1,nfft/2+1);

        	[maxval, idx] = max(t_fftdata);
        	peakamp(p) = maxval;
        	peakfreq(p) = t_freqs(idx);
        	stw(p) = startwin + stp/2;

        	p = p + 1; startwin = startwin + stp/2;

	end

	out.trace_freq = medfilt1(peakfreq(1:end-2),5);
    out.trace_amp = medfilt1(peakamp(1:end-2),5);
    out.trace_tim = stw(1:end-2);
    
	out.trace_slopemean = mean(diff(out.trace_freq));
	out.trace_slopestd = std(diff(out.trace_freq));
	out.trace_slopevar = var(diff(out.trace_freq));
	[out.trace_peakf, indf] = max(out.trace_freq);
	out.trace_peakt = stw(indf);
	out.trace_peakpt = out.trace_peakt / stw(end);


end



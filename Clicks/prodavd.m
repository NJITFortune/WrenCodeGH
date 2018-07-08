function struct = prodavd(old);


struct = old;

for i = 1:length(old.syl);

	sy = struct.song(struct.ind(i,1):struct.ind(i,2));
	sy = sy - mean(sy);

	% get loudness
	struct.loud(i) = sum(abs(sy));
	struct.meanloud(i) = struct.loud(i) / struct.sylen(i);

	% get freq
	fftdata = fft(sy);
	struct.spectrum{i} = abs(real(fftdata(1:round(end/2))));
	stepsize = old.Fs/round(length(sy));
	struct.fftfreqs{i} = stepsize:stepsize:old.Fs/2;
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

% Use Ofer's student analysis here

[q,w,e,r,t,y,u,k,o,p]=deriva(sy,old.Fs);

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


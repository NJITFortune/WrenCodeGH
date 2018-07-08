function out = corrpicker(struct);

Fs = struct.Fs;
	pz = zeros(100,1);

	ref = [pz' struct.freqtrace{4} pz'];

	lref = length(ref);
	zz = zeros(lref,1);
	mh = max(ref);

for i = 1:length(struct.freqtrace);

	tt(1,:) = ref; tt(2,:) = zz;

		addzero = floor( (lref - length(struct.freqtrace{i})) /2);
		az = zeros(addzero, 1);
		compsig = [az' struct.freqtrace{i} az'];
		sf = max(compsig) / mh;
		compsig = compsig / sf;
	max(compsig)
	mh

		tt(2,1:length(compsig)) = tt(2,1:length(compsig)) + compsig;

	out.xc{i} = xcorr( tt(1,:) - mean(tt(1,:)), tt(2,:) - mean(tt(2,:)));
	[out.maxcorr(i) ind] = max(out.xc{i});
	out.percorr(i) = (ind / length(zz)) - 0.5;

	clear tt
end

plot(out.percorr, out.maxcorr, '*');
 

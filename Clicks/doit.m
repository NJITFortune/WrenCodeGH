j=8;
figure(1); hold on;

	for i = 1:length(foo(j).num);
		a = foo(j).num(i);
		plot(out(a).sylen, out(a).trace_peakf, 'r*');
	end;



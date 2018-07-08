function [burstcount isocount totalcounts bf] = burstfractionspont(spikes, isithresh, boundaries);

counter = 0;

startime = boundaries(1);
endtime = boundaries(2);

cells = length(spikes);

for i = 1:cells

burstcount(i) = 0;
isocount(i) = 0;

	reps = length(spikes(i).bospikes);
	[burstspikes isospikes] = burst_extractor( spikes(i).bospikes, isithresh );

	for j = 1:reps

		counter = counter+1;
		burstcount(i) = burstcount(i) + length(burstspikes{j});
		isocount(i) = isocount(i) + length(isospikes{j});

	end

totalspikes(i) = burstcount(i) + isocount(i);
bf(i) = burstcount(i) / totalspikes(i);


clear burstspikes isospikes

end



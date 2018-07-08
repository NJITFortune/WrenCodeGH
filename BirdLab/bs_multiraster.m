function bs_multiraster(spikedata)
% bs_raster(spikedata)
% spikedata is the multi cell array from bs_multiconverter
% e.g. w(1).multi
% Version 8 October 2016 - added multi capability
% Version 18 March 2015, VERSION 1 November 2008, original JAB 7/25/07

clrs(1,:)='b-'; clrs(2,:)='r-'; clrs(3,:)='m-'; 
clrs(4,:)='g-'; clrs(5,:)='c-'; clrs(6,:)='k-'; 
clrs(7,:)='y-';
clrs(8,:)='b-'; clrs(9,:)='r-'; clrs(10,:)='m-'; 
clrs(11,:)='g-'; clrs(12,:)='c-'; clrs(13,:)='k-'; 
clrs(14,:)='y-';

hold on;

for spk = 1:length(spikedata);

for thistrial = 1:length( spikedata(spk).spikes );
   for ss = 1:length( spikedata(spk).spikes{thistrial});
      plot( ones(1,2).*spikedata(spk).spikes{thistrial}(ss), [thistrial-0.3 thistrial+0.3], clrs(spk,:),'LineWidth', 0.5); 
      % plot( ones(1,2).*spiketimes{tt}(ss), [tt-0.3 tt+0.3], 'k', 'LineWidth', 1.5);
   end;
end;

ax = axis;
axis( [ax(1:2) 0 length( spikedata )+1] )



end

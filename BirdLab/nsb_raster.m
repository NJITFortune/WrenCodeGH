function nsb_raster( spiketimes, fig_handle )
% nsb_raster( spiketimes, fig_handle )
%
% SPIKETIMES should be a cell array, containing spike times
% for each trial
% FIG_HANDLE is an optional argument containing the figure number
% in which to plot the raster
%
% JAB 7/25/07

if exist( 'fig_handle', 'var' )
   figure( fig_handle ); clf
else
   figure
end
hold on

for tt = 1:length( spiketimes )
   for ss = 1:length( spiketimes{tt} )
      plot( ones(1,2).*spiketimes{tt}(ss), [tt-0.3 tt+0.3] )
   end
end

ax = axis; axis( [ax(1:2) 0 length( spiketimes )+1] )

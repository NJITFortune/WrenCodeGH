function bs_raster(spiketimes, col)
% bs_raster(spiketimes, col)
% spiketimes is the cell array from bs_converter
% col is the color, e.g. 'b' or 'm'
% Version 18 March 2015, VERSION 1 November 2008, original JAB 7/25/07

if nargin < 2; col = 'k'; end % Default to black tick marks

hold on;

for tt = 1:length( spiketimes )
   for ss = 1:length( spiketimes{tt} )
      plot( ones(1,2).*spiketimes{tt}(ss), [tt-0.3 tt+0.3], col, 'LineWidth', 0.5);
      % plot( ones(1,2).*spiketimes{tt}(ss), [tt-0.3 tt+0.3], 'k', 'LineWidth', 1.5);
   end
end

ax = axis;
axis( [ax(1:2) 0 length( spiketimes )+1] )



end

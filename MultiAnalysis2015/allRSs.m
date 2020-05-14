function [r, max_data] = allRSs( spikes, stimtime, spontime, fig_handle, x, y )
% [Zs, max_data] = allZs( spikes, stimtime, spontime, fig_handle )
%
% SPIKES is a cell array containing spike times for each repetition of the stimulus
% STIMTIMEs define the period of each trial during which stimulation occurred,
%   sent as a vector: [STARTSTIMTIME ENDSTIMTIME]
% SPONTIMEs are the period of each trial to use for calculating spontaneous reponse,
%   in the same format as STIMTIME
% if the optional argument FIG_HANDLE is sent, the Z-plot is generated in the
%   corresponding figure, else no plot is made
%
% return argument Zs is the Z-values at each combination of window size and
%   time offset
% MAX_DATA is a vector containing 3 elements: the maximum Z-score achieved,
%   the window size, and the offset time at which the maximum score occured
% 
% JAB 7/23/07

% try to figure out input arguments
if nargin > 4 & numel( stimtime ) == 1 % probably called wrong, by EF
   startstimtime = stimtime;
   endstimtime = spontime;
   startspontime = fig_handle;
   endspontime = x;
   if nargin == 4
      fig_handle = y;
   end
elseif numel( stimtime ) == 2
   startstimtime = stimtime(1);
   endstimtime = stimtime(2);
   startspontime = spontime(1);
   endspontime = spontime(2);
else
   error( 'check input arguments more carefully, jackass' )
end

% create time vectors
time_step = 0.010; % ms
time_bins = startstimtime:time_step:endstimtime;
% spon_bins = startspontime:time_step:endspontime;

% create Z-window-size vector
min_wind_width = 0.100; % ms
max_wind_width = 2;

wind_step = 0.050;
wind_bins = min_wind_width:wind_step:max_wind_width;

r.RS = zeros( length(wind_bins), length(time_bins) );
    r.Xs = time_bins;
    r.Ys = wind_bins;
    
% for each time bin, during stimulus presentation
for tt = 1:length( time_bins )

   % calculate Z-scores at all usable window sizes
   
   use_wind = wind_bins(wind_bins+time_bins(tt) < endstimtime );
   % length(use_wind)
   for ww = 1:length( use_wind )

      r.RS(ww,tt) = rs( spikes, ...
         time_bins(tt), time_bins(tt) + use_wind(ww),...
         startspontime, endspontime);

   end % for each Z-window size

end % for each time bin

% shift responses into center of matrix
for ww = 2:length( wind_bins )
   shift_size = round( wind_bins(ww)/2 / time_step );
   r.RS(ww,:) = [r.RS(ww,end-shift_size+1:end) r.RS(ww,1:end-shift_size)];
end % for each Z-window size

% find indices of maximum Z-score
    [m,i] = max( r.RS, [], 1 );
    [max_data(1), max_data(3)] = max( m );
    max_data(2) = i(max_data(3));
    
    max_data(3) = r.Xs(max_data(3));
    max_data(2) = r.Ys(max_data(2));

% plot results
if exist( 'fig_handle', 'var' )
   figure( fig_handle ); clf
   ax(1) = subplot(211);
       surf( r.Xs, r.Ys, r.RS )

       xlabel( {'Window center (s)'} )
       ylabel( 'Window size (s)' )
       zlabel( 'Z-score' )

%    % label ticks with correct values
%    xt = get( gca, 'xtick' );
%    xt = xt - (startstimtime + time_step/2) / time_step; % put ticks at (at least) integer values
%    set( gca, 'xtick', xt )
%    set( gca, 'xticklabel', round( xt .* time_step + startstimtime + time_step/2 ) )
%    
%    yt = get( gca, 'ytick' );
%    yt = (yt(1):(yt(2)-yt(1))/2:yt(end)); % subdivide y ticks
%    yt = yt - min_wind_width/wind_step;
%    set( gca, 'ytick', yt )
%    set( gca, 'yticklabel', yt .* wind_step + min_wind_width ) 
   
   % set color and view properties
   grid off
   shading interp
   colormap jet
%    ax_min = -2; % default minimum plot value
%    ax_max = 3; % default max
%    axis( [min( [xt(1) ax(1)] ) max( [xt(end) ax(2)] ) ...
%       min( [yt(1) ax(3)] ) max( [yt(end) ax(4)] ) ...
%       min( [ax(5) ax_min] ) max( [ax(6) ax_max] )] )
%   caxis( [min( [ax(5) ax_min] ) max( [ax(6) ax_max] )] )
   set( gca, 'view', [0 90] )
   
   hold on;
   plot3(max_data(3), max_data(2), max_data(1), 'ko');
   text(max_data(3)+0.2, max_data(2), max_data(1), num2str(max_data(1)));
   text(max_data(3)+0.2, max_data(2)-0.075, max_data(1), num2str(max_data(3)));
   text(max_data(3)+0.2, max_data(2)-0.15, max_data(1), num2str(max_data(2)));
   
   ax(2) = subplot(212);
   bs_raster(spikes);
   linkaxes(ax, 'x');
   %caxis([0 5]);
   
   
end

    function qwe = rs(spks, startim, endtim, sponstart, sponstop)
        stimspks = 0; sponspks = 0;
        for jj = 1:length(spks)
           stimspks = stimspks + length(find(spks{jj} > startim & spks{jj} < endtim));
           sponspks = sponspks + length(find(spks{jj} > sponstart & spks{jj} < sponstop));            
        end
        
       qwe = (stimspks/((endtim-startim)*length(spks))) - (sponspks/((sponstop-sponstart)*length(spks)));
       
    end

end


 



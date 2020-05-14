function [zScore]=nsbZfinder(spiketimes, startstim, endstim, startspon, endspon)

% Find the durations (which are specified in msec) and convert to seconds.
StimDur=(endstim-startstim);
SponDur=(endspon-startspon);

% Loop through each trial and count spikes
for i=1:length(spiketimes)

   StimCount(i) = length( find( spiketimes{i} > startstim & spiketimes{i} < endstim ) );
   SponCount(i) = length( find( spiketimes{i} > startspon & spiketimes{i} < endspon ) );
   % --JB
   
%    % Get the indices for spikes in the area of interest
% 	StimRes=(spiketimes{i} > startstim & spiketimes{i} < endstim) .* spiketimes{i};
% 
%    % Get the indices for spikes in the spontaneous region
% 	SponRes=(spiketimes{i} > startspon & spiketimes{i} < endspon) .* spiketimes{i};
% 
%    % get the actual times and eliminate zeros
%    spikesinStimrange{i}=find(StimRes);
%    spikesinSponrange{i}=find(SponRes);
% 
%    % Get the rates	
% 	StimCount(i)=length(spikesinStimrange{i})/StimDur;
% 	SponCount(i)=length(spikesinSponrange{i})/SponDur;

end

StimRate = StimCount ./ StimDur;
SponRate = SponCount ./ SponDur;
% --JB

% This is the calculation for the Z-score

sm = mean(SponRate); % spont mean rate
sv = var(SponRate);  % spont rate variance
rm = mean(StimRate); % resp mean rate
rv = var(StimRate);  % resp rate variance
cv = cov(StimRate,SponRate); % ?? apparently does something different in Octave...
if size( cv, 1 ) ~= 1, cv = cv(1,2); end % --JB: rate covariance btwn spont and resp

denom = real( sqrt(rv + sv - 2*cv) );

zScore = (rm - sm) / max( [denom 1e-3] ); % to avoid dividing by zero


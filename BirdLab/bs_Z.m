function zScore = bs_Z(spiketimes, startstim, endstim, startspon, endspon)
% zScore = bs_Z(spiketimes, startstim, endstim, startspon, endspon);
% zScore = bs_Z(spiketimes, stimtim);
% Find the durations (which are specified in sec) and convert to seconds.
% Two ways to call it
% VERSION 18 August 2008, 26 March 2015

% if nargin == 2
%     stimtim = startstim;
% 
%     startstim = stimtim(1);
%     endstim= stimtim(2);
%     startspon = stimtim(3);
%     endspon = stimtim(4);
% end

StimDur=(endstim-startstim);
SponDur=(endspon-startspon);

%% Loop through each trial and count spikes
for i=length(spiketimes):-1:1

   StimPass(i) = length( find( spiketimes{i} > startstim & spiketimes{i} < endstim ) );
   SponPass(i) = length( find( spiketimes{i} > startspon & spiketimes{i} < endspon ) );

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
% 	StimPass(i)=length(spikesinStimrange{i})/StimDur;
% 	SponPass(i)=length(spikesinSponrange{i})/SponDur;

end

StimPass = StimPass ./ StimDur;
SponPass = SponPass ./ SponDur;

% --JB
% This is the calculation for the Z-score

sm = mean(SponPass); % spont mean
sv = var(SponPass);  % spont variance
rm = mean(StimPass); % resp mean
rv = var(StimPass);  % resp variance
cv = cov(StimPass,SponPass); % ?? apparently does something different in Octave...

    if size( cv, 1 ) ~= 1, cv = cv(1,2); end % --JB

denom = real( sqrt(rv + sv - 2*cv) );

zScore = (rm - sm) / max( [denom 1e-10] ); % to avoid dividing by zero

end

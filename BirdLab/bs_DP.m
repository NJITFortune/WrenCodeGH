function dPrime = bs_DP( spiketimes1, stimwind1, spiketimes2, stimwind2 )
% dPrime = nsbDPfinder( spiketimes1, stimwind1, spiketimes2, stimwind2 )
% VERSION 18 August 2008
%
% SPIKETIMES are cell arrays containing vectors of the spike times in response
%   to each presentation of the two stimuli
% STIMWINDs are 2x1 vectors containing the [START,STOP] times of the stimulus,
%   to be used in determining which spikes in each trial are stimulus-evoked
%
% JAB 7/27/07

% Find the durations (which are specified in sec) and convert to seconds.
stim1dur = -diff( stimwind1 );
stim2dur = -diff( stimwind2 );

stim1count = zeros( 1, length( spiketimes1 ) );
stim2count = zeros( 1, length( spiketimes2 ) );

reps=min([length(spiketimes1) length(spiketimes2)]);

% Loop through each trial and count spikes
for i=1:reps
   stim1count(i) = length( find( spiketimes1{i} > stimwind1(1) & spiketimes1{i} < stimwind1(2) ) );
   stim2count(i) = length( find( spiketimes2{i} > stimwind2(1) & spiketimes2{i} < stimwind2(2) ) );
end

stim1rate = stim1count ./ stim1dur;
stim2rate = stim2count ./ stim2dur;

% This is the calculation for the d-prime measure

m1 = mean( stim1rate );
m2 = mean( stim2rate );
v1 = var( stim1rate );
v2 = var( stim2rate );

denom = sqrt( v1 + v2 );

dPrime = (m1 - m2) / max( [denom 1e-3] ); % to avoid dividing by zero


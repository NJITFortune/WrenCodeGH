function out = bs_calculator(struct, stimnum, plotwin, burstthresh, tts)    
% out = bs_calculator(spikes, stims, Fs, win, stimtim, threshold);
% struct is the structure at the neuronal unit level
% stimnum is the stimulus ID
% plotwin is what to show - not to be confused with the times for analysis
% tts is option and contains start and end for stimulus followed by same
% for sponteous.
% Version 31 March 2015, Version 21 July 2010, VERSION 18 August 2008

clf

if nargin < 5; 

	bs_physplot(struct, stimnum, plotwin, 200, 0, 0);
%bs_physplot(struct, stimnum, win, binwidth, plt_type, prt)
    
	% Get the click for the start of the stimulus...
	disp('Click to get the start of the stimulus (or the start of the response, if you prefer).');
	[stimstart, ~] = ginput(1);

	subplot(4,1,4); hold on; plot([stimstart stimstart], [-1 1], 'g');

	% Get the click for the end of the stimulus...
	disp('Click to get the end of the stimulus (or the end of the response, if you prefer).');
	[stimend, ~] = ginput(1);

	subplot(4,1,4); hold on; plot([stimend stimend], [-1 1], 'r');

	% Get the end of the spontaneous
	disp('Click to get the end of the SPONTANEOUS activity.');
	[sponend, ~] = ginput(1);

    dur = stimend - stimstart;
	sponstart = sponend - dur;

	out.clicks = [stimstart stimend sponstart sponend];

else
	stimstart = tts(1); stimend = tts(2); sponstart = tts(3); sponend = tts(4);
end

dur = stimend - stimstart;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 6 
	thresh = burstthresh;
else
	% Separate Bursts and Isos

	clf
    histola = bs_isiHist(struct.stim(stimnum).spikes, 60);

	% Get the separation point for burst and isolated spikes
	disp('Click where you want to separate bursts and isolated spikes.');
	[thresh, ~] = ginput(1);
end
    
	[bursts, isos] = bs_extractor(struct.stim(stimnum).spikes, thresh);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out.all.stimcount=0; out.all.sponcount=0;
out.burst.stimcount=0; out.burst.sponcount=0;
out.iso.stimcount=0; out.iso.sponcount=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get total spikes for response and spontaneous
for i = 1:length(struct.stim(stimnum).spikes)
    
out.all.stimcount = out.all.stimcount + length( find (struct.stim(stimnum).spikes{i} > stimstart & struct.stim(stimnum).spikes{i} < stimend) );
out.all.sponcount = out.all.sponcount + length( find (struct.stim(stimnum).spikes{i} > sponstart & struct.stim(stimnum).spikes{i} < sponend) );
out.burst.stimcount = out.burst.stimcount + length( find (bursts{i} > stimstart & bursts{i} < stimend) );
out.burst.sponcount = out.burst.sponcount + length( find (bursts{i} > sponstart & bursts{i} < sponend) );
out.iso.stimcount = out.iso.stimcount + length( find (isos{i} > stimstart & isos{i} < stimend) );
out.iso.sponcount = out.iso.sponcount + length( find (isos{i} > sponstart & isos{i} < sponend) );

end

out.all.stimrate = out.all.stimcount / (dur * length(struct.stim(stimnum).spikes) );
out.all.sponrate = out.all.sponcount / (dur * length(struct.stim(stimnum).spikes) );
% ResponseStrength = ResponseDuringStimulusSpikesPerSecond - SpontSpikesPerSecond
out.all.RespStrength = out.all.stimrate - out.all.sponrate;

%%% Now do the same for burst and isolated spikes

out.burst.stimrate = out.burst.stimcount / (dur * length(struct.stim(stimnum).spikes) );
out.burst.sponrate = out.burst.sponcount / (dur * length(struct.stim(stimnum).spikes) );
out.burst.RespStrength = out.burst.stimrate - out.burst.sponrate;
out.iso.stimrate = out.iso.stimcount / (dur * length(struct.stim(stimnum).spikes) );
out.iso.sponrate = out.iso.sponcount / (dur * length(struct.stim(stimnum).spikes) );
out.iso.RespStrength = out.iso.stimrate - out.iso.sponrate;

% Z-scores

out.all.z = bs_Z(struct.stim(stimnum).spikes, stimstart, stimend, sponstart, sponend);
out.burst.z = bs_Z(bursts, stimstart, stimend, sponstart, sponend);
out.iso.z = bs_Z(isos, stimstart, stimend, sponstart, sponend);

% plot for fun

t = linspace(plotwin(1),plotwin(2),length(struct.stim(stimnum).stim));
subplot(4,1,1); plot(t,struct.stim(stimnum).stim); ylabel('Stimulus')
trimaxes(plotwin);
subplot(4,1,2); bs_raster(struct.stim(stimnum).spikes); ylabel('All spikes')
trimaxes(plotwin);
subplot(4,1,3); bs_raster(bursts); ylabel('Burst spikes')
trimaxes(plotwin);
subplot(4,1,4); bs_raster(isos);   ylabel('Isolated spikes'); 
trimaxes(plotwin);
xlabel('time (sec)');

% add the raw data to this structure... why not?
out.burspikes = bursts;
out.isospikes = isos;
out.allspikes = struct.stim(stimnum).spikes;


function trimaxes(win)
ax = axis;
ax(1:2) = win;
axis(ax)


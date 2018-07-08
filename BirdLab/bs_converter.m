function zf = bs_converter(spiketimes, kbd, stim, limits)
% function zf = bs_converter(SpikeTimes, KeyCodes, Stimulus, limits)
% BS_CONVERTER Converts a Spike2 matlab file into a cell array for use
% with the BS array of analysis tools.
% Version 18 March 2015; VERSION 21 July 2010; VERSION 18 August 2008
% Usage:
% 1. Load File that you have save from Spike2: e.g., "load leftHVCdata.mat"
% 2. Figure out what structures have these data (use "whos" command):
%    Spike times (e.g. leftHVCdata_Ch6)
%    Keyboad codes (usually Ch31 - e.g. leftHVCdata_Ch31)
%    Stimulus recordings (e.g. leftHVCdata_Ch4)
%
% 3. Run this command:
%    e.g. zf = bs_converter(yourspikes_Ch6, yourkeycodes_Ch31, yourstim_Ch4, [10 15]);
%
%    The limits ([10 15])are the pre- and post keycode time (both are positive!!!!!)
%    In this example, it would take 10 seconds before the stimulus onset
%    and 15 seconds after the stimulus onset. Make sure that these are
%    larger than you expect to need.
%
% 4. The script will ask you some questions.  Answer honestly - this is for posterity.
%
% The output structure has the following:
% 
% zf(#).spikes is a cell array with the spike times for each repetition of
% the stimulus. # is the index of the stimulus (e.g. BOS or REV). 
% So, zf(4).spikes{3} would give you the 3 repetition of the 4th stimulus.
%
% zf.Fs is samplerate of stimulus, zf.StimName is a string name of stimulus.
% zf.pretrig/postrig is the time in seconds before/after the onset of stim.
% zf.stim is the stimulus and zf.tim is a time vector for it.

%% Initial Questions

  bird = input('Bird identity: ','s');
  sex = input('Bird sex: ','s');
  dt = input('Experiment date (YYYYMMDD): ','s');
  physsite = input('Recording Site: ','s');
  unitnum = input('Unit Number: ','s');


%% Prep work %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get spiketimes - the user could give us either the structure from Spike2
% or a list of spiketimes. I expect that it will usually be the Spike2 structure.

    if isstruct(spiketimes); 
    		spikeTimes = spiketimes.times; % structure from Spike2
    else
    		spikeTimes = spiketimes; % a list of spike times
    end;

	zf.allspiketimes = spikeTimes; % Add this to the output structure

% Get stimulus codes and times from the "keyboard"

    keycodes = kbd.codes;
	keytimes = kbd.times;

% Get the list keycodes in the file - usually start with 49.

    StimMin = min(keycodes(:,1)); % first stim (usually 49)
    StimMax = max(keycodes(:,1)); % last stim
        StimMinLessOne = StimMin - 1; % Spike2 starts at 49 - but I want to start at 1.
    stimulist = (StimMin-StimMinLessOne):(StimMax-StimMinLessOne); % list of stimuli e.g. 49 50 51 52 now 1 2 3 4
    stimnum = length(stimulist); % number of stimuli that were used in this experiment.

% Get the samplerate of the stimulus

	Fs =  1/stim.interval;

% Set the pre- and post- trigger values (this step is just to make the code
% easier to read).
    
    pretrig = limits(1);
    posttrig = limits(2);

% Give the user some feedback - how many stimuli do we have and what are their numbers.
sprintf('There were a total of %d stimuli.',stimnum)
disp('The stimuli are:')
disp(stimulist)
clf

%% Pre-Allocations for speed
    for i = stimnum:-1:1;
        zf(i).stim = zeros(1,(pretrig+posttrig)*Fs);
        zf(i).tim = 1/Fs:1/Fs:length(zf(i).stim)/Fs; % a time vector 
        zf(i).tim = zf(i).tim - pretrig; % adjust so that zero is the start of the stimulus.
        zf(i).StimNum = i;
        zf(i).pretrig = pretrig;
        zf(i).posttrig = posttrig;
        zf(i).Fs = Fs;
        zf(i).stimend = [];
        zf(i).sex = sex;
        zf(i).birdname = bird;
        zf(i).date = dt;
        zf(i).site = physsite;
        zf(i).unit = unitnum;
    end;


%% This is the big loop for each stimulus in the experiment %%%%%%% 

for i = 1:stimnum;
    
    % Extract the stimulus so that we can show it to the user who will
    % input its name

    % Get the start times for the current stimulus
    Starts = keytimes(keycodes == (StimMinLessOne+stimulist(i)));

    % Set the window around the start of the stimulus to extract the
    % stimulus. We use the second iteration of the stimulus for no
    % particular reason "Starts(2)".
    windowstart = round( ((Starts(2) - pretrig) * Fs) );
    windowend =  round( ((Starts(2) + posttrig) * Fs) );

    stimulus = stim.values(windowstart:windowend); % The stimulus data
    zf(i).stim = stimulus(1:length(zf(i).stim)); 

    % Plot the current stimulus
    figure(1);
    subplot(2,1,1); specgram(zf(i).stim,1024,Fs);
    subplot(2,1,2); plot(zf(i).tim,zf(i).stim);
        xlim([-pretrig posttrig]);

        % Ask the user to click on the end of the stimulus
        sprintf ('Click on the END of the stimulus in BOTTOM panel.')
        [zf(i).stimend, ~] = ginput(1);
        sprintf ('End of stimulus is %0.5g seconds.', zf(i).stimend)

        % Ask the user to identify the stimulus
        sprintf ('Enter a label for this stimulus %d :', stimulist(i))
        zf(i).StimName = input('-----> : ','s');

 
%% This subloop extracts the spikes for each stimulus repetition and puts them into the cell array
    
    zf(i).spikes = cell(length(Starts),1);
    
    for j=1:length(Starts)
        spiketimes = spikeTimes > (Starts(j) - pretrig) & spikeTimes < (Starts(j) + posttrig);
        %         spiketimes = find(spikeTimes > (Starts(j) - pretrig) & spikeTimes < (Starts(j) + posttrig));
        zf(i).spikes{j} = spikeTimes(spiketimes) - Starts(j);
        %        zf(i).spikes{j} = spikeTimes(spiketimes) - Starts(j);
    end
        
end;


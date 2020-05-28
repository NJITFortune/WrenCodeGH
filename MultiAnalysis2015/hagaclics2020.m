function struct = hagaclics2020(sng, Fs)
% out = hagaclics(sng, Fs);
% sng is the sample data - straight from the wav file
% Fs is the sample rate in Hz
%
% This function depends on oscson for plotting.
% This function depends on syldata for analysis.
% This function depends on the SAMIII tree via syldata.
% This function works with slicer for selecting syllable types.
% This function depends on ginputc for clicking (can be reverted to ginput easily).

%% Housekeeping

% Make and time sequence and get the length of the signal.
    tim = 1/Fs:1/Fs:length(sng)/Fs;

% This is the width of the clicking window in seconds. 
    windwid = 2.5;

% This is the overlap in seconds between previous and next windows.
    ovrlp = 0.2;

% High-pass filter the sng

    [b,a] = butter(2,200/(Fs*2), 'high');
    sng = filtfilt(b,a,sng);

%% Initialize BSX and SYLNUM and PRECLICK

% For analysis we start at time 0 (base X or bsx) and syllable #1.
    bsx = 0; sylnum = 1; preclick = 0;

%% Make initial PROGRESS plot
figure(2); clf; hold on;
    plot(tim,sng,'b');
    xlim([0 tim(end)]);

%% Get clicks: Loop until user tells us that this is the end of sng
cntu = 1;

while cntu < 10
       
% Get clicks

    tt = find(tim >= bsx & tim < bsx + windwid); % The time window for clicking

    figure(2); % Update our current position by plotting a vertical magenta line
        plot([tim(tt(end)) tim(tt(end))], [-0.8 0.8], 'm', 'LineWidth', 2);        
   
     % Get clicks with the embedded function clickplotter   
     cts = clickplotter(sng(tt), Fs, preclick);

% We loop to analyze each syllable

        for ss = 1:2:length(cts)-1

            % Get the region of data for the current click from the data
            stt = find(tim >= bsx + cts(ss) & tim < bsx + cts(ss+1));

            % Send the data to our analysis and tracing function
            tmp = syldata(sng(stt),Fs);
            
            % CLEANUP THE DATA
            
            % 1) Find and eliminate zero frequency
            
            if ismember(0,[tmp.trace_freq]) == 1

            % Get the indices where these zeros happened
                zips = find([tmp.trace_freq] == 0);

            % Most of the time there is only one zero value, usually at the
            % end - so we can easily fix that.

                if length(zips) == 1 
                    % If the zero is anywhere but the first value, we'll simply
                    % copy the previous value
                    if zips ~= 1 
                        tmp.trace_freq(zips) = tmp.trace_freq(zips-1);
                    else
                        tmp.trace_freq(zips) = tmp.trace_freq(zips+1);
                    end
                else
                    % But if we have more than a single zero value, we'll need
                    % to get user to help. Here we have the user click each zero.
                    NumberofUserClicks = length(zips);
                    fprintf('Multiple Zeros: N=%i - click each to fix\n', NumberofUserClicks)
            
                    figure(3); clf;
                    %wrengram(sng(tmp.sylind(1):tmp.sylind(2)), Fs);
                    specgram(sng(stt),1024,Fs,[],1000); ylim([-10 5000]);
                    lliimmss = xlim;
                    xlim([0 lliimmss(2)+0.1]);
                    hold on; 
                    plot(tmp.trace_tim, tmp.trace_freq,'*-', 'LineWidth', 2);
                    plot(tmp.trace_tim(zips), tmp.trace_freq(zips), 'ro');
            
                    clks = ginput(NumberofUserClicks);
                    tmp.trace_freq(zips) = clks(:,2);            
                end    

            end
            
% Place into the output structure.
 
                struct(sylnum).sylen = tmp.sylen;
                struct(sylnum).loud = tmp.loud;
                struct(sylnum).meanloud = tmp.meanloud;
                struct(sylnum).spec_p = tmp.spec_p;
                struct(sylnum).spec_f = tmp.spec_f;
                struct(sylnum).spec_pfilt = tmp.spec_pfilt;
                struct(sylnum).spec_peakf = tmp.spec_peakf;
                struct(sylnum).spec_minf = tmp.spec_minf;
                struct(sylnum).spec_maxf = tmp.spec_maxf;
                struct(sylnum).spec_band = tmp.spec_band;
                struct(sylnum).trace_freq = tmp.trace_freq;
                struct(sylnum).trace_amp = tmp.trace_amp;
                struct(sylnum).trace_tim = tmp.trace_tim;
                struct(sylnum).trace_slopemean = tmp.trace_slopemean;
                struct(sylnum).trace_slopestd = tmp.trace_slopestd;
                struct(sylnum).trace_slopevar = tmp.trace_slopevar;
                struct(sylnum).trace_peakf = tmp.trace_peakf;
                struct(sylnum).trace_peakt = tmp.trace_peakt;
                struct(sylnum).trace_peakpt = tmp.trace_peakpt;

        % Do Ofer's analysis

        [q,w,e,r,t,y,u,k,o,p]=deriva(sng(stt) - mean(sng(stt)),Fs);

                struct(sylnum).ofer_mSpecDeriv = q;
                struct(sylnum).ofer_mAM = w;
                struct(sylnum).ofer_mFM = e;
                struct(sylnum).ofer_mEntropy = r;
                struct(sylnum).ofer_mAmp = t;
                struct(sylnum).ofer_Gravity = y;
                struct(sylnum).ofer_mPitchGood = u;
                struct(sylnum).ofer_mPitch = k;
                struct(sylnum).ofer_PitchChose = o;
                struct(sylnum).ofer_PitchWeight = p;
                
                struct(sylnum).ofer_PW = mean(p);
                struct(sylnum).ofer_GC = mean(y);
                
         % add the fundamental components of the syllable

            struct(sylnum).syltim = [tim(stt(1)) tim(stt(end))];
            struct(sylnum).sylind = [stt(1) stt(end)];

            struct(sylnum).year = nfo.year; struct(sylnum).month = nfo.month;
            struct(sylnum).day = nfo.day; struct(sylnum).location = nfo.location;
            struct(sylnum).note = nfo.note; struct(sylnum).wavname = nfo.wavname;

            % PLOT ONTO OUR WINDOW
            
            figure(1); hold on;
            timstr = struct(sylnum).trace_tim + struct(sylnum).syltim(1);
            plot(timstr - bsx, struct(sylnum).trace_freq, 'g');
            hold off;

            % Advance our syllable number now

            sylnum = sylnum + 1;
            
        end
        
        
        % We are done clicking through the window, so let's save temp data
        % and reset bsx
        
        %!% save temporary.mat struct;
               
        bsx = struct(sylnum-1).syltim(2) - ovrlp;
        preclick = ovrlp;
        
        % Ask if we are done
        cntu = inputdlg('Done? 99=yes, 0 or [return] to continue: ');
        cntu = str2num(cntu{:});
        if isempty(cntu) == 1; cntu = 1; end
        
end

%%%%%%%%%%%%%%%%%%% END OF CLICKING


%% After all is said and done, get inter-syllable-intervals

for j = 2:length(struct)
    struct(j).ISI = struct(j).syltim(1) - struct(j-1).syltim(2);
end


function [ clicktimes ] = clickplotter( data, Fs, preCLK )
% clicktimes = clickplotter(data, Fs);
% data is a chunk of raw wren wrecorded data
% Fs is the sample rate in Hz
% This function plots the chunk of song and has the user click all of the
% syllable starts and ends in order. 
% requires OSCSON

    yvals = [100 4500];
    
%% Plot and get clicks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        figure(1); clf; hold on;
        
        figprop = get(gcf,'Position'); 
        set(gcf,'Position',[figprop(1) figprop(2) 1000 400]);

        oscson(data,Fs);

        % Plot prior click in the main window

        plot([preCLK preCLK], yvals, 'm'); 

    % Get clicks for the starts and ends of the syllables
    
        [clicktimes, ~] = ginputc('Color','w','ShowPoints',true,'ConnectPoints',false);        
        
end




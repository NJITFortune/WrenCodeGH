function [md, fd] = dualclicsB(malorig, femorig, Fs, td)
% [maledata femaledata] = dualclics(male, female, Fs);
% mal is the sample data - straight from the Male wav file
% fem is the sample data - straight from the Female wav file
% mal and fem should be exactly the same length
% td is the time difference between the recordings as determined by the
% initial clicks at the start of the recording. If the male recording is
% delayed relative to the female, td is positive. If the female is delayed
% relative to the male, td is negative. Perfectly aligned, use 0.
% Fs is the sample rate in Hz
%
% This function depends on oscson for plotting.
% This function depends on syldata for analysis.
% This function depends on the SAMIII tree via syldata.
% This function depends on sexd for determining the sex of each syllable.
% This function depends on ginputc for clicking (can be reverted to ginput easily).
% This function works with slicer for selecting syllable types.

%% Housekeeping

% Send our data out for pre-handling
    
   [mal, fem, dat, tim] = wrenxc(malorig, femorig, Fs, td);

% This is the width of the clicking window in seconds. 
windwid = 2.5;

% This is the overlap in seconds between previous and next windows.
ovrlp = 0.2;

% We have to do this twice - once for female and once for male

for sex = 1:2
    
    if sex == 1; sng = fem; oth = mal; end
    if sex == 2; sng = mal; oth = fem; end

%% Initialize BSX and SYLNUM

% For the new analysis we start at time 0 (base X or bsx) and syllable #1.
bsx = 0; sylnum = 1; preclick = 0;

%% Loop through windows until user tells us that this is the end of sng
cntu = 1;

while cntu < 20
       
% Get clicks

        tt = find(tim >= bsx & tim < bsx + windwid);

        if preclick == 0 
            %figure(2); close(2); figure(2); % Don't be a hater!
            %figprop = get(gcf,'Position'); 
            %set(gcf,'Position',[figprop(1) figprop(2)-400 1000 400]);
            %oscsonB(sng(tt), oth(tt), Fs, [-100 50]); %plots Figure 2; 
            disp('Hit Enter to Move to next segment')
            cts = clickplotterB(sng(tt),oth(tt), Fs); % plots Figure 1;
        else
%             figure(2); 
%             figprop = get(gcf,'Position'); 
%             set(gcf,'Position',[figprop(1) figprop(2) 1000 400]);
%             oscson(oth(tt), Fs, [-100 50]); 
%             hold on; plot([ovrlp ovrlp], [100 5500], 'm'); hold off;
            disp('Hit Enter to Move to next segment')
            cts = clickplotterB(sng(tt), oth(tt),Fs, preclick);
        end
   
        if mod(length(cts),2) == 1 % TRY AGAIN LOSER
                shit = input('Odd numbers of clicks do not work, clickturd. ');
        end

% We loop to analyze each syllable

        for ss = 1:2:length(cts)-1

            % Get the region of data for the current click from the data
            stt = find(tim >= bsx + cts(ss) & tim < bsx + cts(ss+1));

            % Send the data to our analysis and tracing function
            temptemp = syldata(sng(stt),Fs);
 
                struct(sylnum).sylen = temptemp.sylen;
                struct(sylnum).loud = temptemp.loud;
                struct(sylnum).meanloud = temptemp.meanloud;
                struct(sylnum).spec_p = temptemp.spec_p;
                struct(sylnum).spec_f = temptemp.spec_f;
                struct(sylnum).spec_pfilt = temptemp.spec_pfilt;
                struct(sylnum).spec_peakf = temptemp.spec_peakf;
                struct(sylnum).spec_minf = temptemp.spec_minf;
                struct(sylnum).spec_maxf = temptemp.spec_maxf;
                struct(sylnum).spec_band = temptemp.spec_band;
                struct(sylnum).trace_freq = temptemp.trace_freq;
                struct(sylnum).trace_tim = temptemp.trace_tim;
                struct(sylnum).trace_slopemean = temptemp.trace_slopemean;
                struct(sylnum).trace_slopestd = temptemp.trace_slopestd;
                struct(sylnum).trace_slopevar = temptemp.trace_slopevar;
                struct(sylnum).trace_peakf = temptemp.trace_peakf;
                struct(sylnum).trace_peakt = temptemp.trace_peakt;
                struct(sylnum).trace_peakpt = temptemp.trace_peakpt;

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

%            if sex == 1; struct(sylnum).wavfile = nfo.Femwav; end;
%            if sex == 2; struct(sylnum).wavfile = nfo.Malwav; end;

            % PLOT ONTO OUR WINDOW
            figure(1); 
            subplot(2,1,1); hold on;
            timstr = struct(sylnum).trace_tim + struct(sylnum).syltim(1);
            plot(timstr - bsx, struct(sylnum).trace_freq, 'k','LineWidth',3);
            hold off;

            % Advance our syllable number now

            sylnum = sylnum + 1;
            
        end
                
        % We are done clicking through the window, so let's save temp data
        % and reset bsx
             
        bsx = struct(sylnum-1).syltim(2) - ovrlp;
        preclick = ovrlp;
        
        % Ask if we are done
        
        cntu = input('Done? 99=yes, 0 or [return] to continue: ');
        
        if isempty(cntu) == 1; cntu = 1; end
        cntu = cntu + 1;

end


% After all is said and done, get inter-syllable-intervals

    for j = 1:length(struct)-1
        struct(j).ISI = struct(j+1).syltim(1) - struct(j).syltim(2);
    end


    if sex == 1; fd = struct; clear struct; end
    if sex == 2; md = struct; clear struct; end
    
end

%%%%%%%%%%%%%%%%%%% END OF CLICKING

%% Post Clicking cleanup

if (length(fd) ~= length(md))
    shit = input('Number of Male and Female Syllables do not match, clickturd. ');
end

% Sex determination 

for j = length(fd):-1:1

    strt = min([fd(j).sylind(1) md(j).sylind(1)]);
    stp = max([fd(j).sylind(2) md(j).sylind(2)]);

    sx = length(find(dat(strt:stp) > 0)) - length(find(dat(strt:stp) < 0));

    if sx < 0; fd(j).sex = 0; md(j).sex = 0; end
    if sx > 0; fd(j).sex = 1; md(j).sex = 1; end
    
end
    
end          
   

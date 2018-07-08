function struct = fivesecrule(sng, Fs, old)
% fivesecrule
% This function depends on oscson for plotting.
% This function depends on syldata for analysis.
% This function depends on the SAMIII tree via syldata.
%
% out = fivesecrule(sng, Fs, [old]);
% sng is the sample data
% Fs is the sample rate in Hz
% old is the "in process" structure

%% Is this the first time, or are we continuing our analysis?

% Test to see if we are starting with a new file, or restarting an old file.
if nargin == 2; 

    newanal = 1;
% Get some basic information before we start
%    nfo.year = input('Enter the year, e.g. 2011: ');
%    nfo.month = input('Enter the month, e.g. 09: ');
%    nfo.day = input('Enter the day of the month, e.g. 07: ');
%    nfo.location{:} = input('Enter the location in single quotes: ');
%    nfo.note{:} = input('Enter a note in single quotes: ');
nfo.year=2011;nfo.month=12;nfo.day=25;nfo.location={'Yanayacu'}; nfo.note={'Fabulous'};

else

    newanal = 2;
    struct = old; % copy old data into our output structure

end;

%% Basic housekeeping before we start

% Make and time sequence and get the length of the signal.

tim = 1/Fs:1/Fs:length(sng)/Fs;
maxtim = tim(end);

% Initialize some variables.  yvals is for plotting verticle lines.
adv = 1; yvals = [0 5000]; 

% This is the width of the clicking window in seconds. 
windwid = 5;

% This is the overlap between the previous and the next window in seconds.
mvfd = windwid - 0.5;

%% Setup for clicking. If new, set variables. If old, get current status.

if newanal == 1; bsx = 0; sylnum = 1; end;

if newanal == 2; 
    
    bsx = old(end).bsx; 
    sylnum = length(old); 
    
    nfo.year = old(1).year; nfo.month = old(1).month; nfo.day = old(1).day;
    nfo.location = old(1).location; nfo.note = old(1).note;


end;

%% If old, let's plot the last (up to 15 seconds) of clicked data.

if newanal == 2;
    
    wnds = floor (bsx / 5);
        if wnds > 3; wnds = 3; end;
    
        
end;





%% Start Clicking

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% START OF CLICKING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The main loop of the program - we go until we've done the whole 
% length of time of the sample. It would be better if we could stop and
% pick up where we left off...

figure(1);

while bsx < tim(end);

	% get our indices for the time window to display
	tt = find(tim >= bsx & tim <= bsx+windwid);

	% adv is a variable that tells us to refresh the display... only 
	% necessary at the beginning and when we have a click that is in 
	% the last 0.5 seconds of the window.

	if adv == 1; 

		close (1); figure(1); % oscson(sng(tt),Fs); 

		figprop = get(gcf,'Position'); 
        	set(gcf,'Position',[figprop(1) figprop(2) 1000 400]);
		subplot(1,4,[1 2 3]);
		oscson(sng(tt),Fs);


	% Plot prior clicks in the main window - very helpful

		if sylnum > 1; 
            figure(1); hold on;
            if newanal == 1;
            plot([xclick(sylnum-1,2)-bsx xclick(sylnum-1,2)-bsx], yvals, 'm'); 
            end;
            if newanal == 2;
            plot([xclick(sylnum,2)-bsx xclick(sylnum,2)-bsx], yvals, 'm'); 
            end;
            hold off; 
		end;
	end;

	% Here we set adv to zero so that we don't advance... 
	adv = 0;

	% get two clicks for the start and end of the syllable
	clicksOK = 0;
    figure(1); subplot(1,4,[1 2 3]);
	[x y] = ginput(2);

    
% tssyl is the actual clicked duration of the syllable
% pssyl is with a 0.1 sec buffer on both sides for plotting  

tssyl = find (tim >= bsx+x(1) & tim < bsx+x(2));
pssyl = find (tim >= bsx+x(1)-0.1 & tim < bsx+x(2)+0.1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Try to plot what the user just clicked in a separate window with
    % 0.5 seconds on either side. We adjust the plot if the data is too 
    % close to the start or to the end of the signal.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure(1);

		if x(1)+bsx > 0.5 & x(2)+bsx < tim(end)-0.5; 
			samp = find (tim >= bsx+x(1)-0.1 & tim < bsx+x(2)+0.1);
	                subplot(1,4,4); 
			oscson(sng(samp), Fs);
			magentaline = length(samp)/Fs - 0.1; 
			hold on;
               plot([0.1 0.1], yvals, 'g', 'LineWidth',2);
               plot([magentaline magentaline], yvals, 'm', 'LineWidth',2);
			hold off;
		end;
		if x(1)+bsx <= 0.5;
			samp = find (tim < bsx+x(2)+0.1);
                        subplot(1,4,4); 
			oscson(sng(samp), Fs);
			magentaline = length(samp)/Fs - 0.1; 
			hold on;
               plot([x(1) x(1)], yvals, 'g', 'LineWidth',2);
               plot([magentaline magentaline], yvals, 'm', 'LineWidth',2);
			hold off;
		end;
        if x(1)+bsx > 0.5 & x(2)+bsx > tim(end)-0.1;
                 samp = find (tim >= bsx+x(1)-0.1);
                 subplot(1,4,4); 
                 oscson(sng(samp), Fs);
                 hold on;
                    plot([0.1 0.1], yvals, 'g', 'LineWidth',2);
                    magentaline = 0.5 + (x(2) - x(1));
                    plot([magentaline magentaline], yvals, 'm', 'LineWidth',2);
                 hold off;
        end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % The user gets the final say - ask whether the clicks are OK.
        	clicksOK = input('If bad type 0 or hit return; otherwise syl number or 99 to end: ');

% This is our signal, typing 99, to end the loop
        if clicksOK == 99; bsx = bsx + 99; 
        else; 
            
            
            if isempty(clicksOK) == 1; clicksOK = 0; end;

            if clicksOK == 0; honestly = input('Redo the click??? 1 yes 0 no ');
                if honestly == 0; clicksOK = input('Enter the syllable number: '); end;
            end;
            
	% If OK, then we add it to the final list
		if clicksOK ~= 0;

        xclick(sylnum,1) = x(1) + bsx; xclick(sylnum,2) = x(2) + bsx;

            
%%%%%%%%%%%%%%%%%%%%%
% Analyze and Trace
%%%%%%%%%%%%%%%%%%%%% 
 
 % execute complete analysis on the syllable chunk           

        anatrac = syldata(sng(tssyl),Fs);
 
 % add the missing components of our structure 

            anatrac.syl = [tim(tssyl(1)), tim(tssyl(end))];
            anatrac.ind = [tssyl(1) tssyl(end)];
            anatrac.category = 0;
            anatrac.bsx = bsx; 
            anatrac.xclick(1) = xclick(sylnum,1);
            anatrac.xclick(2) = xclick(sylnum,2);
            anatrac.year = nfo.year; anatrac.month = nfo.month;
            anatrac.day = nfo.day; anatrac.location = nfo.location;
            anatrac.note = nfo.note;
            struct(sylnum) = anatrac;           
            struct(sylnum).category = clicksOK;
            
            % save temporary.mat struct ; 

            if struct(sylnum).category < 8; fignum = 3; end;
            if struct(sylnum).category > 8; fignum = 4; end;

            figure(fignum); 
            if sylnum == 1; 
                subplot(8,4,1); 
                figprop = get(gcf,'Position'); 
                set(gcf,'Position',[figprop(1) figprop(2) 500 800]);
		axis 'off';
            end;
            
            subwndw = length(find([struct.category] == struct(sylnum).category)); 
            if subwndw > 4; a=floor(subwndw/4); subwndw = subwndw - a*4; end;
            subwndw = subwndw + ((struct(sylnum).category -1) * 4);
		if fignum == 4; subwndw = subwndw - 32; end;
            subplot(8,4,subwndw); 

            oscson(sng(pssyl),Fs); xlim([-0.2 0.5]);
		axis 'off';	
			if x(2) >= 2.5; 
				bsx = bsx + mvfd; 
				adv = 1;
			end;

			if x(2) < 2.5;
				figure(1); subplot(1,4,[1 2 3]); hold on; 
               			freqtim = anatrac.freqtim + x(1);
                % Plot the freq trace on top of the spectrogram in blue
				plot(freqtim, anatrac.freqtrace, 'b');
		% Plot the limits of the recent clicks
				plot([x(1) x(1)], yvals, 'g');
				plot([x(2) x(2)], yvals, 'm');
		% Plot the number of the syllable
				midsyl = x(1) + ((x(2) - x(1))/2);
				set(gcf,'DefaultTextColor','white')
				text(midsyl, 4000, num2str(clicksOK));
				set(gcf,'DefaultTextColor','black')
				hold off;
			end;

			clear x y anatrac;
			sylnum = sylnum + 1;
        end;
        
% This is our signal, typing 99, to end the loop
        if clicksOK == 99; bsx = bsx + 99; end;

        % close(2); % This was our review window. the collection windows will persist
        
end;
end;

%%%%%%%%%%%%%%%%%%% END OF CLICKING

%%

% might we be able to include these in 'struct' ???
%data.song = sng;
%data.Fs = Fs;
%data.time = tim;




% After all is said and done, get inter-syllable-intervals

	for j = 2:length(struct);
        	struct(j).ISI = struct(j).syl(1) - struct(j-1).syl(2);
    end;

%% the end
end

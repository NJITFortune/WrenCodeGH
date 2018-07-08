function struct = hagaclics(sng, Fs)
% hagaclics
% This function depends on oscson for plotting.
% This function depends on the SAMIII tree via syldata.m.
% out = hagaclics(sng, Fs);
%

% Get some basic information before we start
%struct.date.year = input('Enter the year, e.g. 2011: ');
%struct.date.month = input('Enter the month, e.g. 09: ');
%struct.date.day = input('Enter the day of the month, e.g. 07: ');
%struct.info.location{:} = input('Enter the location in single quotes: ');
%struct.info.note{:} = input('Enter a note in single quotes: ');

% Close all of the open figures.  That can be annoying! But it MUST be done...
% close all;

% Make and time sequence and get the length of the signal.
tim = 1/Fs:1/Fs:length(sng)/Fs;
maxtim = tim(end);

% Initialize some variables.  yvals is for plotting verticle lines.
adv = 1; bsx = 0; st = 1; yvals = [0 5000]; sylnum = 1;

% This is the width of the clicking window in seconds. 
windwid = 3;

% This is the overlap between the previous and the next window in seconds.
mvfd = windwid - 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% START OF CLICKING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The main loop of the program - we go until we've done the whole 
% length of time of the sample. It would be better if we could stop and
% pick up where we left off...

figure(30);

while bsx < tim(end);

	% get our indices for the time window to display
	tt = find(tim >= bsx & tim <= bsx+windwid);

	% adv is a variable that tells us to refresh the display... only 
    % necessary at the beginning and when we have a click that is in 
    % the last 0.5 seconds of the window.
	if adv == 1; 
		close (30); figure(30); oscson(sng(tt),Fs); 
		figprop = get(gcf,'Position'); 
        set(gcf,'Position',[figprop(1) figprop(2) 800 400]);

		if st > sylnum; 
            figure(30); hold on; 
            plot([xclick(st-1,2)-bsx xclick(st-1,2)-bsx], yvals, 'm'); 
            hold off; 
		end;
	end;

	% Here we set adv to zero so that we don't advance... 
	adv = 0;

	% get two clicks for the syllable
	clicksOK = 0;
    figure(30);
	[x y] = ginput(2);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Plot what the user just clicked in a separate window with
	% 0.5 seconds on either side.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
tssyl = find (tim >= bsx+x(1) & tim < bsx+x(2));
pssyl = find (tim >= bsx+x(1)-0.1 & tim < bsx+x(2)+0.1);

    figure(31);

		if x(1)+bsx > 0.5 & x(2)+bsx < tim(end)-0.5; 
			samp = find (tim >= bsx+x(1)-0.5 & tim < bsx+x(2)+0.5);
			oscson(sng(samp), Fs);
			magentaline = length(samp)/Fs - 0.5; 
			hold on;
               plot([0.5 0.5], yvals, 'g', 'LineWidth',2);
               plot([magentaline magentaline], yvals, 'm', 'LineWidth',2);
			hold off;
		end;
		if x(1)+bsx <= 0.5;
			samp = find (tim < bsx+x(2)+0.5);
			oscson(sng(samp), Fs);
			magentaline = length(samp)/Fs - 0.5; 
			hold on;
               plot([x(1) x(1)], yvals, 'g', 'LineWidth',2);
               plot([magentaline magentaline], yvals, 'm', 'LineWidth',2);
			hold off;
		end;
        if x(1)+bsx > 0.5 & x(2)+bsx > tim(end)-0.5;
                 samp = find (tim >= bsx+x(1)-0.5);
                 oscson(sng(samp), Fs);
                 hold on;
                    plot([0.5 0.5], yvals, 'g', 'LineWidth',2);
                    magentaline = 0.5 + (x(2) - x(1));
                    plot([magentaline magentaline], yvals, 'm', 'LineWidth',2);
                 hold off;
        end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % The user gets the final say - ask whether the clicks are OK.
        	clicksOK = input('If good type 0 or hit return; bad type 1; if OK and done with sample 99: ');
		
            if isempty(clicksOK) == 1; clicksOK = 0; end;

	% If OK, then we add it to the final list
		if clicksOK ~= 1;

            
        %%%%%%%%%%%%%%%%%%%%%
        % Analyze and Trace
        %%%%%%%%%%%%%%%%%%%%%
        % 
        % (sylnum)
        % 
            
        anatrac = syldata(sng(tssyl),Fs);

            % We will plot the trace on top of the spectrogram closeup
            figure(31); hold on;
            freqtim = anatrac.freqtim + 0.5;
            plot(freqtim, anatrac.freqtrace, 'g'); 
            
                        
            struct(sylnum).category = input('Please categorize this syllable #: ');
            pltwndw = struct(sylnum).category;  
            
            figure(pltwndw); 
               subwndw = length(find([struct.category] == struct(sylnum).category));
                if subwndw > 9; a=floor(subwndw/9); subwndw = subwndw - a*9; end;
               subplot(3,3,subwndw);
                oscson(sng(pssyl),Fs); xlim([-0.2 0.5]); 
                
            
%HUH% samp = find (tim >= bsx+x(1)-0.5 & tim < bsx+x(2)+0.5);
%HUH% 		oscson(sng(samp), Fs
            

			xclick(st,1) = x(1) + bsx; xclick(st,2) = x(2) + bsx;
		
			if x(2) >= 2.5; 
				bsx = bsx + mvfd; 
				adv = 1;
			end;
			if x(2) < 2.5;
				figure(30); hold on; 
				plot([x(1) x(1)], yvals, 'g');
				plot([x(2) x(2)], yvals, 'm');
				hold off;
			end;

			clear x y; 
			st = st + 1; % Obsolete?
			sylnum = sylnum + 1;
        end
        
% This is our signal, typing 99, to end the loop
        if clicksOK == 99; bsx = bsx + 99; end;

        close(31); % This was our review window. the collection windows will persist
        
end

%%%%%%%%%%%%%%%%%%% END OF CLICKING


struct.song = sng;
struct.Fs = Fs;
struct.time = tim;



%%%%%%%%%%%%%%%%%%% BEGIN Basic Syllable-by-Syllable Analyses

for i = 1:length(xclick);

	struct.syl(i,1) = xclick(i,1);
	struct.syl(i,2) = xclick(i,2);

% Perhaps the most useful bits - the durations and the indices

	struct.sylen(i) = diff(struct.syl(i,:));
	tmpa = find(tim >= xclick(i,1) & tim <= xclick(i,2));
	struct.ind(i,:) = [tmpa(1) tmpa(end)];

% get loudness - which is kinda ridiculous, but kinda not

	sy = struct.song(struct.ind(i,1):struct.ind(i,2));
	sy = sy - mean(sy);
	struct.loud(i) = sum(abs(sy));
	struct.meanloud(i) = struct.loud(i) / struct.sylen(i);

end

%%%%%%%%%%%%%%%%%%% END Basic Syllable-by-Syllable Analyses


% After all is said and done, get inter-syllable-intervals

	for j = 1:length(struct.syl)-1;
        	struct.ISI(j) = struct.syl(j+1,1) - struct.syl(j,2);
	end

% the end
end

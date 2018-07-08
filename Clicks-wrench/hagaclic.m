function struct = hagaclic(sng, Fs, old)
% hagaclic
% This function depends on oscson for plotting.
% This function depends on syldata for analysis.
% This function depends on the SAMIII tree via syldata.
%
% out = hagaclic(sng, Fs, [old]);
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
    
end;


%% OUR CLICKING WINDOW

figure(1);
figprop = get(gcf,'Position'); 
set(gcf,'Position',[figprop(1) figprop(2) 1000 400]);    

% get our indices for the time window to display
    tt = find(tim >= bsx & tim <= bsx+windwid);

    oscson(sng(tt),Fs); set(gcf,'HitTest','off');
    [x y] = getpts(gcf);

%%  REVIEW CLICKS ROUTINE

for i = 1:2:length(x);
        i
                close(1); figure(1);
	            figprop = get(gcf,'Position'); 
                set(gcf,'Position',[figprop(1) figprop(2) 1000 400]);
                subplot(1,4,[1 2 3]);
                oscson(sng(tt),Fs);
                hold on;
                plot([x(i) x(i)], yvals, 'g'); 
                plot([x(i+1) x(i+1)], yvals, 'm'); 
                hold off;
    
        %% Plot in the small window

		if x(i)+bsx > 0.5 & x(i+1)+bsx < tim(end)-0.5; 
			samp = find (tim >= bsx+x(i)-0.1 & tim < bsx+x(i+1)+0.1);
	                figure(1); subplot(1,4,4); 
			oscson(sng(samp), Fs);
			magentaline = length(samp)/Fs - 0.1; 
			hold on;
               plot([0.1 0.1], yvals, 'g', 'LineWidth',2);
               plot([magentaline magentaline], yvals, 'm', 'LineWidth',2);
			hold off;
		end;
		if x(i)+bsx <= 0.5;
			samp = find (tim < bsx+x(i+1)+0.1);
                        figure(1); subplot(1,4,4); 
			oscson(sng(samp), Fs);
			magentaline = length(samp)/Fs - 0.1; 
			hold on;
               plot([x(i) x(i)], yvals, 'g', 'LineWidth',2);
               plot([magentaline magentaline], yvals, 'm', 'LineWidth',2);
			hold off;
		end;
        if x(i)+bsx > 0.5 && x(i+1)+bsx > tim(end)-0.1;
                 samp = find (tim >= bsx+x(i)-0.1);
                        figure(1); subplot(1,4,4); 
                 oscson(sng(samp), Fs);
                 hold on;
                    plot([0.1 0.1], yvals, 'g', 'LineWidth',2);
                    magentaline = 0.5 + (x(i+1) - x(i));
                    plot([magentaline magentaline], yvals, 'm', 'LineWidth',2);
                 hold off;
        end;

clicksOK = input('If clicks OK, hit return or enter 0; FIX 1st type 1, 2nd type 2, both 3: ');
    if isempty(clicksOK) == 1; clicksOK = 0; end;
    
if clicksOK == 0;
    i
    % commit data here
    
end;

if clicksOK ~= 0;
    i
    figure(2);
    		if x(i)+bsx > 0.5 & x(i+1)+bsx < tim(end)-0.5; 
			samp = find (tim >= bsx+x(i)-0.1 & tim < bsx+x(i)+0.1);
			oscson(sng(samp), Fs);
			magentaline = length(samp)/Fs - 0.1; 
			hold on;
               plot([0.1 0.1], yvals, 'g', 'LineWidth',2);
               plot([magentaline magentaline], yvals, 'm', 'LineWidth',2);
			hold off;
		end;
		if x(i)+bsx <= 0.5;
			samp = find (tim < bsx+x(i+1)+0.1);
			oscson(sng(samp), Fs);
			magentaline = length(samp)/Fs - 0.1; 
			hold on;
               plot([x(i) x(i)], yvals, 'g', 'LineWidth',2);
               plot([magentaline magentaline], yvals, 'm', 'LineWidth',2);
			hold off;
		end;
        if x(i)+bsx > 0.5 & x(i+1)+bsx > tim(end)-0.1;
                 samp = find (tim >= bsx+x(i)-0.1);
                 oscson(sng(samp), Fs);
                 hold on;
                    plot([0.1 0.1], yvals, 'g', 'LineWidth',2);
                    magentaline = 0.5 + (x(i+1) - x(i));
                    plot([magentaline magentaline], yvals, 'm', 'LineWidth',2);
                 hold off;
        end;
    
    if clicksOK == 3; zz = 2; else; zz = 1; end;

    [xx yy] = ginput(zz);
        
end;

if clicksOK ~= 0; close(2); end;

end;




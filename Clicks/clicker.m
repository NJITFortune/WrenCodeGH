function [starts ends] = clicker(signal, Fs);
% Usage: [starts ends] = clicker(signal, Fs);
% This is not intended to be called by the user, but
% rather by other scripts.
% Click in pairs to get the starts and ends of syllable events.

clicking = 1; % Lets us know when to stop clicking in the while loop
yvals = [0 5000]; % Y-values for lines that we will draw on spectrogram

% Our loop - click until the user tells us they got it right...
while clicking ~= 0;

% If there is a mistake, we need to redo the clicking... that is
% what "contin" is for...
contin = 1;

% Make a figure with different dimensions - 800 wide by 400 tall
figure(1); figprop = get(gcf,'Position');
set(gcf,'Position',[figprop(1) figprop(2) 800 400]);

% This is for plotting wren signals. Perhaps we could alter this
% script here if we wanted to plot bat signals or other
oscson(signal, Fs);

% Get the clicks... continues until user hits "Return" button.
[x y] = ginput;

% Check that we have an even number of clicks - start/end pairs
	if mod(length(x),2) == 1;
		fprintf('You fool, you clicked an odd number of times.\n');
		clicking = 1;
		contin = 0;
	end

% If we have an even number of pairs, then plot the clicks for user review.
	if contin == 1;

		hold on;
	% Plot the starts (green) and ends (magenta).	
		for k = 2:2:length(x);
			plot([x(k-1) x(k-1)], yvals, 'g', 'LineWidth',2);
			plot([x(k) x(k)], yvals, 'm', 'LineWidth',2);
		end;
	% The user gets the final say
	clicking = input('If good type 0; bad type 1: ');

	end

% The thing can get funky if you don't close the figure
close(1);

end

% Our output...
	starts = x(1:2:end-1);
	ends = x(2:2:end);


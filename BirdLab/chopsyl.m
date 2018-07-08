function [syls nops] = chopsyl(signal,Fs, freq);

tim = 1/Fs : 1/Fs : length(signal)/Fs;

% Remove any DC offset
signal = signal - mean(signal);

% Get the envelope using Hilbert
envelope = abs(hilbert(signal));

% Low-pass filter the envelope
% FREQ IS KEY - a good starting value is 100
        order = 7; low = freq; 
        Wnlow = low*2/Fs;
        [bb,aa] = butter(order,Wnlow,'low');
        envelope = 100 * filtfilt(bb,aa,envelope);

while xx != 1
	% Plot the envelope and get a click
	figure(1); hold off;
	semilogy(envelope);
	[x,thresh] = ginput(1);

	hold on; plot([tim(1) tim(end)], [thresh thresh], 'm'); hold off;

	xx = input('Type 1 if OK, any other number if not'); 
end

% Now chop up the signal

above = find(envelope > thresh);

% Make a list of zeros that is the length of the signal.. 
zz = zeros(length(signal), 1);

% Now make every value that was above threshold equal to 1.
zz(above) = 1;

% To get the starts and ends we will use diff.
% 'diff' takes the difference between adjacent values...
% For example sample(2) - sample(1)...
yy = diff(zz);

% Start times will equal 1 and ends -1
starts = find( yy == 1 );
ends = find(yy == -1);

% Plot the starts and ends
figure(2); plot(tim,signal);
hold on;
length(starts)
for i = 1:length(starts)
plot([tim(starts(i)) tim(starts(i))], [-1 1],'g');
end
for j = 1:length(ends)
plot([tim(ends(j)) tim(ends(j))],[-1 1],'r');
end
hold off;

% Now this is simple!  We use a loop - the 'for' command
% to get each syllable.  
% syllable{1} will have the entire 1st syllable, and
% syllable{2} will have the enture 2nd syllable, etc.

for i = 1:length(starts)
  syls{i} = signal(starts(i):ends(i));
  timmy{i} = tim(starts(i):ends(i));
end

% Now we will get the silent parts between syllables...
% This can be important because sometimes the time
% between signals is an independent signal.  

% The procedure is almost the same as for the 
% syllables as above.  What is different is that instead
% of copying the signal from signal, we instead make
% pure silence by putting in a flat signal with an
% amplitude of 0.   

% for each of the ends (which is the start of each silence)
for j = 2:length(ends)

% Make a list of zeros with the length of the interval
% We get the length of the interval by getting the difference
% between the end of the silence, which is the 'start' of 
% the next syllable, and the start of the silence, which is
% the 'end' of the previous syllable.  Complicated!?!?

  nops{j-1} = zeros( (starts(j)-ends(j-1)), 1);

end


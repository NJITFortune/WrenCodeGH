function syls = slicer(in)
% syls = slicer(in);
% This is an interface for selecting syllables from frequency tracers
% created using "hagaclics" (via syldata).

%% Housekeeping

% Set some colors
% colrl(1,:)='r-'; colrl(2,:)='b-'; colrl(3,:)='m-'; colrl(4,:)='g-'; 
% colrl(5,:)='c-'; colrl(6,:)='k-'; colrl(7,:)='y-'; 
numcolors = 50;
    clrs = hsv(numcolors); % 50 colors for plotting up to 20 different fishes
%    oclrs = hsv(numcolors); % 20 colors for plotting up to 20 different fishes
%    clrs = zeros(numcolors, 3);
        % figure(3); clf; hold on; for kk=1:20; plot(kk,kk, '*', 'MarkerEdgeColor', oclrs(kk,:)); end;
%    shuff = [11 2 7 13 18 1 16 12 10 17 3 8 14 19 20 4 16 5 15 9];
%    for i=1:20; clrs(i,:) = oclrs(shuff(i),:); end; 


% This is the number of syllables that the user has passed us.
% ll = 1:length(in);

btch = 50; % How many we do per batch or...
sbtch = 1;  % Start with the first syllable

% Start with the first syllable, no selected syllables.
sylidx = 1;
slctd = 0; 

%% Loop to take care of all of the syllables - until all are selected
while slctd < length(in)

    % We can only handle so many at a time, which is determined by btch
    ll = sbtch:sbtch+btch;
    % This advance by btch could go beyond the total number of syllables - let's prevent this. 
    if ll(end) > length(in); ll=ll(1):length(in); end;


%% Cycle through until we get all of them in the current plot.
while length(ll) > 0
    
    % Plot the current 'll' traces - starts with all, but goes down each iteration
    figure(1); clf; 
    for i = 1:length(ll)
        plot(in(ll(i)).trace_tim, in(ll(i)).trace_freq); hold on;
    end;

% Get two clicks for a line (which is actually a box for the analysis)

figure(1); hold on;
[xx, yy] = ginputc(2,'Color','k','ShowPoints',true,'ConnectPoints',true);
figure(1); hold on; 
    plot(xx,yy,'r*-');
    xx = sort(xx); yy = sort(yy);
    pause(0.2);
    
% Determine which traces go into that line and 'save' them into syls.num

j = 1;

for i = 1:length(ll)
    % Find which ones are in the x region
    xWins = find(in(ll(i)).trace_tim > xx(1) & in(ll(i)).trace_tim < xx(2));

    if isempty(xWins) == 0 
        % Find which ones are also in the y region
        yWins = find(in(ll(i)).trace_freq(xWins) > yy(1) & in(ll(i)).trace_freq(xWins) < yy(2));
            if isempty(yWins) == 0
                % These are our winners!
                syls(sylidx).num(j) = ll(i); j=j+1;
            end;
    end;

    clear xWins yWins; 

end;

j=j-1;

% % % % Plot the selected traces in figure 2.

% % % if j >= 1  % Had better be unless the user is a stupid idiot
% % % 
% % %     for k = 1:length([syls(sylidx).num])
% % %         
% % %     figure(2); clf; hold on;
% % %     
% % %  %       colidx = sylidx;
% % %  %       if colidx > 7; colidx = mod(colidx,7); end;
% % %  %       if colidx == 0; colidx = 1; end;
% % %     
% % %     plot(in(syls(sylidx).num(k)).trace_tim,in(syls(sylidx).num(k)).trace_freq, '-', 'Color', clrs(sylidx,:));
% % %     end;
% % % 
% % % end;

% Now reset the number of items after removing the selected syllables
    ll = setdiff(ll,[syls(sylidx).num]);
    slctd = slctd + length([syls(sylidx).num]);
    sylidx = sylidx + 1;
    
end;
    
    sbtch = sbtch + btch + 1;

end;

%% Lets do a simple sort on each syllable type and grab the first of each

for j = length(syls):-1:1
    syls(j).num = sort(syls(j).num);
    tmpsyl(j) = syls(j).num(1);
end

% Now we sort that list
tmpsylist = sort(tmpsyl);

% And arrange according to first appearance.
for j = 1:length(tmpsyl)
    p = find(tmpsyl == tmpsylist(j));
    psyl(j).num = syls(p).num;
end;
syls = psyl;

%% Set the "id" value for the main structure tst
%
% for j=1:length(syls); for i=1:length(syls(j).num); tst(syls(j).num(i)).id = j; end; end;


%% Plots all of syllable type 2. Change k as necessary. Max 16 syllable reps
% (4x4 subplots)
%
% k=2; 
% for i=1:length([syls(k).num]);subplot(4,4,i);
% wrengram(a(tst(syls(k).num(i)).sylind(1):tst(syls(k).num(i)).sylind(2)),Fs); 
% end;

%% This plots all of the stick figures on a continuous plot... syls is
% syllable list from this script, tst is from hagaclics
%
%  colrl(1,:)='r-'; colrl(2,:)='b-'; colrl(3,:)='m-'; colrl(4,:)='g-'; 
%  colrl(5,:)='c-'; colrl(6,:)='k-'; colrl(7,:)='y-'; colrl(8,:)='r-';
%  colrl(9,:)='b-'; colrl(10,:)='m-'; colrl(11,:)='g-';
%
%for k=1:length(syls);
%    figure(3); hold on;
%    for i=1:length([syls(k).num]);
%        plot(tst(syls(k).num(i)).trace_tim + tst(syls(k).num(i)).syltim(1), tst(syls(k).num(i)).trace_freq,colrl(k));
%    end;
%end;

%figure(2);
%for j=1:min([16 length(syls)]);
%	subplot(4,4,j); hold on;
%        
%	for i = 1:length([syls(j).num]);
%        	plot(in(syls(j).num(i)).trace_tim,in(syls(j).num(i)).trace_freq,'k');
%            ylim([500 5500]); xx = in(syls(j).num(1)).trace_tim(end);
%            xlim([-0.05 xx+0.05]);
%	end;
%end;



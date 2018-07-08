function out = rcncl(dataM, slicedata)
% This takes the output from the "new" slicer 'foo'
% to allow you to combine syllables.

% This is the syllable

listola = 1:length(slicedata);

k=1; newgroups = 1;

while length(listola) > 0;

close all; % This is an annoying step, don't you think?

%% For the first syllable in the list "listola", plot all exemplars
for i=1:length([slicedata(listola(k)).num]); 
    
    % plot the first window in black - this is our selected syllable
	subplot(4,7,1); hold on; 
	plot(dataM(slicedata(listola(k)).num(i)).trace_tim,dataM(slicedata(listola(k)).num(i)).trace_freq,'k'); 
    
    % plot the subsequent windows in yellow - our selected syllable
	for j=2:28;
		subplot(4,7,j); hold on;
		plot(dataM(slicedata(listola(k)).num(i)).trace_tim,dataM(slicedata(listola(k)).num(i)).trace_freq,'g');

	end;

	maxlen(i) = dataM(slicedata(listola(k)).num(i)).trace_tim(end);
    
end;

    lenmax = max(maxlen);
    clear maxlen;
    
% Plot end lines at the start and end of the syllable
for j=1:28; subplot(4,7,j); 
    plot([0 0],[500 4500], 'm');
	xlim([-0.05 lenmax + 0.05]); ylim([500 4500]); 
	plot([lenmax lenmax],[500 4500], 'm');
end;

% Plot other syllables on top of that in red
for j=2:min([28 length(listola)]);
	subplot(4,7,j); hold on;
	for i = 1:length([slicedata(listola(j)).num]);
        	plot(dataM(slicedata(listola(j)).num(i)).trace_tim,dataM(slicedata(listola(j)).num(i)).trace_freq,'r');
	end;
    title(j,'FontSize',16);
end;

%% Now we manage the syllables
yy = input('Join?? ');

    if isempty(yy) == 1; 
        yy = 0; 
            tofu = [slicedata(listola(k)).num];
            for qq = 1:length(tofu); out(newgroups).num(qq) = tofu(qq); end;        
        newgroups = newgroups+1;
        listola = listola(2:end); 
    end;

    if yy(1) > 0; 
        yy = sort(yy);
        for ly = length(yy):-1:1;
        
            tofu = [slicedata(listola(k)).num slicedata(listola(yy(ly))).num]; 
            for qq = 1:length(tofu); out(newgroups).num(qq) = tofu(qq); end;
        
            listola = listola(find(listola ~= listola(yy(ly))));
        end;
        
        listola = listola(find(listola ~= listola(k)));
        newgroups = newgroups+1;
    end;            
        
    %listola = listola(find(listola ~= listola(k)))

end;

%% Now we want to sort - our slicing is somewhat random relative to the
% actual order of syllables.  So take the first syllable of each group
% and start with the first syllable and go from there.

for j = length(out):-1:1;
    out(j).num = sort(out(j).num);
    tmpsyl(j) = out(j).num(1);
end;
tmpsylist = sort(tmpsyl);
for j = 1:length(tmpsyl);
    p = find(tmpsyl == tmpsylist(j));
    psyl(j).num = out(p).num;
end;
out = psyl;


%% Now sex determine those suckers
%
%    for pp = 1:length(out);
%        % sexd( datam, dataf, sylnum, dat )
%        sex(pp) = sexd(dataM, dataF, out(pp).num, dF);
%        for qq = length([out(pp).num]):-1:1
%            out(pp).sex(qq) = sex(pp);
%        end;
%    end;


%% Review what we've done.
figure(2);
for j=1:min([16 length(out)]);
	subplot(4,4,j); hold on;
        tmp = length([out(j).num]);
        title(tmp,'FontSize',16);
	for i = 1:length([out(j).num]);
        	plot(dataM(out(j).num(i)).trace_tim,dataM(out(j).num(i)).trace_freq,'k');
            ylim([500 5500]); xx = dataM(out(j).num(1)).trace_tim(end); xlim([-0.05 xx+0.05]);
	end;
end;



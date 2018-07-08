function out = wrencleaner(in, dat, Fs)
% out = wrencleaner(in)
% Where in is the structure from hagaclics and out is the "cleaned up"
% version. This focusses on fixing the most common errors in
% the syllable trace.

out = in;

%% Loop through each syllable and check and correct the most common errors
for i = 1:length(in);
    
%% Check for zeros - we should never have a zero frequency

    if ismember(0,[in(i).trace_freq]) == 1;

        % Get the indices where these zeros happened
        zips = find([in(i).trace_freq] == 0);

        % Most of the time there is only one zero value, usually at the
        % end - so we can easily fix that.

        if length(zips) == 1; 
            % If the zero is anywhere but the first value, we'll simply
            % copy the previous value
            if zips ~= 1; 
                out(i).trace_freq(zips) = in(i).trace_freq(zips-1);
            else
                out(i).trace_freq(zips) = in(i).trace_freq(zips+1);
            end;
        else
            % But if we have more than a single zero value, we'll need
            % to get user to help. Here we have the user click each zero.
            sprintf('Multiple Zeros');
            NumberofUserClicks = length(zips);
            
            figure(1); 
            hold off;
            wrengram(dat(in(i).sylind(1):in(i).sylind(2)), Fs);
            hold on; 
            plot(in(i).trace_tim, in(i).trace_freq,'*-');
            plot(in(i).trace_tim(zips), in(i).trace_freq(zips), 'ro');
            
            clks = ginput(NumberofUserClicks);
            out(i).trace_freq(zips) = clks(:,2);            
        end;    

    end;    

%% Check for changes in frequency that are too big to be real
    
    % Find FM slopes of greater than 500Hz/sample.
    spks = find(abs(diff(out(i).trace_freq)) > 500);

    % If we find a large FM event, we need to fix it.
    if (isempty(spks) == 0);

        % Fill in between areas. We ask the user to click and we go
        % forwards and backwards from there.
            
        figure(2); hold off; 
        plot(in(i).trace_tim, out(i).trace_freq,'b-*');
        hold on;
        plot(in(i).trace_tim(spks), out(i).trace_freq(spks), 'r-*');

        bifur = ginput(1);
        bifur = bifur(1);
        
        % Do the right side of the click.
        % This code sucks - probably good enough, but not a good general
        % strategy.
        gtbifur = find (in(i).trace_tim(spks) > bifur);
        if (isempty(gtbifur) == 0);
            if (length(gtbifur) == 1); 
                out(i).trace_freq(spks(gtbifur):end) = out(i).trace_freq(spks(gtbifur)-1);
            end;
            if (length(gtbifur) == 2);
                out(i).trace_freq(spks(gtbifur(1)):spks(gtbifur(2))) = ...
                    mean([out(i).trace_freq(spks(gtbifur(1)-1)), out(i).trace_freq(spks(gtbifur(2))+1)]);
            end;    
            if (length(gtbifur) == 3);
                out(i).trace_freq(spks(gtbifur(1)):spks(gtbifur(2))) = ...
                    mean([out(i).trace_freq(spks(gtbifur(1))), out(i).trace_freq(spks(gtbifur(2))+1)]);
                out(i).trace_freq(spks(gtbifur(3)):end) = out(i).trace_freq(spks(gtbifur(3))-1);
            end;            
        end;
        
        % Do the left side of the click.
        % This code sucks - probably good enough, but not a good general
        % strategy.
        ltbifur = find (in(i).trace_tim(spks) < bifur);
        if (isempty(ltbifur) == 0);
            if (length(ltbifur) == 1); 
                out(i).trace_freq(1:spks(ltbifur)) = out(i).trace_freq(spks(ltbifur)+1);
            end;
            if (length(ltbifur) == 2);
                out(i).trace_freq(spks(ltbifur(1)):spks(ltbifur(2))) = ...
                    mean([out(i).trace_freq(spks(ltbifur(1))), out(i).trace_freq(spks(ltbifur(2))+1)]);
            end;    
            if (length(ltbifur) == 3);
                out(i).trace_freq(spks(ltbifur(2)):spks(ltbifur(3))) = ...
                    mean([out(i).trace_freq(spks(ltbifur(2))), out(i).trace_freq(spks(ltbifur(3))+1)]);
                out(i).trace_freq(spks(1:ltbifur(1))) = out(i).trace_freq(spks(ltbifur(1))+1);
            end;            
        end;
        clear gtbifur ltbifur bifur
    end;    
end;


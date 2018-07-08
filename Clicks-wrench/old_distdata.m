function anal = distdata(foo, out)
% Usage: function anal = distdata(foo, out)
% Where foo and out are from distance clicking
% Depends on the scripts 

%% Setup

% Colors for plotting
clrs(1,:)='k*'; clrs(2,:)='r*'; clrs(3,:)='c*'; clrs(4,:)='g*'; clrs(5,:)='y*'; 

% Get the sample rate of the signal
Fs = 1 / (out.ftim(2) - out.ftim(1)); 

% Sanity check - do both animals have the same number of syllables?
if length(foo.m) ~= length(foo.f)
    fprintf('Male and Female syllable counts do not match \n');
    % CALL SCRIPT FOR SYLLABLE NUMBER RESOLUTION
end;

% Show the user the clicked data - female on top
fprintf('Plot data and label. \n');
figure(1); clf;

% Specgrams
    xa(1) = subplot(211); specgram(out.fc.fem, 2048, Fs, [], 2000); caxis([0 30]);
    xa(2) = subplot(212); specgram(out.mc.mal, 2048, Fs, [], 2000); caxis([0 30]);
    
    linkaxes(xa, 'xy'); ylim([200 4500]); colormap('HOT');     

% Plot the starts and ends on the female specgram
    figure(1); subplot(211); hold on;
        for i=1:length(foo.f)
            plot([foo.f(i).syltim(1) foo.f(i).syltim(1)], [500 5000], 'g');
            plot([foo.f(i).syltim(2) foo.f(i).syltim(2)], [500 5000], 'r');
            text(foo.f(i).syltim(1) + 0.010, 4000, num2str(i), 'Color', 'w');
            famp(i) = sum(foo.f(i).trace_amp);
        end

% Plot the starts and ends on the male specgram
    figure(1); subplot(212); hold on;
        for i=1:length(foo.m)
           plot([foo.m(i).syltim(1) foo.m(i).syltim(1)], [500 5000], 'g');
           plot([foo.m(i).syltim(2) foo.m(i).syltim(2)], [500 5000], 'r');
           text(foo.m(i).syltim(1) + 0.010, 4000, num2str(i), 'Color', 'w');
           mamp(i) = sum(foo.m(i).trace_amp);
        end

%% Use amplitude to figure out who is male and who is female.

% These will be the list of syllable indices for male and for female
malsylidx = []; femsylidx = [];

% We will mark the male and female syllables as we go along, so call the figure
figure(1); 

% This loop goes through each syllable in the female structure (could have
% used male structure if we wanted - we already checked to make sure that
% both were the same length)
    for i=1:length(foo.f) 
        if famp(i) - mamp(i) > 0 % Female signal louder than male signal
            femsylidx(end+1) = i; % Then this is a female syllable
                % Mark both subplots
                subplot(211); text(foo.f(i).syltim(1) + 0.010, 4000, num2str(i), 'Color', 'm');
                subplot(212); text(foo.m(i).syltim(1) + 0.010, 4000, num2str(i), 'Color', 'm');
        end;
        if famp(i) - mamp(i) < 0 % Male signal louder than female signal
            malsylidx(end+1) = i; % Then this is a male sylable
                % Mark both subplots
                subplot(211); text(foo.f(i).syltim(1) + 0.010, 4000, num2str(i), 'Color', 'b');
                subplot(212); text(foo.m(i).syltim(1) + 0.010, 4000, num2str(i), 'Color', 'b');
        end;
    end;

%% Use kmeans to segregate syllables produced by the female and male birds.

% At the moment we are using only the syllable lengths (sylen) and peak frequencies (trace_peakf).

% Using a temporary variable ff as an array of datums for kmeans
    ff = [];
    ff(:,1) = [foo.f(femsylidx).sylen] * 10000;
    ff(:,2) = [foo.f(femsylidx).trace_peakf];

% Using a temporary variable mm as an array of datums for kmeans
    mm = [];
    mm(:,1) = [foo.m(malsylidx).sylen] * 10000;
    mm(:,2) = [foo.m(malsylidx).trace_peakf];

% Ask the user how many clusters of female and male syllables to make.
figure(2); clf; % Ploting the syllable length (sylen) against peak frequency (trace_peakf).
    subplot(121); plot([foo.f(femsylidx).sylen], [foo.f(femsylidx).trace_peakf], 'm*'); ylim([500 6000]);
    subplot(122); plot([foo.m(malsylidx).sylen], [foo.m(malsylidx).trace_peakf], 'b*'); ylim([500 6000]);

    numFclusts = input('How many female syllable types? ');
    numMclusts = input('How many male syllable types? ');

% ASK ABOUT COMPLETE HAND LABELLING HERE IN CASE TOO FEW EXAMPLES FOR Kmeans %%%%%%%%

fprintf('Using k-means to segregate syllable types. \n');
% Kmeans for female
femtypes = kmeans(ff, numFclusts, 'Distance', 'cityblock', 'Replicates', 1000);
    for kk = 1:numFclusts
        femtype{kk} = femsylidx(femtypes == kk);
            if length(femtype{kk}) > 1
                femintervals{kk} = diff(femtype{kk});
            end;
    end;

% Kmeans for male
maltypes = kmeans(mm, numMclusts, 'Distance', 'cityblock', 'Replicates', 1000);
    for kk = 1:numMclusts
        maltype{kk} = malsylidx(maltypes == kk);
            if length(maltype{kk}) > 1
                malintervals{kk} = diff(maltype{kk});
            end;
    end;

% Replot the data but with colors for each kmeans group
    figure(2); clf;

    subplot(121); hold on; ylim([500 5000]);
        for i = 1:length(femtype)
            plot([foo.f(femtype{i}).sylen], [foo.f(femtype{i}).trace_peakf], clrs(i,:)); ylim([500 6000]);
        end;

    subplot(122); hold on; ylim([500 5000]);
        for i = 1:length(maltype)
            plot([foo.m(maltype{i}).sylen], [foo.m(maltype{i}).trace_peakf], clrs(i,:)); ylim([500 6000]);
        end;

% Plot syllable traces on top of each other.

% These are the concatenated difference in indices for use with calculating the mode
% intervals between repetitions of a particular syllable type.
femintcon = []; malintcon = []; 

figure(4); clf; % Female syllables
    for jj = 1:length(femtype)
        femints{jj} = diff([femtype{jj}]);
        femintcon = [femintcon femints{jj}];
        for kk = 1:length(femtype{jj})
            subplot(1,length(femtype),jj); hold on; ylim([200 5000]);
            plot(foo.f(femtype{jj}(kk)).trace_tim, foo.f(femtype{jj}(kk)).trace_freq, 'm');
        end
    end
    
figure(5); clf; % Male syllables
    for jj = 1:length(maltype)
        malints{jj} = diff([maltype{jj}]);
        malintcon = [malintcon malints{jj}];
        for kk = 1:length(maltype{jj})
            subplot(1,length(maltype),jj); hold on; ylim([200 5000]);
            plot(foo.m(maltype{jj}(kk)).trace_tim, foo.m(maltype{jj}(kk)).trace_freq, 'b');
        end
    end;

%% Guess at intervals that are wrong to identify misclassified syllables.  

% Most of the time all of the intervals should be the same as the mode

    modint = mode([malintcon femintcon]); % These are the concatenated intervals from above

%% FEMALE SYLLABLES FIX
fprintf('Testing female syllable intervals. \n');

    wrongshortfem = {}; wronglongfem = {}; % These will have the indices of the starting syllable for bad intervals

for pp = 1:length(femints) % Cycle through each syllable type for female
    wronglongfem{pp} = femtype{pp}([femints{pp}] > modint); % These have only long bad intervals

    if ~isempty(wronglongfem{pp}) 
        fprintf('Found %i long intervals in female syllable %i . \n', length(wronglongfem{pp}), pp);
        fprintf(' The syllables are: '); femtype{pp}
        fprintf('\n The intervals are: '); femints{pp}    
    end 
%     if ~isempty(wrongshortfem{pp}) 
%         fprintf('Found %i short intervals in syllable %i. \n', length(wrongshortfem{pp}), pp); 
%         fprintf(' The syllables are: '); femtype{pp}
%         fprintf('\n The intervals are: '); femints{pp}
%         
%     end
end

% Let's try and fix syllables

% Did we find any bad intervals? Let's start with the longs...
if ~isempty(wronglongfem) % Then there is a bad long interval in female data
    fprintf('found a wronglongfem. \n'); 
    for jj = 1:length(wronglongfem) 
        if ~isempty(wronglongfem{jj}) %identifies probelm syllable type-
            fprintf('found a problem in syllable type %i. \n', jj);
            for kk = 1:length(wronglongfem{jj}) 
               if femints{jj}(femints{jj} ~= modint) == modint * 2 %finds double interval, the one we need to fix-
                fprintf('We have a winner - double interval. \n');
                transfersyl = (wronglongfem{jj}(kk) + modint); 
                
                % Eliminate the syllable from incorrect syllable type
                
                    for qq = 1:length(femtype)
                        femtype{qq} = femtype{qq}(find(femtype{qq} ~= transfersyl));
                    end
                    
                % Insert the syllables into the correct syllable type and make ascending order
                    femtype{jj} = sort([femtype{jj} transfersyl]); 
               end    
            end
        end
    end
end

%% MALE SYLLABLES FIX
fprintf('Testing male syllable intervals. \n');

    wrongshortmal = {}; wronglongmal = {}; % These will have the indices of the starting syllable for bad intervals

for pp = 1:length(malints) % Cycle through each syllable type for female
    wronglongmal{pp} = maltype{pp}([malints{pp}] > modint); % These have only long bad intervals

    if ~isempty(wronglongmal{pp}) 
        fprintf('Found %i long intervals in male syllable %i . \n', length(wronglongmal{pp}), pp);
        fprintf(' The syllables are: '); maltype{pp} 
        fprintf('\n The intervals are: '); malints{pp}        
    end 
end

if ~isempty(wronglongmal) % Then there is a bad long interval in male data
    fprintf('found a wronglongmal. \n'); 
    for jj = 1:length(wronglongmal) 
        if ~isempty(wronglongmal{jj}) %identifies probelm syllable type-
            fprintf('found a problem in syllable type %i. \n', jj);
            for kk = 1:length(wronglongmal{jj}) 
               if malints{jj}(malints{jj} ~= modint) == modint * 2 %finds double interval, the one we need to fix-
                fprintf('We have a winner - double interval in MALE. \n');
                transfersyl = (wronglongmal{jj}(kk) + modint); 
                
                % Eliminate the syllable from incorrect syllable type
                
                    for qq = 1:length(maltype)
                        maltype{qq} = maltype{qq}(maltype{qq} ~= transfersyl);
                    end
                    
                % Insert the syllables into the correct syllable type and make ascending order
                    maltype{jj} = sort([maltype{jj} transfersyl]); 
                
               end    
            end
        end
    end
end

%% plot and cleanup

figure(3); clf;
    subplot(121); hold on; ylim([500 6000]);
        for i = 1:length(femtype)
            plot([foo.f(femtype{i}).sylen], [foo.f(femtype{i}).trace_peakf], clrs(i,:)); ylim([500 6000]);
        end;

    subplot(122); hold on; ylim([500 6000]);
        for i = 1:length(maltype)
            plot([foo.m(maltype{i}).sylen], [foo.m(maltype{i}).trace_peakf], clrs(i,:)); ylim([500 6000]);
        end;


anal.femsyls = femsylidx;
%anal.femnum =
%anal.fem

anal.malsyls = malsylidx;
%anal.malnum = 


%% Award-Winning Code
%    wrongshortfem{pp} = femtype{pp}([femints{pp}] < modint); % These have only bad short intervals

%     figure(6); clf;
% 
% if length(wrongshortfem) > 1
%     for i = 1:length(wrongshortfem)
%         
%         if ~isempty(wrongshortfem{i} > 1)
%             
%         startsyl = min(wrongshortfem{i}) - 1; 
%             if startsyl < 1; startsyl = 1; end;
%         endsyl = max(wrongshortfem{i}) + 1; 
%             if endsyl > femtype{i}(end); endsyl = femtype{i}(end); end;
% 
%         for jj = startsyl:endsyl
%             subplot(length(wrongshortfem), 1+endsyl-startsyl, jj - startsyl + 1);
%             ylim([500 5000]);
%             plot(foo.f(jj).trace_tim, foo.f(jj).trace_freq, 'm');
%         end
%         end
%     end
% end
% 
% figure(7); clf;
%  if length(wronglongfem) > 1
%     for i = 1:length(wronglongfem)
%         
%         if ~isempty(wronglongfem{i} > 1)
%             
%         startsyl = min(wronglongfem{i}) - 1; 
%             if startsyl < 1; startsyl = 1; end;
%         endsyl = max(wronglongfem{i}) + 1; 
%             if endsyl > femtype{i}(end); endsyl = femtype{i}(end); end;
% 
%         for jj = startsyl:endsyl
%             subplot(length(wronglongfem), 1+endsyl-startsyl, jj - startsyl + 1);
%             ylim([500 5000]);
%             plot(foo.f(jj).trace_tim, foo.f(jj).trace_freq, 'm');
%         end
%         end
%     end
% end


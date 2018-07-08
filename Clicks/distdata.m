function anal = distdata(foo, out)
% Usage: function anal = distdata(foo, out)
% Where foo and out are from distance clicking
% Depends on the scripts 
% anal.mamp / famp are amplitudes of syllables at the male microphone / female microphone by index.
% anal.femtype / maltype are cell arrays, each array has the syllable indices for that syllable type.

%% Setup

% Sanity check - do both animals have the same number of syllables?
if length(foo.m) ~= length(foo.f)
    fprintf('Male and Female syllable counts do not match! \n');
    fprintf('Please reclick this entry and try again. \n');
end

% Colors and plotting styles
clrs(1,:)='k*'; clrs(2,:)='r*'; clrs(3,:)='c*'; clrs(4,:)='g*'; clrs(5,:)='m*'; 

% Get the sample rate of the signal
    Fs = 1 / (out.ftim(2) - out.ftim(1)); 

% Show the user the clicked data - female on top 
    fprintf('Plot data and label. \n');
    figure(1); clf;

% Specgrams
    xa(1) = subplot(211); specgram(out.fc.fem, 2048, Fs, [], 2000); caxis([-15 30]);
    xa(2) = subplot(212); specgram(out.mc.mal/2, 2048, Fs, [], 2000); caxis([-15 30]); % ADDED /2
    
    linkaxes(xa, 'xy'); ylim([200 4500]); colormap('HOT');     

% Plot the starts and ends on the female specgram and calculate the amplitude famp
    figure(1); subplot(211); hold on;
        for i=1:length(foo.f)
            plot([foo.f(i).syltim(1) foo.f(i).syltim(1)], [500 5000], 'g');
            plot([foo.f(i).syltim(2) foo.f(i).syltim(2)], [500 5000], 'r');
            text(foo.f(i).syltim(1) + 0.010, 4000, num2str(i), 'Color', 'w');
            famp(i) = sum(foo.f(i).trace_amp);
        end

% Plot the starts and ends on the male specgram and calculate the amplitude mamp
    figure(1); subplot(212); hold on;
        for i=1:length(foo.m)
           plot([foo.m(i).syltim(1) foo.m(i).syltim(1)], [500 5000], 'g');
           plot([foo.m(i).syltim(2) foo.m(i).syltim(2)], [500 5000], 'r');
           text(foo.m(i).syltim(1) + 0.010, 4000, num2str(i), 'Color', 'w');
           mamp(i) = sum(foo.m(i).trace_amp);
        end

        
%% Use amplitude to figure out who is male and who is female.

% Normalize the amplitudes at each microphone (microphones can be set to radically 
% different gains... and the placement of microphones also affects the raw amplitudes of signals.
    
    anal.mamp = mamp / max(mamp);
    anal.famp = famp / max(famp);
        mamp = anal.mamp; famp = anal.famp;

% These will be the list of syllable indices for male and for female
    malsylidx = []; femsylidx = [];

% We will mark the male and female syllables as we go along, so call the figure
    figure(1); 

% This loop goes through each syllable in the female structure (could have
% used male structure if we wanted - we already checked to make sure that
% both were the same length)

    for i=1:length(foo.f) % For each syllable...
        
        if famp(i) - mamp(i) > 0 % Female signal louder than male signal
            femsylidx(end+1) = i; % Then this is a female syllable
                % Label both subplots with syllable numbers
                subplot(211); text(foo.f(i).syltim(1) + 0.010, 4000, num2str(i), 'Color', 'm'); 
                subplot(212); text(foo.m(i).syltim(1) + 0.010, 4000, num2str(i), 'Color', 'm');
        end;
        
        if famp(i) - mamp(i) < 0 % Male signal louder than female signal
            malsylidx(end+1) = i; % Then this is a male sylable
                % Label both subplots with syllable numbers
                subplot(211); text(foo.f(i).syltim(1) + 0.010, 4000, num2str(i), 'Color', 'b');
                subplot(212); text(foo.m(i).syltim(1) + 0.010, 4000, num2str(i), 'Color', 'b');
        end;
    end;

    % Recheck low amplitude syllables because the automatic sex determination is often wrong % figure(27); clf; hold on; plot(famp, '-*m'); plot(mamp, '-*b'); % Show the relative amplitudes of the syllables    

    for jj = 1:length(famp)
        lowampthresh = 0.3; % This threshold for low amplitude syllables may need to be reconsidered
        if (anal.famp(jj) <  lowampthresh && anal.mamp(jj) < lowampthresh)
           % Then we have a low amplitude syllable and the sex determination can be innaccurate
           
           fprintf('Syllable %i is very low amplitude. \n', jj);
           mf = input('Is it Male (1) or Female (2)?: ');
           
           currsex = 2*length(find(femsylidx == jj)) + length(find(malsylidx == jj)); % 2 if female, 1 if male
           
           if mf > currsex % We need to change from male to female
                malsylidx = malsylidx(malsylidx ~= jj);
                femsylidx(end+1) = jj; femsylidx = sort(femsylidx);
           end
           if mf < currsex % We need to change from female to male
                femsylidx = femsylidx(femsylidx ~= jj);
                malsylidx(end+1) = jj; malsylidx = sort(malsylidx);               
           end
           
        end
    end
    
    % ADD PLOTTING OF FIXED SYLLABLES
    
    % Ask if there are mis-classified syllables and fix them
 
    fprintf('Sex classification OK, or change syllable? If ok, enter 99.\n');
    gonogo = input('If not OK, give syllable number you want to change: ');

    while gonogo < 90
    
        currsex = 2*length(find(femsylidx == gonogo)) + length(find(malsylidx == gonogo)); % 2 if female, 1 if male
        mf = input('Is it Male (1) or Female (2)?: ');

           if mf > currsex % We need to change from male to female
                malsylidx = malsylidx(malsylidx ~= gonogo);
                femsylidx(end+1) = gonogo; femsylidx = sort(femsylidx);
           end
           if mf < currsex % We need to change from female to male
                femsylidx = femsylidx(femsylidx ~= gonogo);
                malsylidx(end+1) = gonogo; malsylidx = sort(malsylidx);               
           end
    gonogo = input('Done?  (99) Or another change? (syllable number): ');
    
    end
    
%% Segregation 1: For duets with fewer than 16 syllables, we will do hand segregation.    
    
    if length(foo.m) < 16 
        
        femtype = {}; maltype = {};
        
    fprintf('This is a short duet, so you get to do hand sorting! \n');
% A figure with all of the syllable traces
    figure(2); clf; xxa(1) = subplot(211); hold on; xxa(2) = subplot(212); hold on;
        for jj=1:length(foo.m)
            subplot(211); plot(foo.f(jj).trace_tim + foo.f(jj).syltim(1), foo.f(jj).trace_freq, 'm');
                plot([foo.f(jj).syltim(1) foo.f(jj).syltim(1)], [1000 4000], 'g'); 
            subplot(212); plot(foo.m(jj).trace_tim + foo.m(jj).syltim(1), foo.m(jj).trace_freq, 'b');
                plot([foo.m(jj).syltim(1) foo.m(jj).syltim(1)], [1000 4000], 'g'); 
            if find(femsylidx == jj)
                subplot(211); text(foo.f(jj).syltim(1) + 0.050, 750, num2str(jj), 'Color', 'm');
                subplot(212); text(foo.m(jj).syltim(1) + 0.050, 750, num2str(jj), 'Color', 'm');
            end
            if find(malsylidx == jj)
                subplot(211); text(foo.f(jj).syltim(1) + 0.050, 750, num2str(jj), 'Color', 'b');
                subplot(212); text(foo.m(jj).syltim(1) + 0.050, 750, num2str(jj), 'Color', 'b');
            end            
        end
        linkaxes(xxa, 'x');
        
        femlistofcategories = []; mallistofcategories = []; tmpfemsylist = []; tmpmalsylist = [];


        for jj = 0:5:length(foo.m) % Removed -5
            
            figure(3); clf; qwe(1) = subplot(211); hold on; qwe(2) = subplot(212); hold on; 

 
            for kk = 1:min([5 length(foo.m)-jj])
                if length(find(femsylidx == jj+kk)) == 1 % Female syllable
                subplot(211); plot(foo.f(jj+kk).trace_tim + foo.f(jj+kk).syltim(1), foo.f(jj+kk).trace_freq, 'm-', 'LineWidth', 2);
                tmpfemsylist = [tmpfemsylist jj+kk];
                end
                if length(find(malsylidx == jj+kk)) == 1 % Male syllable                
                subplot(212); plot(foo.m(jj+kk).trace_tim + foo.m(jj+kk).syltim(1), foo.m(jj+kk).trace_freq, 'b-', 'LineWidth', 2);
                tmpmalsylist = [tmpmalsylist jj+kk];
                end
            end
            linkaxes(qwe, 'x');
            
                asdf = input('Enter FEMALE syllables categories as a list e.g. [1 2 3 1 2]: ');
                    femlistofcategories = [femlistofcategories asdf];
                asdf = input('Enter MALE syllables categories as a list e.g. [1 2 1 2 1]: ');
                    mallistofcategories = [mallistofcategories asdf];
        
        end
            
        % Build the syllable type structures
        % Female
            for pp = 1:max(femlistofcategories)
                femtype{pp} = tmpfemsylist(femlistofcategories == pp);
            end
            % Female
            for pp = 1:max(mallistofcategories)
                maltype{pp} = tmpmalsylist(mallistofcategories == pp);
            end

        
        
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
    end

    end % End of short duet analysis
    
%% Segregation 2: For duets with 16 or more syllables, we will try and use kmeans to segregate syllables.
    if length(foo.m) >= 16 
% At the moment we are using only the syllable lengths (sylen) and peak frequencies (trace_peakf).

% Using a temporary variable ff as an array of datums for kmeans for female
% syllables
    ff = [];
    ff(:,1) = [foo.f(femsylidx).sylen] * 10000;
    ff(:,2) = [foo.f(femsylidx).trace_peakf];

% Using a temporary variable mm as an array of datums for kmeans for male
% syllables
    mm = [];
    mm(:,1) = [foo.m(malsylidx).sylen] * 10000;
    mm(:,2) = [foo.m(malsylidx).trace_peakf];

% Ask the user how many clusters of female and male syllables to make.
figure(5); clf; % Ploting the syllable length (sylen) against peak frequency (trace_peakf).
    subplot(121); plot([foo.f(femsylidx).sylen], [foo.f(femsylidx).trace_peakf], 'm*'); ylim([500 5000]);
    subplot(122); plot([foo.m(malsylidx).sylen], [foo.m(malsylidx).trace_peakf], 'b*'); ylim([500 5000]);

    numFclusts = input('How many female syllable types? ');
    numMclusts = input('How many male syllable types? ');

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

    subplot(221); hold on; ylim([500 5000]);
        for i = 1:length(femtype)
            plot([foo.f(femtype{i}).sylen], [foo.f(femtype{i}).trace_peakf], clrs(i,:)); ylim([500 5000]);
        end;

    subplot(222); hold on; ylim([500 5000]);
        for i = 1:length(maltype)
            plot([foo.m(maltype{i}).sylen], [foo.m(maltype{i}).trace_peakf], clrs(i,:)); ylim([500 5000]);
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

    
    end
    
    
%% Guess at intervals that are wrong to identify misclassified syllables.  

% Most of the time all of the intervals should be the same as the mode

    modint = mode([malintcon femintcon]); % These are the concatenated intervals from above

fprintf('Testing female syllable intervals. \n');
    somethingwentwrongFEMALE = [];
    wronglongfem = {}; % These will have the indices of the starting syllable for bad intervals

for pp = 1:length(femints) % Cycle through each syllable type for female
    wronglongfem{pp} = femtype{pp}([femints{pp}] > modint); % These have only long bad intervals

    if ~isempty(wronglongfem{pp}) 
        fprintf('Found %i long intervals in female syllable %i . \n', length(wronglongfem{pp}), pp);
        fprintf(' The syllables are: '); femtype{pp}
        fprintf('\n The intervals are: '); femints{pp}    
        somethingwentwrongFEMALE = 1;
    end 
%     if ~isempty(wrongshortfem{pp}) 
%         fprintf('Found %i short intervals in syllable %i. \n', length(wrongshortfem{pp}), pp); 
%         fprintf(' The syllables are: '); femtype{pp}
%         fprintf('\n The intervals are: '); femints{pp}
%         
%     end
end

% Let's try and fix syllables-code below to be replaced with 

% Did we find any bad intervals? Let's start with the longs...
if ~isempty(somethingwentwrongFEMALE) % Then there is a bad long interval in female data
    fprintf('found a wronglongfem. \n'); 
    for jj = 1:length(wronglongfem) 
        if ~isempty(wronglongfem{jj}) %identifies problem syllable type-
            fprintf('found a problem in syllable type %i. \n', jj);
            for kk = 1:length(wronglongfem{jj}) 
               if femints{jj}(femints{jj} ~= modint) == modint * 2 %finds double interval, the one we need to fix-
                fprintf('We have a winner - double interval. \n');
                transfersyl = (wronglongfem{jj}(kk) + modint); 
                
                % Eliminate the syllable from incorrect syllable type
                
                    for qq = 1:length(femtype)
                        femtype{qq} = femtype{qq}(femtype{qq} ~= transfersyl);
                    end
                    
                % Insert the syllables into the correct syllable type and make ascending order
                    femtype{jj} = sort([femtype{jj} transfersyl]); 
               end    
            end
        end
    end
end

fprintf('Testing male syllable intervals. \n');
    somethingwentwrongMALE = [];
    wrongshortmal = {}; wronglongmal = {}; % These will have the indices of the starting syllable for bad intervals

for pp = 1:length(malints) % Cycle through each syllable type for male
    wronglongmal{pp} = maltype{pp}([malints{pp}] > modint); % These have only long bad intervals

    if ~isempty(wronglongmal{pp}) 
        fprintf('Found %i long intervals in male syllable %i . \n', length(wronglongmal{pp}), pp);
        fprintf(' The syllables are: '); maltype{pp} 
        fprintf('\n The intervals are: '); malints{pp}  
        somethingwentwrongMALE = 1;
    end 
end

if ~isempty(somethingwentwrongMALE) % Then there is a bad long interval in male data
    fprintf('found a wronglongmal. \n'); 
    for jj = 1:length(wronglongmal) 
        if ~isempty(wronglongmal{jj}) %identifies problem syllable type-
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

figure(5);

    subplot(223); hold on; ylim([500 5000]);
        for i = 1:length(femtype)
            plot([foo.f(femtype{i}).sylen], [foo.f(femtype{i}).trace_peakf], clrs(i,:)); ylim([500 5000]);
        end;

    subplot(224); hold on; ylim([500 5000]);
        for i = 1:length(maltype)
            plot([foo.m(maltype{i}).sylen], [foo.m(maltype{i}).trace_peakf], clrs(i,:)); ylim([500 5000]);
        end;
        

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

%% Final review of syllable types.

    trF = foo.f(1).syltim(1);
    trM = foo.m(1).syltim(1);

    % A figure with all of the syllable traces
    figure(2); clf; xxa(1) = subplot(211); hold on; xxa(2) = subplot(212); hold on;

    for jj=1:length(foo.m)
                %plot([foo.f(jj).syltim(1) foo.f(jj).syltim(1)], [1000 4000], 'g'); 
                %plot([foo.m(jj).syltim(1) foo.m(jj).syltim(1)], [1000 4000], 'g'); 
                
            if find(femsylidx == jj)
                subplot(211); plot(foo.f(jj).trace_tim + foo.f(jj).syltim(1) -trF, foo.f(jj).trace_freq, 'm', 'LineWidth', 5);
                subplot(212); plot(foo.m(jj).trace_tim + foo.m(jj).syltim(1) -trM, foo.m(jj).trace_freq, 'm', 'LineWidth', 3);
                subplot(211); text(foo.f(jj).syltim(1) + 0.050 -trF, 750, num2str(jj), 'Color', 'm');
                subplot(212); text(foo.m(jj).syltim(1) + 0.050 -trM, 750, num2str(jj), 'Color', 'm');
            end
            if find(malsylidx == jj)
                subplot(211); plot(foo.f(jj).trace_tim + foo.f(jj).syltim(1) -trF, foo.f(jj).trace_freq, 'b', 'LineWidth', 3);
                subplot(212); plot(foo.m(jj).trace_tim + foo.m(jj).syltim(1) -trM, foo.m(jj).trace_freq, 'b', 'LineWidth', 5);
                subplot(211); text(foo.f(jj).syltim(1) + 0.050 -trF, 750, num2str(jj), 'Color', 'b');
                subplot(212); text(foo.m(jj).syltim(1) + 0.050 -trM, 750, num2str(jj), 'Color', 'b');
            end
            subplot(211); 
                plot(foo.f(jj).trace_tim + foo.f(jj).syltim(1) -trF, foo.f(jj).trace_freq, 'k', 'LineWidth', 1);
            subplot(212); 
                plot(foo.m(jj).trace_tim + foo.m(jj).syltim(1) -trM, foo.m(jj).trace_freq, 'k', 'LineWidth', 1);
                        
    end
        linkaxes(xxa, 'x');

        for yy = 1:length(femtype)
            for zz = 1:length(femtype{yy})
               figure(2); subplot(211); text(foo.f(femtype{yy}(zz)).syltim(1) + 0.050 -trF, 3750, num2str(yy), 'Color', 'k'); 
            end
        end
        for yy = 1:length(maltype)
            for zz = 1:length(maltype{yy})
               figure(2); subplot(212); text(foo.m(maltype{yy}(zz)).syltim(1) + 0.050 -trM, 3750, num2str(yy), 'Color', 'k'); 
            end
        end
        drawnow;
        close(4); close(5);
        
%% And we put the data into our output structure... probably should have done that much earlier.

anal.maltype = maltype;
anal.femtype = femtype;

fprintf('Be sure to save!  e.g. save filename.mat foo out anal \n');





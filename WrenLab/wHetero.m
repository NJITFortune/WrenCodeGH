function out = wHetero(in)
% Usage: Calculates response strength to solo and duet syllables.
% Relies on rs, a nested function below, to calculate Response Strength.
% Load the Chronic data structure first:
% load ChronicCompleat2017p.mat (OLD)
% load ChronicCompleat2018a.mat (Current as of 4-Aug-2018)

%% Setup

% How many bins do we want for our cycle? numsteps is into 360 degrees

    numsteps = 20;
    extrasteps = 10;

% Set up the degrees
    degreestep = 360 / numsteps;
    degreebase = -extrasteps*degreestep:degreestep:numsteps*degreestep+degreestep*(extrasteps-1);

% For time analysis    
    windowdur = 0.020; % millisecond window (meaning +/- msec from boundary)
    prepostwindows = 21; % half before and half after
    windowtims = -0.200:0.020:0.200;
        
% These are our descriptions of the syllables

[msolosyls, mduetsyls, fsolosyls, fduetsyls, spon] = wData;

% Initialize the bins for each segment of the cycle (and beyond!)
    femheterodegbins = zeros(1, numsteps+(2*extrasteps));
    malheterodegbins = zeros(1, numsteps+(2*extrasteps));

    malautodegbin = malheterodegbins; 
    femautodegbin = malheterodegbins;
    %femsolodegbin = femdegreebin; 
    %malsolodegbin = maldegreebin;
    
    fheterodeg(1).bins = femheterodegbins; 
    mheterodeg(1).bins = malheterodegbins;
    
    fheterotim(1).bins = zeros(1,prepostwindows); 
    mheterotim(1).bins = zeros(1,prepostwindows);
    malautotimbin(1).bins = zeros(1,prepostwindows); 
    femautotimbin(1).bins = zeros(1,prepostwindows); 

    fsolo(1).bins = femheterodegbins; 
    msolo(1).bins = malheterodegbins;
    mautodeg(1).bins = malheterodegbins; 
    fautodeg(1).bins = femheterodegbins;

    mspondeg = []; 
    fspondeg = []; 
    fspontim = [];
    mspontim = [];
    
    msolospon = []; 
    fsolospon = [];
    mautospondeg = []; 
    fautospondeg = [];    
    mautospontim = []; 
    fautospontim = [];    
    femtimbin = zeros(1,prepostwindows);
    maltimbin = zeros(1,prepostwindows);
    


    
%% Cycle through each pair    
for curpair = 1:length(spon) % Cycle for each pair
    
% DUET MALE Syllables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:length(mduetsyls{curpair}) % For each male duet syllable
        
        idx = length(fheterodeg)+1; % Add into the next entry in the structure
        cursylstart = in(curpair*2).syl(mduetsyls{curpair}(j)).tim(1); % Start time of current syllable
        cursylend = in(curpair*2).syl(mduetsyls{curpair}(j)).tim(2); % End time of current syllable
        curdur = cursylend - cursylstart; % Duration of the entire syllable
        curstepdur = curdur / numsteps; % Duration of our segment (normalizes to 360 degrees)

        % Cycle through each degree segment  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for k = -extrasteps:numsteps+extrasteps-1
            
            tmp = 0; spontmp = 0; autotmp = 0; sponautotmp = 0; 
            % This picks a random window within the spontaneous window with
            % the same duration as the syllable
                sponstart = spon(1,curpair) + ((abs(spon(1,curpair) - spon(2,curpair)) - curstepdur) * rand);
                sponend = sponstart + curstepdur;
            
            % HETEROGENOUS (male syllables, female spikes)    
            for i=1:4 % 4 electrodes in a tetrode always (CHRONIC DATA ONLY)
                % Simply sum up the number of spikes in the window.
                tmp = tmp + length(find(in(curpair*2).Cspikes{i} > cursylstart + curstepdur*k ...
                    & in(curpair*2).Cspikes{i} < cursylstart + curstepdur*(k+1)));
                
                femheterodegbins(k+extrasteps+1) = femheterodegbins(k+extrasteps+1) + tmp;
                
                spontmp = spontmp + length(find(in(curpair*2).Cspikes{i} > sponstart ...
                    & in(curpair*2).Cspikes{i} < sponend));
            end
            
                fheterodeg(idx).bins(k+extrasteps+1) = tmp;
                fspondeg(end+1) = spontmp;
            
            % AUTOGENOUS (male syllables, male spikes)
            for i=1:4 % 4 electrodes in a tetrode always
                autotmp = autotmp + length(find(in((curpair*2)-1).Cspikes{i} > cursylstart + curstepdur*k ...
                    & in((curpair*2)-1).Cspikes{i} < cursylstart + curstepdur*(k+1)));
                
                malautodegbin(k+extrasteps+1) = malautodegbin(k+extrasteps+1) + autotmp;
                
                sponautotmp = sponautotmp + length(find(in((curpair*2)-1).Cspikes{i} > sponstart ...
                    & in((curpair*2)-1).Cspikes{i} < sponend));
            end
                    
                mautodeg(idx).bins(k+extrasteps+1) = autotmp;
                mautospondeg(end+1) = sponautotmp;

        end  % END OF DEGREE ANALYSIS
        
        % Cycle through time-based analysis  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for k = 1:prepostwindows % 20 segments (10 before and 10 after) around start and end syllables
            tmp = 0; spontmp = 0; autotmp = 0; sponautotmp = 0; 
            % This picks a random window within the spontaneous window with
            % the same duration as the syllable
                sponstart = spon(1,curpair) + ((abs(spon(1,curpair) - spon(2,curpair)) - windowdur) * rand);
                sponend = sponstart + windowdur;
                
            % HETEROGENOUS (male syllables, female spikes)    
            for i=1:4 % 4 electrodes in a tetrode always (CHRONIC DATA ONLY)
                % Simply sum up the number of spikes in the window.
                % fembin is the sum of all (across duets) 
                tmp = tmp + length(find(in(curpair*2).Cspikes{i} > cursylstart-(prepostwindows*(windowdur/2))+windowdur*(k-1) ...
                    & in(curpair*2).Cspikes{i} < cursylstart-(prepostwindows*(windowdur/2))+windowdur*k));
                
                femtimbin(k) = femtimbin(k) + tmp;
               
                spontmp = spontmp + length(find(in(curpair*2).Cspikes{i} > sponstart ...
                    & in(curpair*2).Cspikes{i} < sponend));
            end
            
                fheterotim(idx).bins(k) = tmp;
                fspontim(end+1) = spontmp;
            
            % AUTOGENOUS (male syllables, male spikes)
            for i=1:4 % 4 electrodes in a tetrode always
                autotmp = autotmp + length(find(in((curpair*2)-1).Cspikes{i} > (cursylstart-(prepostwindows*windowdur/2)) + curstepdur*(k-1) ...
                    & in((curpair*2)-1).Cspikes{i} < (cursylstart-(prepostwindows*windowdur/2)) + curstepdur*k));

                malautotimbin(k) = malautotimbin(k) + autotmp;
                
                sponautotmp = sponautotmp + length(find(in((curpair*2)-1).Cspikes{i} > sponstart ...
                    & in((curpair*2)-1).Cspikes{i} < sponend));
            end
                    
                mautotim(idx).bins(k) = autotmp;
                mautospontim(end+1) = sponautotmp;
        
        end
        
        
    end % End of male duet syllables
        
% DUET FEMALE Syllables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:length(fduetsyls{curpair}) % Female duet syllables
        
        idx = length(mheterodeg)+1;
        cursylstart = in((curpair*2)-1).syl(fduetsyls{curpair}(j)).tim(1);
        cursylend = in((curpair*2)-1).syl(fduetsyls{curpair}(j)).tim(2);
        curdur = cursylend - cursylstart;
        curstepdur = curdur / numsteps;
        
        for k = -extrasteps:numsteps+extrasteps-1    
            
            tmp = 0; spontmp = 0;  autotmp = 0; sponautotmp = 0;
            sponstart = spon(1,curpair) + ((abs(spon(1,curpair) - spon(2,curpair)) - curstepdur) * rand);
            sponend = sponstart + curstepdur;
            
            % HETEROGENOUS (female syllables, male spikes)    
            for i=1:4 % 4 electrodes in a tetrode always
                malheterodegbins(k+extrasteps+1) = malheterodegbins(k+extrasteps+1) + length(find(in((curpair*2)-1).Cspikes{i} > cursylstart + curstepdur*k ...
                    & in((curpair*2)-1).Cspikes{i} < cursylstart + curstepdur*(k+1)));
                tmp = tmp + length(find(in((curpair*2)-1).Cspikes{i} > cursylstart + curstepdur*k ...
                    & in((curpair*2)-1).Cspikes{i} < cursylstart + curstepdur*(k+1)));
                spontmp = spontmp + length(find(in(curpair*2).Cspikes{i} > sponstart ...
                    & in(curpair*2).Cspikes{i} < sponend));
            end
            
                mheterodeg(idx).bins(k+extrasteps+1) = tmp;
                mspondeg(end+1) = spontmp;
                
            % AUTOGENOUS (female syllables, female spikes)
            for i=1:4 % 4 electrodes in a tetrode always
                femautodegbin(k+extrasteps+1) = femautodegbin(k+extrasteps+1) + length(find(in(curpair*2).Cspikes{i} > cursylstart + curstepdur*k ...
                    & in(curpair*2).Cspikes{i} < cursylstart + curstepdur*(k+1)));
                autotmp = autotmp + length(find(in(curpair*2).Cspikes{i} > cursylstart + curstepdur*k ...
                    & in(curpair*2).Cspikes{i} < cursylstart + curstepdur*(k+1)));
                sponautotmp = sponautotmp + length(find(in(curpair*2).Cspikes{i} > sponstart ...
                    & in(curpair*2).Cspikes{i} < sponend));
            end
                    
                fautodeg(idx).bins(k+extrasteps+1) = autotmp;
                fautospondeg(end+1) = sponautotmp;                
                
        end        
        
        % Cycle through time-based analysis  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for k = 1:prepostwindows % 20 segments (10 before and 10 after) around start and end syllables
            tmp = 0; spontmp = 0; autotmp = 0; sponautotmp = 0; 
            % This picks a random window within the spontaneous window with
            % the same duration as the syllable
                sponstart = spon(1,curpair) + ((abs(spon(1,curpair) - spon(2,curpair)) - windowdur) * rand);
                sponend = sponstart + windowdur;
                
            % HETEROGENOUS (female syllables, male spikes)    
            for i=1:4 % 4 electrodes in a tetrode always (CHRONIC DATA ONLY)
                % Simply sum up the number of spikes in the window.
                % fembin is the sum of all (across duets) 
                tmp = tmp + length(find(in((curpair*2)-1).Cspikes{i} > cursylstart-(prepostwindows*(windowdur/2))+windowdur*(k-1) ...
                    & in((curpair*2)-1).Cspikes{i} < cursylstart-(prepostwindows*(windowdur/2))+windowdur*k));
                
                maltimbin(k) = maltimbin(k) + tmp;
               
                spontmp = spontmp + length(find(in((curpair*2)-1).Cspikes{i} > sponstart ...
                    & in((curpair*2)-1).Cspikes{i} < sponend));
            end
            
                mheterotim(idx).bins(k) = tmp;
                mspontim(end+1) = spontmp;
            
            % AUTOGENOUS (female syllables, female spikes)
            for i=1:4 % 4 electrodes in a tetrode always
                autotmp = autotmp + length(find(in(curpair*2).Cspikes{i} > (cursylstart-(prepostwindows*windowdur/2)) + curstepdur*(k-1) ...
                    & in(curpair*2).Cspikes{i} < (cursylstart-(prepostwindows*windowdur/2)) + curstepdur*k));

                femautotimbin(k) = femautotimbin(k) + autotmp;
                
                sponautotmp = sponautotmp + length(find(in(curpair*2).Cspikes{i} > sponstart ...
                    & in(curpair*2).Cspikes{i} < sponend));
            end
                    
                fautotim(idx).bins(k) = autotmp;
                fautospontim(end+1) = sponautotmp;
        
        end
        
    end % End of female duet syllables
    
    
    
%     % SOLO Syllables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for j = 1:length(msolosyls{curpair}) % Male solo syllables
%         
%         idx = length(f)+1;
%         cursylstart = in(curpair*2).syl(msolosyls{curpair}(j)).tim(1);
%         cursylend = in(curpair*2).syl(msolosyls{curpair}(j)).tim(2);
%         curdur = cursylend - cursylstart;
%         curstepdur = curdur / numsteps;
%         
%         for k = -extrasteps:numsteps+extrasteps-1
%             tmp = 0; spontmp = 0;
%             sponstart = spon(1,curpair) + ((abs(spon(1,curpair) - spon(2,curpair)) - curstepdur) * rand);
%             sponend = sponstart + curstepdur;
%             for i=1:4 % 4 electrodes in a tetrode always
%                 femsolobin(k+extrasteps+1) = femsolobin(k+extrasteps+1) + length(find(in(curpair*2).Cspikes{i} > cursylstart + curstepdur*k ...
%                     & in(curpair*2).Cspikes{i} < cursylstart + curstepdur*(k+1)));
%                 tmp = tmp + length(find(in(curpair*2).Cspikes{i} > cursylstart + curstepdur*k ...
%                     & in(curpair*2).Cspikes{i} < cursylstart + curstepdur*(k+1)));
%                 spontmp = spontmp + length(find(in(curpair*2).Cspikes{i} > sponstart ...
%                     & in(curpair*2).Cspikes{i} < sponend));
%             end
%             fsolo(idx).bins(k+extrasteps+1) = tmp;
%             fsolospon(end+1) = spontmp;
%         end        
%         
%     end % End of male duet syllables    
%     
%     for j = 1:length(fsolosyls{curpair}) % Female solo syllables
%         
%         idx = length(m)+1;
%         cursylstart = in((curpair*2)-1).syl(fsolosyls{curpair}(j)).tim(1);
%         cursylend = in((curpair*2)-1).syl(fsolosyls{curpair}(j)).tim(2);
%         curdur = cursylend - cursylstart;
%         curstepdur = curdur / numsteps;
%         
%         for k = -extrasteps:numsteps+extrasteps-1         
%             tmp = 0; spontmp = 0;
%             sponstart = spon(1,curpair) + ((abs(spon(1,curpair) - spon(2,curpair)) - curstepdur) * rand);
%             sponend = sponstart + curstepdur;
%             for i=1:4 % 4 electrodes in a tetrode always
%                 malsolobin(k+extrasteps+1) = malsolobin(k+extrasteps+1) + length(find(in((curpair*2)-1).Cspikes{i} > cursylstart + curstepdur*k ...
%                     & in((curpair*2)-1).Cspikes{i} < cursylstart + curstepdur*(k+1)));
%                 tmp = tmp + length(find(in((curpair*2)-1).Cspikes{i} > cursylstart + curstepdur*k ...
%                     & in((curpair*2)-1).Cspikes{i} < cursylstart + curstepdur*(k+1)));
%                 spontmp = spontmp + length(find(in(curpair*2).Cspikes{i} > sponstart ...
%                     & in(curpair*2).Cspikes{i} < sponend));
%             end
%             msolo(idx).bins(k+extrasteps+1) = tmp;
%             msolospon(end+1) = spontmp;
%         end        
%         
%     end % End of female solo syllables    
    
end % curpair (cycle through spons)


    out.mspondeg = mspondeg;
    out.fspondeg = fspondeg;
    out.malbindeg = malheterodegbins;
    out.fembindeg = femheterodegbins;

    out.mspontim = mspontim;
    out.fspondeg = fspontim;
    
    out.malbintim = maltimbin;
    out.fembintim = femtimbin;
            figure(5); clf; plot(windowtims, out.fembintim, '*-m');
            hold on; ot(windowtims, out.malbintim, '*-b');
    
    out.mautospondeg = mautospondeg;
    out.malautobindeg = malautodegbin;
    out.fautospondeg = fautospondeg;
    out.fautobindeg = femautodegbin;

%     out.msolospon = msolospon;
%     out.fsolospon = fsolospon;
%     out.malsolobin = malsolobin;
%     out.femsolobin = femsolobin;
        
    
    guessfspon = sum(fspondeg) / (numsteps+extrasteps);
    guessmspon = sum(mspondeg) / (numsteps+extrasteps);
    
    figure(1); clf; % Separate plots for HETEROGENOUS
    
        subplot(121); plot(degreebase, out.fembindeg, '*-m'); hold on;
            plot([0, 0], [1, max(out.fembindeg)], 'k-'); plot([360, 360], [1, max(out.fembindeg)], 'k-');
            plot([0, 360], [guessfspon, guessfspon], 'r-');
            
        subplot(122); plot(degreebase, out.malbindeg, '*-b'); hold on;   
            plot([0, 0], [1, max(out.malbindeg)], 'k-'); plot([360, 360], [1, max(out.malbindeg)], 'k-'); 
            plot([0, 360], [guessmspon, guessmspon], 'c-');
            
    figure(2); clf; % Single plot for HETEROGENOUS
            plot(degreebase, out.fembindeg/max(out.fembindeg), '*-m'); 
            hold on;
            plot(degreebase, out.malbindeg/max(out.malbindeg), '*-b');

            plot([0, 0], [0, 1], 'k-'); plot([360 360], [0 1], 'k-');
            
            plot([0 360], [guessfspon/max(out.fembindeg), guessfspon/max(out.fembindeg)], 'r-');
            plot([0 360], [guessmspon/max(out.malbindeg), guessmspon/max(out.malbindeg)], 'c-');
            
    figure(3); clf; % Separate plots for AUTOGENOUS
        subplot(121); plot(degreebase, out.fautobindeg, '*-m'); hold on;
            plot([0, 0], [1, max(out.fautobindeg)], 'k-'); plot([360, 360], [1, max(out.fautobindeg)], 'k-');
            plot([0 360], [guessfspon, guessfspon], 'r-');
            
        subplot(122); plot(degreebase, out.malautobindeg, '*-b'); hold on;   
            plot([0 0], [1 max(out.malautobindeg)], 'k-'); plot([360 360], [1 max(out.malautobindeg)], 'k-'); 
            plot([0 360], [guessmspon, guessmspon], 'c-');
            
    figure(4); clf; % Single plot for AUTOGENOUS
            plot(degreebase, out.fautobindeg/max(out.fautobindeg), '*-m'); 
            hold on;
            plot(degreebase, out.malautobindeg/max(out.malautobindeg), '*-b');

            plot([0, 0], [0, 1], 'k-'); plot([360 360], [0 1], 'k-');
            
            plot([0 360], [guessfspon/max(out.fautobindeg), guessfspon/max(out.fautobindeg)], 'r-');
            plot([0 360], [guessmspon/max(out.malautobindeg), guessmspon/max(out.malautobindeg)], 'c-');

            
            
allsteps = length(fheterodeg(1).bins);
% for jj = 1:allsteps     
%     for kk = 2:length(f)
%         z(kk-1) = f(kk).bins(jj);
%     end        
%     out.medianf(jj) = median(z); 
%     out.stdf(jj) = std(z);
%     out.modef(jj) = mode(z);
% end
% allsteps = length(m(1).bins);
% for jj = 1:allsteps     
%     for kk = 2:length(m)
%         z(kk-1) = m(kk).bins(jj);
%     end        
%     out.medianm(jj) = median(z); 
%     out.stdm(jj) = std(z);
%     out.modem(jj) = mode(z);
% end


% figure(2); clf;
%         subplot(121); plot(degreebase, out.medianf, '*-m'); hold on;
%             %plot([0 0], [0 2], 'k-'); plot([360 360], [0 2], 'k-'); 
%         subplot(122); plot(degreebase, out.medianm, '*-b'); hold on;   
%             %plot([0 0], [0 3], 'k-'); plot([360 360], [0 3], 'k-'); 
% 
% figure(3); clf;
%         subplot(121); plot(degreebase, out.modef, '*-m'); hold on;
%             %plot([0 0], [0 2], 'k-'); plot([360 360], [0 2], 'k-'); 
%         subplot(122); plot(degreebase, out.modem, '*-b'); hold on;   
%             %plot([0 0], [0 3], 'k-'); plot([360 360], [0 3], 'k-'); 
            
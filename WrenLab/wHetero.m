function out = wHetero(in, padding)
% Usage: Calculates response strength to solo and duet syllables.
% Relies on rs, a nested function below, to calculate Response Strength.
% Load the Chronic data structure first:
% load ChronicCompleat2017p.mat (OLD)
% load ChronicCompleat2018a.mat (Current as of 4-Aug-2018)

% 'pad' is a critical variable - it is the shift in the window around the
% clicked boundaries of the syllables for the calculation of RS. 
% The boundaries of each syllable are used for this analysis. This is
% inherently problematic as we expect pre-motor activity in awake animals to occur PRIOR to
% the sound and auditory activity in urethane-anesthetized animals to occur AFTER the sound.
% We are comfortable with using a value of '0' it is a rather unbiased. Change
% the value of padding in seconds (e.g. 0.020 or -0.030) to look at the effects on the results.
pad = 0.000; 

% analpad = 0.050;
numsteps = 20;
extrasteps = 10;

degreestep = 360 / numsteps;
degreebase = -extrasteps*degreestep:degreestep:numsteps*degreestep+degreestep*(extrasteps-1);

% The user can specify the padding via an argin for convenience.
if nargin == 2; pad = padding; end

[msolosyls, mduetsyls, fsolosyls, fduetsyls, spon] = wData;

    fembin = zeros(1, numsteps+(2*extrasteps));
    malbin = zeros(1, numsteps+(2*extrasteps));
    femsolobin = fembin; malsolobin = malbin;
    f(1).bins = fembin;
    m(1).bins = malbin;
    mspon = []; fspon = []; msolospon = []; fsolospon = [];
    fsolo(1).bins = fembin;
    msolo(1).bins = malbin;

for curpair = 1:length(spon) % Cycle for each pair
    
    % DUET Syllables
    for j = 1:length(mduetsyls{curpair}) % Male duet syllables
        
        idx = length(f)+1;
        cursylstart = in(curpair*2).syl(mduetsyls{curpair}(j)).tim(1);
        cursylend = in(curpair*2).syl(mduetsyls{curpair}(j)).tim(2);
        curdur = cursylend - cursylstart;
        curstepdur = curdur / numsteps;
        
        for k = -extrasteps:numsteps+extrasteps-1
            tmp = 0; spontmp = 0;
            sponstart = spon(1,curpair) + ((abs(spon(1,curpair) - spon(2,curpair)) - curstepdur) * rand);
            sponend = sponstart + curstepdur;
            for i=1:4 % 4 electrodes in a tetrode always
                fembin(k+extrasteps+1) = fembin(k+extrasteps+1) + length(find(in(curpair*2).Cspikes{i} > cursylstart + curstepdur*k ...
                    & in(curpair*2).Cspikes{i} < cursylstart + curstepdur*(k+1)));
                tmp = tmp + length(find(in(curpair*2).Cspikes{i} > cursylstart + curstepdur*k ...
                    & in(curpair*2).Cspikes{i} < cursylstart + curstepdur*(k+1)));
                spontmp = spontmp + length(find(in(curpair*2).Cspikes{i} > sponstart ...
                    & in(curpair*2).Cspikes{i} < sponend));
            end
            f(idx).bins(k+extrasteps+1) = tmp;
            fspon(end+1) = spontmp;
        end        
        
    end % End of male duet syllables
        
    for j = 1:length(fduetsyls{curpair}) % Female duet syllables
        
        idx = length(m)+1;
        cursylstart = in((curpair*2)-1).syl(fduetsyls{curpair}(j)).tim(1);
        cursylend = in((curpair*2)-1).syl(fduetsyls{curpair}(j)).tim(2);
        curdur = cursylend - cursylstart;
        curstepdur = curdur / numsteps;
        
        for k = -extrasteps:numsteps+extrasteps-1         
            tmp = 0; spontmp = 0;
            sponstart = spon(1,curpair) + ((abs(spon(1,curpair) - spon(2,curpair)) - curstepdur) * rand);
            sponend = sponstart + curstepdur;
            for i=1:4 % 4 electrodes in a tetrode always
                malbin(k+extrasteps+1) = malbin(k+extrasteps+1) + length(find(in((curpair*2)-1).Cspikes{i} > cursylstart + curstepdur*k ...
                    & in((curpair*2)-1).Cspikes{i} < cursylstart + curstepdur*(k+1)));
                tmp = tmp + length(find(in((curpair*2)-1).Cspikes{i} > cursylstart + curstepdur*k ...
                    & in((curpair*2)-1).Cspikes{i} < cursylstart + curstepdur*(k+1)));
                spontmp = spontmp + length(find(in(curpair*2).Cspikes{i} > sponstart ...
                    & in(curpair*2).Cspikes{i} < sponend));
            end
            m(idx).bins(k+extrasteps+1) = tmp;
            mspon(end+1) = spontmp;
        end        
        
    end % End of female duet syllables
    
    % SOLO Syllables
    for j = 1:length(msolosyls{curpair}) % Male solo syllables
        
        idx = length(f)+1;
        cursylstart = in(curpair*2).syl(msolosyls{curpair}(j)).tim(1);
        cursylend = in(curpair*2).syl(msolosyls{curpair}(j)).tim(2);
        curdur = cursylend - cursylstart;
        curstepdur = curdur / numsteps;
        
        for k = -extrasteps:numsteps+extrasteps-1
            tmp = 0; spontmp = 0;
            sponstart = spon(1,curpair) + ((abs(spon(1,curpair) - spon(2,curpair)) - curstepdur) * rand);
            sponend = sponstart + curstepdur;
            for i=1:4 % 4 electrodes in a tetrode always
                femsolobin(k+extrasteps+1) = femsolobin(k+extrasteps+1) + length(find(in(curpair*2).Cspikes{i} > cursylstart + curstepdur*k ...
                    & in(curpair*2).Cspikes{i} < cursylstart + curstepdur*(k+1)));
                tmp = tmp + length(find(in(curpair*2).Cspikes{i} > cursylstart + curstepdur*k ...
                    & in(curpair*2).Cspikes{i} < cursylstart + curstepdur*(k+1)));
                spontmp = spontmp + length(find(in(curpair*2).Cspikes{i} > sponstart ...
                    & in(curpair*2).Cspikes{i} < sponend));
            end
            fsolo(idx).bins(k+extrasteps+1) = tmp;
            fsolospon(end+1) = spontmp;
        end        
        
    end % End of male duet syllables    
    
    for j = 1:length(fsolosyls{curpair}) % Female solo syllables
        
        idx = length(m)+1;
        cursylstart = in((curpair*2)-1).syl(fsolosyls{curpair}(j)).tim(1);
        cursylend = in((curpair*2)-1).syl(fsolosyls{curpair}(j)).tim(2);
        curdur = cursylend - cursylstart;
        curstepdur = curdur / numsteps;
        
        for k = -extrasteps:numsteps+extrasteps-1         
            tmp = 0; spontmp = 0;
            sponstart = spon(1,curpair) + ((abs(spon(1,curpair) - spon(2,curpair)) - curstepdur) * rand);
            sponend = sponstart + curstepdur;
            for i=1:4 % 4 electrodes in a tetrode always
                malsolobin(k+extrasteps+1) = malsolobin(k+extrasteps+1) + length(find(in((curpair*2)-1).Cspikes{i} > cursylstart + curstepdur*k ...
                    & in((curpair*2)-1).Cspikes{i} < cursylstart + curstepdur*(k+1)));
                tmp = tmp + length(find(in((curpair*2)-1).Cspikes{i} > cursylstart + curstepdur*k ...
                    & in((curpair*2)-1).Cspikes{i} < cursylstart + curstepdur*(k+1)));
                spontmp = spontmp + length(find(in(curpair*2).Cspikes{i} > sponstart ...
                    & in(curpair*2).Cspikes{i} < sponend));
            end
            msolo(idx).bins(k+extrasteps+1) = tmp;
            msolospon(end+1) = spontmp;
        end        
        
    end % End of female solo syllables    
    
end % curpair (cycle through spons)


    out.mspon = mspon;
    out.fspon = fspon;
    out.malbin = malbin;
    out.fembin = fembin;

    guessfspon = sum(fspon) / (numsteps+extrasteps);
    guessmspon = sum(mspon) / (numsteps+extrasteps);
    
    figure(1); clf; 
        subplot(121); plot(degreebase, out.fembin, '*-m'); hold on;
            plot([0 0], [1 140], 'k-'); plot([360 360], [1 140], 'k-');
            plot([0 360], [guessfspon, guessfspon], 'r-');
        subplot(122); plot(degreebase, out.malbin, '*-b'); hold on;   
            plot([0 0], [1 180], 'k-'); plot([360 360], [1 180], 'k-'); 
            plot([0 360], [guessmspon, guessmspon], 'c-');
            
    figure(2); clf; 
            plot(degreebase, out.fembin/max(out.fembin), '*-m'); 
            hold on;
            plot(degreebase, out.malbin/max(out.malbin), '*-b');
            plot([0 0], [0 1], 'k-'); plot([360 360], [0 1], 'k-');
            plot([0 360], [guessfspon/max(out.fembin), guessfspon/max(out.fembin)], 'r-');
            plot([0 360], [guessmspon/max(out.malbin), guessmspon/max(out.malbin)], 'c-');
            
            
            
allsteps = length(f(1).bins);
for jj = 1:allsteps     
    for kk = 2:length(f)
        z(kk-1) = f(kk).bins(jj);
    end        
    out.medianf(jj) = median(z); 
    out.stdf(jj) = std(z);
    out.modef(jj) = mode(z);
end
allsteps = length(m(1).bins);
for jj = 1:allsteps     
    for kk = 2:length(m)
        z(kk-1) = m(kk).bins(jj);
    end        
    out.medianm(jj) = median(z); 
    out.stdm(jj) = std(z);
    out.modem(jj) = mode(z);
end


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
            
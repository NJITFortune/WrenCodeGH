% load Wren2015acute_versionS.mat

% Road is 1-120
% Rock is 121-656
% Lucky is 657-1125

% Examples:
% Get the indices for Male Road to stim 1
% intersect(intersect([idx.Road], [idx.M]), find([idx.stimidx] == 1))

% Get the indices for Female Road to stim 1
% intersect(intersect([idx.Road], [idx.F]), find([idx.stimidx] == 1))

%% Get max Z for each neuron
    fstd(length(idx.Funit)).Zs = []; % Preallocation for speed
    mstd(length(idx.Munit)).Zs = []; % Preallocation for speed

% Females
fneuronmaxZ = zeros(1,length(idx.Funit));
fneuronminZ = ones(1,length(idx.Funit))*10;

    for j=1:length(idx.Funit) 
        for k=1:length(idx.Funit{j})
           fneuronmaxZ(j) = max([fneuronmaxZ(j), w(idx.Funit{j}(k)).Z]);
           fneuronminZ(j) = min([fneuronminZ(j), w(idx.Funit{j}(k)).Z]);
           fstd(j).Zs(end+1) = w(idx.Funit{j}(k)).Z;
        end
    end
    
    % Males
mneuronmaxZ = zeros(1,length(idx.Munit));
mneuronminZ = ones(1,length(idx.Munit))*10;

    for j=1:length(idx.Munit)
        for k=1:length(idx.Munit{j})
           mneuronmaxZ(j) = max([mneuronmaxZ(j), w(idx.Munit{j}(k)).Z]);
           mneuronminZ(j) = min([mneuronminZ(j), w(idx.Munit{j}(k)).Z]);
           mstd(j).Zs(end+1) = w(idx.Munit{j}(k)).Z;
        end
    end


% What I found - there were some neurons that showed no variation in Z
% value across stimuli.  Let's eliminate those as non-responders. What is
% cool is that the stronger the max Z, the wider the range of Zs to
% different stimuli, suggesting selectivity!

goodMaleneurons = find(mneuronmaxZ-mneuronminZ > 0.1);
goodFemaleneurons = find(fneuronmaxZ-fneuronminZ > 0.1);

figure(1); clf; title('Max versus STD Zs'); hold on;
for j=1:length(goodFemaleneurons)
    dotsize = 1*length(idx.Funit{goodFemaleneurons(j)});
%    plot(fneuronmaxZ(goodFemaleneurons(j)), fneuronmaxZ(goodFemaleneurons(j))-fneuronminZ(goodFemaleneurons(j)), 'mo', 'MarkerSize', dotsize); 
    plot(fneuronmaxZ(goodFemaleneurons(j)), std(fstd(goodFemaleneurons(j)).Zs), 'mo', 'MarkerSize', dotsize); 
end
for j=1:length(goodMaleneurons)
    dotsize = 1*length(idx.Funit{goodMaleneurons(j)});
%    plot(mneuronmaxZ(goodMaleneurons(j)), mneuronmaxZ(goodMaleneurons(j))-mneuronminZ(goodMaleneurons(j)), 'bd', 'MarkerSize', dotsize);
    plot(mneuronmaxZ(goodMaleneurons(j)), std(mstd(goodMaleneurons(j)).Zs), 'bd', 'MarkerSize', dotsize); 
end
%plot([0,6], [0,6], 'k-')
plot([0,6], [0,2], 'k-')
xlabel('Max Z score');
% ylabel('Z score range');
ylabel('Standard deviation, Zs');

% SfN2019-Zranges.eps

%% Do the pairs

% Road is w 1-120
% Funits up to 10, Munits up to 10
% BV stim 7

RoadList = [1 2 3 4 5 6 8 9 10 11 12]; % BV is missing
stimRoad(2000).Zf = []; stimRoad(2000).Zm = [];
for j = idx.Road % For every Road entry    
    if w(j).sex == 'F'   
        stimRoad(w(j).stimidx).Zf(end+1) = w(j).Z;
    else
        stimRoad(w(j).stimidx).Zm(end+1) = w(j).Z;        
    end    
end

RockList = [13 14 15 16 17 18 19];
stimRock(2000).Zf = []; stimRock(2000).Zm = [];
for j = idx.Rock % For every Rock entry    
    if w(j).sex == 'F'        
        stimRock(w(j).stimidx).Zf(end+1) = w(j).Z;
    else
        stimRock(w(j).stimidx).Zm(end+1) = w(j).Z;        
    end    
end

LuckyList = [20 21 22 23 24 25 26 27 28 29 30 31 32];
stimLucky(2000).Zf = []; stimLucky(2000).Zm = [];
for j = idx.Lucky % For every Lucky entry    
    if w(j).sex == 'F'        
        stimLucky(w(j).stimidx(1)).Zf(end+1) = w(j).Z;
    else
        stimLucky(w(j).stimidx(1)).Zm(end+1) = w(j).Z;        
    end    
end

figure(2); clf; title('Correlation between Female and Male Zs');
xlabel('Male Z');
ylabel('Female Z')
hold on;
allZs = [];
for k=1:length(RoadList)
    plot(mean(stimRoad(RoadList(k)).Zm), mean(stimRoad(RoadList(k)).Zf), 'bo', 'MarkerSize', 8);
    allZs(end+1,:) = [mean(stimRoad(RoadList(k)).Zm), mean(stimRoad(RoadList(k)).Zf)];
end
for k=1:length(RockList)
    plot(mean(stimRock(RockList(k)).Zm), mean(stimRock(RockList(k)).Zf), 'rs', 'MarkerSize', 8);
    allZs(end+1,:) = [mean(stimRock(RockList(k)).Zm), mean(stimRock(RockList(k)).Zf)];
end
for k=1:length(LuckyList)
    plot(mean(stimLucky(LuckyList(k)).Zm), mean(stimLucky(LuckyList(k)).Zf), 'md', 'MarkerSize', 8);
    allZs(end+1,:) = [mean(stimLucky(LuckyList(k)).Zm), mean(stimLucky(LuckyList(k)).Zf)];
end
xlim([-0.5, 3]); ylim([-0.5 3]);
plot([-0.5, 3], [-0.5, 3], 'k-');
% SfN2019-Zcorrelation.eps

figure(3); clf; hold on; title('Correlation between Female and Male Zs');
for j=1:length(allZs); plot(allZs(j,1), allZs(j,2), 'ko'); end
xlim([-0.5, 3]); ylim([-0.5 3]);

mdl = fitlm(allZs(:,1), allZs(:,2))


%% Do we want to do something with BV?  That seems useful.



Zthresh = 0.5;

% Get list of valid stims (should be same for both males and females)
stims = unique(idx.stimidx(intersect([idx.Lucky], [idx.F])));

% Get unique pairings
c = combnk(1:length(stims), 2); % c(p,:)

FdPrime = []; MdPrime = [];

for zz = 1:length(c)
    
Ftmp1spikes = {}; Ftmp2spikes = {};
Mtmp1spikes = {}; Mtmp2spikes = {};

Fi = intersect(intersect([idx.Lucky], [idx.F]), find([idx.stimidx] == stims(c(zz,1)) ));
Fii = intersect(intersect([idx.Lucky], [idx.F]), find([idx.stimidx] == stims(c(zz,2)) ));
Mi = intersect(intersect([idx.Lucky], [idx.M]), find([idx.stimidx] == stims(c(zz,1)) ));
Mii = intersect(intersect([idx.Lucky], [idx.M]), find([idx.stimidx] == stims(c(zz,2)) ));

% Get the length of the stimuli from the first instance
Fistimend =  w(Fi(1)).stimend; Fiistimend =  w(Fii(1)).stimend;
Mistimend =  w(Mi(1)).stimend; Miistimend =  w(Mii(1)).stimend;

% Females first
for j=1:length(Fi)
    if w(Fi(j)).Z > Zthresh
    for k = 1:length(w(Fi(j)).spikes)
        Ftmp1spikes{end+1} = w(Fi(j)).spikes{k};
    end
    end
end
for j=1:length(Fii)
    if w(Fii(j)).Z > Zthresh
    for k = 1:length(w(Fii(j)).spikes)
        Ftmp2spikes{end+1} = w(Fii(j)).spikes{k};
    end
    end
end

% Now the males
for j=1:length(Mi)
    if w(Mi(j)).Z > Zthresh
    for k = 1:length(w(Mi(j)).spikes)
        Mtmp1spikes{end+1} = w(Mi(j)).spikes{k};
    end
    end
end
for j=1:length(Mii)
    if w(Mii(j)).Z > Zthresh
    for k = 1:length(w(Mii(j)).spikes)
        Mtmp2spikes{end+1} = w(Mii(j)).spikes{k};
    end
    end
end

% length(Ftmp1spikes)
% length(Ftmp2spikes)
% length(Mtmp1spikes)
% length(Mtmp2spikes)

FdPrime(end+1) = bs_DP( Ftmp1spikes, [0 Fistimend], Ftmp2spikes, [0 Fiistimend]);
MdPrime(end+1) = bs_DP( Mtmp1spikes, [0 Mistimend], Mtmp2spikes, [0 Miistimend]);

clear Ftmp1spikes Ftmp2spikes Mtmp1spikes Mtmp2spikes Fi Fii Mi Mii k j 

end

 
%% Let's just do the Z-score thing

zThresh = 1; Rockstims = [];

% For each Rock neuron
for j=1:length(idx.Rock)
    a = zeros(1,length(s));
    for k=1:length(s) 
        a(k) = ~isempty(strfind(s(k).stimname, w(idx.Rock(j)).stimname));
    end
        Rockstims(end+1) = find(a);
end


%% Set the stim IDXs for all
for j=1:length(w)
    a = zeros(1,length(s));
    for k=1:length(s) 
        a(k) = ~isempty(strfind(s(k).stimname, w(j).stimname));
    end
    w(j).stimidx = find(a);
end

%% Plot duet switch that occurs in the 'fieldduet' stimulus

figure(4); clf; subplot(211); % Saved as fielduet-switch-sfn2019poster.eps
    specgram(s(13).stim, 2048, s(13).Fs, [], 2024); ylim([500 4000]);
    caxis([18 45]); xlim([19 26]); colormap(flipud(gray));
    

%% XCORR Business

sig1 = w(10).hst.spers / max(w(10).hst.spers); sig1 = sig1-mean(sig1);
sig2 = w(16).hst.spers / max(w(16).hst.spers); sig2 = sig2-mean(sig2);
sig3 = w(64).hst.spers / max(w(64).hst.spers); sig3 = sig3-mean(sig3);

%pxy = cpsd(w(10).hst.spers, w(64).hst.spers); 
%[pxy,f] = cpsd(w(10).hst.spers, w(64).hst.spers,[],[],[],w(64).hst.Fs);
[pxy,f] = cpsd(sig1, sig3,[],[],[],w(64).hst.Fs);
figure(5); clf; plot(f, real(pxy), 'o-'); xlim([0 8]); %% Peak at 2 Hz for alternation of syllables

figure(6); clf; plot(xcorr(sig1, sig2), 'm');
hold on; plot(xcorr(sig1, sig3), 'b'); % 

% 127M 149F
sig1 = w(127).hst.spers / max(w(127).hst.spers); sig1 = sig1-mean(sig1);
sig2 = w(149).hst.spers / max(w(149).hst.spers); sig2 = sig2-mean(sig2);
hold on; plot(xcorr(sig1, sig2), 'k');

figure(7); clf; 
    subplot(211); hold on; 
    plot(w(10).hst.tim, w(10).hst.spers/max(w(10).hst.spers), 'b');
    plot(w(64).hst.tim, w(64).hst.spers/max(w(64).hst.spers), 'm')
    subplot(212); hold on; 
    plot(w(127).hst.tim, w(127).hst.spers, 'b');
    plot(w(149).hst.tim, w(149).hst.spers, 'm')

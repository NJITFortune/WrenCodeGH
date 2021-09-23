function out = ClickyClicky(dataPath, pairname)
% Usage: out = ClickyClicky(dataPath, pairname)
% Where dataPath might be '/Users/daynf/Documents/WrenData/postBirds-2/'
% and pairname might be 'post'
%
% relies on dualclicsB.m, slicerB.m, ginputc.m, syldataB.m, the SAMIII directory 
% Can be used with DualClickQuick.m
%

%% Load data

% Get male wav data
% dataPath = '/Users/daynf/Documents/WrenData/postBirds-2/';

%[mfilename,mpathname]=uigetfile('*.wav', 'Select Male WAV file');
[mfilename,mpathname]=uigetfile([dataPath, '/', pairname, 'Male*.wav'], 'Select Male WAV file');
[MrawData,mFs] = audioread([mpathname mfilename]);

% Get female wav data
%[ffilename,fpathname]=uigetfile('*.wav', 'Select female WAV file');
% automatically find matching female file
    fpathname=strrep(mpathname,'Male','Female'); 
    ffilename=strrep(mfilename,'Male','Female');
[FrawData,fFs] = audioread([fpathname ffilename]);
fprintf('Loaded female data %s \n', ffilename);

if mFs ~= fFs
    fprintf('Male sample rate %i does not match female sample rate %i', mFs, fFs);
    return
end

    out.Fs = mFs;
    out.mFile = mfilename; 
    out.fFile = ffilename;
    out.pairname = pairname;
    
        s = strsplit(mfilename(1:length(mfilename)-4),'_');        
            timestamp = str2num(s{end});
            out.distance = str2num(regexprep(s{1}, '\D', ''));
            out.day = str2num(s{2});
            out.month = str2num(s{3});
            out.year = str2num(s{4});
            out.location = s{5};
            out.timestamp = timestamp;
            
    FN = [pairname, num2str(out.distance),'m',num2str(timestamp),'.mat'];
    tempFN = ['tmpfile', num2str(out.distance),'m',num2str(timestamp),'.mat'];

%% Filter and normalize raw recording data

[b,a] = butter(5, 240 / (out.Fs/2), 'high'); % Highpass
% [d,c] = butter(5, 10000 / (out.Fs/2), 'low'); % Lowpass

    out.maleMic = filtfilt(b,a,MrawData); % Filter
    out.maleMic = 0.99 * (out.maleMic / (max(abs(out.maleMic)))); % Normalize
    out.femMic = filtfilt(b,a,FrawData); % Filter
    out.femMic = 0.99 * (out.femMic / (max(abs(out.femMic)))); % Normalize

%% Run the dualclics script

[clicked_m, clicked_f] = dualclicsB(out.maleMic, out.femMic, out.Fs, 0);

if length(clicked_m) ~= length(clicked_f) % This should never ever happen
    fprintf('Male Mic syllable count %i does not match Female Mic syllable count %i', length(clicked_m), length(clicked_f));
    return
end

% Save temporary data in case you fuck up slicing
    cd(dataPath) 
    save(tempFN,'clicked_m','clicked_f')

%% Concatonate the data

for ff = length(clicked_f):-1:1 % For each female syllable
    
    out.fsyl(ff).sylen = clicked_f(ff).sylen; % Length of seconds of the syllable (not necessary because given by syltim(2)-syltim(1)
    out.fsyl(ff).syltim = clicked_f(ff).syltim; % Beginning and end of each syllable
    out.fsyl(ff).sylidx = clicked_f(ff).sylind; % Beginning and end of each syllable
    
    out.fsyl(ff).traceFreq = clicked_f(ff).trace_freq; % Beginning and end of each syllable
    out.fsyl(ff).traceTim = clicked_f(ff).trace_tim; % Beginning and end of each syllable
    out.fsyl(ff).trace.slopemean = clicked_f.trace_slopemean;
    out.fsyl(ff).trace.slopestd = clicked_f.trace_slopestd;
    out.fsyl(ff).trace.slopevar = clicked_f.trace_slopevar;
    out.fsyl(ff).trace.peakfreq = clicked_f.trace_peakf;
    
    out.fsyl(ff).fftPower = clicked_f(ff).spec_p; % fft over entire syllable (power values)
    out.fsyl(ff).fftFreqs = clicked_f(ff).spec_f; % frequency values that for fftPower
    out.fsyl(ff).fft.peakf = clicked_f.spec_peakf;
    out.fsyl(ff).fft.peakf = clicked_f.spec_peakf;
    out.fsyl(ff).fft.minf = clicked_f.spec_minf;
    out.fsyl(ff).fft.maxf = clicked_f.spec_maxf;
    out.fsyl(ff).fft.band = clicked_f.spec_band;
    
end % End of female section

for mm = length(clicked_m):-1:1 % For each male syllable
    
    out.msyl(mm).sylen = clicked_m(mm).sylen; % Length of seconds of the syllable (not necessary because given by syltim(2)-syltim(1)
    out.msyl(mm).syltim = clicked_m(mm).syltim; % Beginning and end of each syllable
    out.msyl(mm).sylidx = clicked_m(mm).sylind; % Beginning and end of each syllable
    
    out.msyl(mm).traceFreq = clicked_m(mm).trace_freq; % Beginning and end of each syllable
    out.msyl(mm).traceTim = clicked_m(mm).trace_tim; % Beginning and end of each syllable
    out.msyl(mm).trace.slopemean = clicked_m(mm).trace_slopemean;
    out.msyl(mm).trace.slopestd = clicked_m(mm).trace_slopestd;
    out.msyl(mm).trace.slopevar = clicked_m(mm).trace_slopevar;
    out.msyl(mm).trace.peakfreq = clicked_m(mm).trace_peakf;
    
    out.msyl(mm).mmtPower = clicked_m(mm).spec_p; % mmt over entire syllable (power values)
    out.msyl(mm).mmtFreqs = clicked_m(mm).spec_f; % frequency values that for mmtPower
    out.msyl(mm).fft.peakf = clicked_m(mm).spec_peakf;
    out.msyl(mm).fft.peakf = clicked_m(mm).spec_peakf;
    out.msyl(mm).fft.minf = clicked_m(mm).spec_minf;
    out.msyl(mm).fft.maxf = clicked_m(mm).spec_maxf;
    out.msyl(mm).fft.band = clicked_m(mm).spec_band;
    
end % End of male section

save(tempFN, 'out', 'clicked_m', 'clicked_f') %save temp.mat clicked_m clicked_f out % JUST IN CASE

%% Slicer time - identify syllables and who sang what

% Have the user sort syllables on the basis of the traces
hopeandpray = 1;

while hopeandpray ~=0
    
    msyls = slicerB(clicked_m);
    fsyls = slicerB(clicked_f);
    
    if (length(fsyls) ~= length(msyls))
        fprintf('Did not make the same number of syllable types, clickturd - you have to do this again. \n ');
    end
    
    hopeandpray = abs(length(fsyls) - length(msyls));
end

% Quality control - make sure syllable sequence is the same between
% microphones. This is an unlikely event.

if sum(find([msyls.num] ~= [fsyls.num])) ~= 0
    fprintf('Syllable order does not match between microphones.\n');
    
    currMaleSequence = zeros(1,max([msyls.num]));
    for jk = length(msyls):-1:1
        currMaleSequence(msyls(jk).num) = jk;
    end
    
    currFemaleSequence = zeros(1,max([fsyls.num]));
    for jk = length(fsyls):-1:1
        currFemaleSequence(fsyls(jk).num) = jk;
    end
    
    % Make a plot
    figure(27); clf; 
    subplot(211); title(FN,'Interpreter','none'); specgram(out.maleMic, 1024, out.Fs); ylim([200 5200]);
    caxis([-10 40]); colormap(flipud(gray));
    hold on;
    for j=1:length(out.msyl)
        plot([out.msyl(j).syltim(1) out.msyl(j).syltim(1)], [500 4500], 'g', 'LineWidth', 3);
        plot([out.msyl(j).syltim(2) out.msyl(j).syltim(2)], [500 4500], 'm', 'LineWidth', 3);
        text(out.msyl(j).syltim(1)+0.1, 4000, num2str(currMaleSequence(j)));
    end
    subplot(212); specgram(out.femMic, 1024, out.Fs); ylim([200 5200]);
    caxis([-10 40]); colormap(flipud(gray));
    hold on;
    for j=1:length(out.fsyl)
        plot([out.fsyl(j).syltim(1) out.fsyl(j).syltim(1)], [500 4500], 'g', 'LineWidth', 3);
        plot([out.fsyl(j).syltim(2) out.fsyl(j).syltim(2)], [500 4500], 'm', 'LineWidth', 3);
        text(out.fsyl(j).syltim(1)+0.1, 4000, num2str(currFemaleSequence(j)));
    end
    
%    neworder = inputdlg('Enter proper order: ', 'Error Will Robinson');
    prompt = 'Enter proper order: '; dlgtitle = 'Error Will Robinson'; definput = {' '}; opts.WindowStyle = 'normal';
    neworder = inputdlg(prompt, dlgtitle, [1 50], definput, opts);
    maxsyl = max(str2num(neworder{1}));
    for jj = maxsyl:-1:1
        msyls(jj).num = find(str2num(neworder{1}) == jj);
        fsyls(jj).num = find(str2num(neworder{1}) == jj);
    end
    
    close(27);
end


%% Assign sex

maleMicAmp = 0; femMicAmp = 0;

for pp = 1:length(msyls)
    
    for jj = 1:length(msyls(pp).num)
        maleMicAmp = maleMicAmp + sum(abs(out.maleMic(out.msyl(msyls(pp).num(jj)).sylidx(1):out.msyl(msyls(pp).num(jj)).sylidx(2))));
        femMicAmp = femMicAmp + sum(abs(out.femMic(out.fsyl(fsyls(pp).num(jj)).sylidx(1):out.fsyl(fsyls(pp).num(jj)).sylidx(2))));
    end
    
    if maleMicAmp > femMicAmp % This is a male syllable
        for q = 1:length([msyls(pp).num])
            out.msyl(msyls(pp).num(q)).sexsyltype = pp;
            out.msyl(fsyls(pp).num(q)).sex = 'M';
            out.fsyl(fsyls(pp).num(q)).sexsyltype = pp;
            out.fsyl(fsyls(pp).num(q)).sex = 'M';
        end
    end
    if maleMicAmp < femMicAmp % This is a female syllable
        for q = 1:length([msyls(pp).num])
            out.msyl(msyls(pp).num(q)).sexsyltype = 50+pp;
            out.msyl(fsyls(pp).num(q)).sex = 'F';
            out.fsyl(fsyls(pp).num(q)).sexsyltype = 50+pp;
            out.fsyl(fsyls(pp).num(q)).sex = 'F';
        end
    end
    
end

currMaleSequence = zeros(1,max([msyls.num]));
for jk = length(msyls):-1:1
    currMaleSequence(msyls(jk).num) = jk;
end
currFemaleSequence = zeros(1,max([fsyls.num]));
for jk = length(fsyls):-1:1
    currFemaleSequence(fsyls(jk).num) = jk;
end

scr_siz = get(0,'ScreenSize') ;
figure(3); clf; 
    set(gcf, 'Position', [50 50 scr_siz(3)-200 scr_siz(4)-200])

subplot(211); specgram(out.maleMic, 1024, out.Fs); ylim([200 5200]);
caxis([-10 40]); colormap(flipud(gray));
hold on;
for j=1:length(out.msyl)
    plot([out.msyl(j).syltim(1) out.msyl(j).syltim(1)], [500 4500], 'g', 'LineWidth', 2);
    plot([out.msyl(j).syltim(2) out.msyl(j).syltim(2)], [500 4500], 'r', 'LineWidth', 2);
    if out.msyl(j).sex == 'M'
        text(out.msyl(j).syltim(1)+0.1, 4000, num2str(currMaleSequence(j)), 'Color','cyan','FontSize', 14);
    elseif out.msyl(j).sex == 'F'
        text(out.msyl(j).syltim(1)+0.1, 4000, num2str(currMaleSequence(j)), 'Color','magenta','FontSize', 14);
    end
end
text(0.5, 1000, 'Male Microphone');

subplot(212); specgram(out.femMic, 1024, out.Fs); ylim([200 5200]);
caxis([-10 40]); colormap(flipud(gray));
hold on;
for j=1:length(out.fsyl)
    plot([out.fsyl(j).syltim(1) out.fsyl(j).syltim(1)], [500 4500], 'g', 'LineWidth', 2);
    plot([out.fsyl(j).syltim(2) out.fsyl(j).syltim(2)], [500 4500], 'r', 'LineWidth', 2);
    if out.fsyl(j).sex == 'M'
        text(out.fsyl(j).syltim(1)+0.1, 4000, num2str(currFemaleSequence(j)), 'Color','cyan','FontSize', 14);
    elseif out.fsyl(j).sex == 'F'
        text(out.fsyl(j).syltim(1)+0.1, 4000, num2str(currFemaleSequence(j)), 'Color','magenta','FontSize', 14);
    end
end
text(0.5, 1000, 'Female Microphone');

%% Final opportunity to fix the data

fprintf('Last chance to fix numbers and sex! \n');
    prompt = 'Enter 1 to fix.'; dlgtitle = 'IS THE SEQUENCE CORRECT?'; definput = {' '}; opts.WindowStyle = 'normal';
    yn = inputdlg(prompt, dlgtitle, [1,50], definput, opts);
    % yn = inputdlg('Enter 1 to fix.', 'IS THE SEQUENCE CORRECT?');
if ~isempty(yn{1})
    if str2num(yn{1}) == 1
        sexstring = [];
         for j=1:length(out.fsyl)
             sexstring = [sexstring out.fsyl(j).sexsyltype];
         end
        
        fprintf('Current syllables: %s \n', num2str(sexstring)) ;
        prompt = 'num2str(sexstring)'; dlgtitle = 'Enter new numbers!'; definput = {' '}; opts.WindowStyle = 'normal';
        neworder = inputdlg(prompt, dlgtitle, [1 50], definput, opts);
        %tmpnewstr = inputdlg(num2str(sexstring),'Enter new numbers!');
        newstr = str2num(tmpnewstr{1});
        
        while length(newstr) ~= length(out.fsyl)
            fprintf('WRONG NUMBER OF SYLLABLES, fucktard.\n');
            tmpnewstr = inputdlg('Enter new numbers: ', 'Numbers', [1,50], 'WindowStyle', 'normal');
            newstr = str2num(tmpnewstr{1});
        end
        
        for j = 1:length(newstr)
            out.fsyl(j).sexsyltype = newstr(j);
            out.msyl(j).sexsyltype = newstr(j);
            if newstr(j) > 50
                out.fsyl(j).sex = 'F';
                out.msyl(j).sex = 'F';
            end
            if newstr(j) < 50
                out.fsyl(j).sex = 'M';
                out.msyl(j).sex = 'M';
            end
        end
    end
end
%% One final question

    vv = inputdlg('One final question... 1 vision, 2 no vision', 'Vision?', [1,5], 'WindowStyle', 'normal');
    out.vision = str2num(vv{1});
    

%% Cleanup

save(FN,'out')
fprintf('File %s saved.\n', FN)

% Move files to completed dir
if exist([mpathname,'clicked'],'dir') == 0 
    mkdir([mpathname,'clicked'])
end

if exist([fpathname,'clicked'],'dir') == 0 
    mkdir([fpathname,'clicked'])
end

    movefile([mpathname,mfilename],[mpathname,'clicked/',mfilename])
    movefile([fpathname,ffilename],[fpathname,'clicked/',ffilename])
fprintf('wav files moved to %s.\n', [mpathname,'clicked/',mfilename]);
    




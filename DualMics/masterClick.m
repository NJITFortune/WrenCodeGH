function out = masterClick(dist, vis, loc, day, mon, yr)
% Usage: out = masterClick(dist, vis, loc, day, mon, yr)

%% Prepartions

    out.dist = dist; out.vis = vis; out.loc = loc;
    out.day = day; out.mon = mon; out.yr = yr;

% Get male wav data
    [filename,pathname]=uigetfile('*.wav', 'Select Male WAV file');
        [MrawData,mFs] = audioread([pathname filename]);
% Get female wav data
    [filename,pathname]=uigetfile('*.wav', 'Select female WAV file');
        [FrawData,fFs] = audioread([pathname filename]);

if mFs ~= fFs
   fprintf('Male sample rate %i does not match female sample rate %i', mFs, fFs); 
   return
end

    out.Fs = mFs;    
    
% Filter and normalize raw recording data

    [b,a] = butter(5, 240 / (out.Fs/2), 'high'); % Highpass
    % [d,c] = butter(5, 10000 / (out.Fs/2), 'low'); % Lowpass
    
    out.maleMic = filtfilt(b,a,MrawData); % Filter
        out.maleMic = 0.99 * (out.maleMic / (max(abs(out.maleMic)))); % Normalize
    out.femaleMic = filtfilt(b,a,FrawData); % Filter
        out.femaleMic = 0.99 * (out.femaleMic / (max(abs(out.femaleMic)))); % Normalize
        
%% Run the dualclics script
    [clicked_m, clicked_f] = dualclics(out.maleMic, out.femaleMic, out.Fs, 0);
    
    save temp.mat clicked_m clicked_f out % JUST IN CASE

if length(clicked_m) ~= length(clicked_f) % This should never ever happen
   fprintf('Male syllable count %i does not match female syllable count %i', length(clicked_m), length(clicked_f)); 
   return
end
        
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

    save temp.mat clicked_m clicked_f out % JUST IN CASE

%% Slicer time - identify syllables and who sang what

% Have the user sort syllables on the basis of the traces
hopeandpray = 1;
while hopeandpray ~=0
    
        msyls = slicer(clicked_m);
        fsyls = slicer(clicked_f);

        if (length(fsyls) ~= length(msyls))
            fprintf('Did not make the same number of syllable types, clickturd - you have to do this again. \n ');
        end

        hopeandpray = length(fsyls) - length(msyls);
end

% Quality control - make sure syllable sequency is the same between
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
    figure; clf;
    subplot(211); specgram(out.maleMic, 1024, out.Fs); ylim([200 5200]); 
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
    
neworder = input('Enter proper order: ');
    maxsyl = max(neworder);
    for jj = maxsyl:-1:1
        msyls(jj).num = find(neworder == jj);
        fsyls(jj).num = find(neworder == jj);
    end

end

%% Assign sex

maleMicAmp = 0; femMicAmp = 0;

    for pp = 1:length(msyls)
     
        for jj = 1:length(msyls(pp).num)
          maleMicAmp = maleMicAmp + sum(abs(out.   out.msyl(msyls(pp).num(jj)).sylidx(1): ))
      
        for q = 1:length([msyls(p).num])
            out.msyl(msyls(p).num(q)).sexsyltype = p;
        end
        for q = 1:length([fsyls(p).num])
            out.fsyl(fsyls(p).num(q)).sexsyltype = p;
        end
  end



%% Fix the data


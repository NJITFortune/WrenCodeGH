function out = distantwren
% Usage out = distantwren(f, m, Fs)
% f is Female signal, m is Male signal, Fs is samplerate in Hz
% NOW USES UIGETFILE

%% Get data, make preparations, and plot

% Female
[femfile, fempath, ~] = uigetfile('*.wav', 'Pick the FEMALE wav file');
ff = fullfile(fempath, femfile);
[f, Fs] = audioread(ff);

    if length(f(1,:)) > 1;
        figure(1); clf; 
        subplot(211); plot(f(:,1), 'm'); ylim([-1.1 1.1]);
        subplot(212); plot(f(:,2), 'm'); ylim([-1.1 1.1]);
    
        pickchan = input('Which channel is better? ');    
        f = f(:,pickchan);  f = f / max(abs(f)); 
        figure(1); clf; plot(f, 'm');
    end;

% Male
[malfile, malpath, ~] = uigetfile('*.wav', 'Pick the MALE wav file');
mm = fullfile(malpath, malfile);
[m, Fs] = audioread(mm);

    if length(m(1,:)) > 1;
        figure(1); clf; 
        subplot(211); plot(m(:,1), 'b'); ylim([-1.1 1.1]);
        subplot(212); plot(m(:,2), 'b'); ylim([-1.1 1.1]);
    
        pickchan = input('Which channel is better? ');    
        m = m(:,pickchan);  m = m / max(abs(m));  
        figure(1); clf; plot(m, 'b');
    end;

% Build the Hilbert functions
[b,a] = butter(3, 2*(500/Fs), 'high');
[d,c] = butter(3, 2*(20/Fs), 'low');
hf = filtfilt(d,c,abs(hilbert(filtfilt(b,a,f))));
    hf = hf - min(hf); hf = hf / max(hf); 
hm = filtfilt(d,c,abs(hilbert(filtfilt(b,a,m))));
    hm = hm - min(hm); hm = hm / max(hm); 

% Make our time base
ftim = 1/Fs:1/Fs:length(f)/Fs;
mtim = 1/Fs:1/Fs:length(m)/Fs;

% Plot original signals
figure(2); clf;
%ax(1) = subplot(311); specgram(f,2048,Fs); colormap('HOT'); ylim([0 6000]);
%hold on; plot(ftim(1:10:end), hf(1:10:end), 'm', 'LineWidth', 2); hold off;
%ax(2) = subplot(312); specgram(m,2048,Fs); colormap('HOT'); ylim([0 6000]);
%hold on; plot(mtim(1:10:end), hm(1:10:end), 'b', 'LineWidth', 2); hold off;
%ax(3) = subplot(313); plot(ftim(1:10:end), hf(1:10:end), 'm', 'LineWidth', 2);
%subplot(313); 
plot(ftim(1:10:end), hf(1:10:end), 'm', 'LineWidth', 2);
hold on; plot(mtim(1:10:end), hm(1:10:end), 'b', 'LineWidth', 2); 
%linkaxes(ax,'x');

%% Extract Duet Segments
% We need to process each duet alone - sample may have more than one
numduets = input('How many Duets do you see in this recording? ');

% If there is only one duet, we use the whole sample
    out(1).male = m; out(1).mtim = mtim;
    out(1).female = f; out(1).ftim = ftim;
    out(1).hm = hm; out(1).hf = hf;

% If there are more than one duet, we have the user click
    if numduets > 1; % If more than one duet, we need to click
        sprintf('Click start and end of each duet in bottom panel. ');
        for j = 2:2:numduets*2;
            figure(2); hold on;
            [xx(j-1), ~] = ginput(1); plot([xx(j-1) xx(j-1)], [0 1], 'g', 'LineWidth', 2);
            [xx(j), ~] = ginput(1); plot([xx(j) xx(j)], [0 1], 'r', 'LineWidth', 2);
            ftt = find(ftim > xx(j-1) & ftim < xx(j));
            mtt = find(mtim > xx(j-1) & mtim < xx(j));
            out(j/2).mtim = mtim(mtt); out(j/2).male = m(mtt);
            out(j/2).ftim = ftim(ftt); out(j/2).female = f(ftt);
            out(j/2).hf = hf(ftt); out(j/2).hm = hm(mtt);
        end;
    end;

    pause(1);

%% Process Duets
for i = 1:numduets;

    % Adjust female envelope
    figure(1); clf;
    semilogy(out(i).ftim, out(i).hf, 'm');
    xlim([min(out(i).ftim) max(out(i).ftim)]);
    axf = axis;
    disp('Click just below the baseline.');
    [~, axf(3)] = ginput(1);
    axis(axf); axisf(:,i) = axf;
    
    pause(1);
    
    % Adjust male envelope
    semilogy(out(i).mtim, out(i).hm, 'b');
    xlim([min(out(i).mtim) max(out(i).mtim)]);
    axm = axis;
    disp('Click just below the baseline.');
    [~, axm(3)] = ginput(1);
    axis(axm); axism(:,i) = axm;
    
    pause(1);

% Get the first syllable    
   figure(1); clf; % Female recording   
    subplot(211); semilogy(out(i).ftim, out(i).hf, 'm'); axis(axisf(:,i));
    subplot(212); specgram(out(i).female, 2048, Fs, [], 2000); colormap('HOT'); ylim([500 6000]); caxis([-30 20]);
        
    disp('Click precisely at the base of the first syllable in this recording near the female.');
    figure(1); subplot(211);
    [fx(i), ~] = ginput(1);
    
  figure(2); clf; % Male recording   
    subplot(211); semilogy(out(i).mtim, out(i).hm, 'b'); axis(axism(:,i));
    subplot(212); specgram(out(i).male, 2048, Fs, [], 2000); colormap('HOT'); ylim([500 6000]); caxis([-30 20]);

  disp('Click precisely at the base of the first syllable in this recording near the male.');
    figure(2); subplot(211);
    [mx(i), ~] = ginput(1);    
    
% Pull out the sample for this duet    
    ttf{i} = find(out(i).ftim >= (fx(i)-1)); % Get 1 second before the first click
    ttm{i} = find(out(i).mtim >= (mx(i)-1)); % Get 1 second before the first click
    
    minlen = min([length(ttf{i}) length(ttm{i})]); % Go to the end of the sample
    ttf{i} = ttf{i}(1:minlen);
    ttm{i} = ttm{i}(1:minlen);


fm = input('Was the first syllable a male (0) or a female (1) syllable?? ');

    
    if fm == 1; % Female first
        out(i).fc.fem = out(i).female(ttf{i});
        out(i).fc.mal = out(i).male(ttm{i});
        out(i).fc.ftim = out(i).ftim(ttf{i});
        out(i).fc.hf = out(i).hf(ttf{i});
        out(i).fc.hm = out(i).hm(ttm{i});
        out(i).fc.atim = 1/Fs:1/Fs:length(ttf{i})/Fs;
    end;
    if fm == 0; % Male first
        out(i).mc.fem = out(i).female(ttf{i});
        out(i).mc.mal = out(i).male(ttm{i});
        out(i).mc.mtim = out(i).mtim(ttm{i});
        out(i).mc.hm = out(i).hm(ttm{i});
        out(i).mc.hf = out(i).hf(ttf{i});
        out(i).mc.atim = 1/Fs:1/Fs:length(ttf{i})/Fs;
    end;
    
    if fm == 1; % Female was first, now get the male
    % Now we need to get the male start            
    
    figure(2); clf; % First round, male recording
        subplot(211); semilogy(out(i).mtim, out(i).hm, 'b'); axis(axism(:,i));
        subplot(212); specgram(out(i).male, 2048, Fs, [], 2000); colormap('HOT'); ylim([500 6000]); caxis([-30 20]);

        disp('Click precisely at the base of the first MALE syllable in this recording.');
        figure(2); subplot(211);
        [nmx(i), ~] = ginput(1);
        
    figure(1); clf; % Second round, fenale recording   
        subplot(211); semilogy(out(i).ftim, out(i).hf, 'm'); axis(axisf(:,i));
        subplot(212); specgram(out(i).female, 2048, Fs, [], 2000); colormap('HOT'); ylim([500 6000]); caxis([-30 20]);
        
        disp('Click precisely at the base of the first MALE syllable in this recording.');
        figure(1); subplot(211);
        [nfx(i), ~] = ginput(1);

        
    bactine = 1 + min([nmx(i)-mx(i), nfx(i)-fx(i)]);
    tttf{i} = find(out(i).ftim >= (nfx(i)-bactine));
    tttm{i} = find(out(i).mtim >= (nmx(i)-bactine));
    
    minlen = min([length(tttf{i}) length(tttm{i})]);
    tttf{i} = tttf{i}(1:minlen);
    tttm{i} = tttm{i}(1:minlen);
    
    out(i).mc.fem = out(i).female(tttf{i});
    out(i).mc.mal = out(i).male(tttm{i});
    out(i).mc.mtim = out(i).mtim(tttm{i});
    out(i).mc.hm = out(i).hm(tttm{i});
    out(i).mc.hf = out(i).hf(tttf{i});
    out(i).mc.atim = 1/Fs:1/Fs:length(tttm{i})/Fs;
    end;

    if fm == 0; % Male was first, now get the female
    % Now we need to get the female start            
    figure(1); clf; % First round, female recording   
        subplot(211); semilogy(out(i).ftim, out(i).hf, 'm'); axis(axisf(:,i));
        subplot(212); specgram(out(i).female, 2048, Fs, [], 2000); colormap('HOT'); ylim([500 6000]); caxis([-30 20]);
        
        disp('Click precisely at the base of the first FEMALE syllable in this recording.');
        figure(1); subplot(211);
        [nfx(i), ~] = ginput(1);
    
    figure(2); clf; % Second round, male recording
        subplot(211); semilogy(out(i).mtim, out(i).hm, 'b'); axis(axism(:,i));
        subplot(212); specgram(out(i).male, 2048, Fs, [], 2000); colormap('HOT'); ylim([500 6000]); caxis([-30 20]);

        disp('Click precisely at the base of the first FEMALE syllable in this recording.');
        figure(2); subplot(211);
        [nmx(i), ~] = ginput(1);

    bactine = 1 + min([nmx(i)-mx(i), nfx(i)-fx(i)]);
    tttf{i} = find(out(i).ftim >= (nfx(i)-bactine));
    tttm{i} = find(out(i).mtim >= (nmx(i)-bactine));
    
    minlen = min([length(tttf{i}) length(tttm{i})]);
    tttf{i} = tttf{i}(1:minlen);
    tttm{i} = tttm{i}(1:minlen);
    
    out(i).fc.fem = out(i).female(tttf{i});
    out(i).fc.mal = out(i).male(tttm{i});
    out(i).fc.ftim = out(i).mtim(tttm{i});
    out(i).fc.hf = out(i).hf(tttf{i});
    out(i).fc.hm = out(i).hm(tttm{i});
    out(i).fc.atim = 1/Fs:1/Fs:length(tttf{i})/Fs;
    end;

    
    
    
% Fabulous Figure 3 shows all of the processed and raw data. 
% 411 is Female specgram, 412 is Female-aligned envelope
% 413 us Male-aligned envelope, and 414 is Male specgram
figure(1); clf;  
    subplot(311); % Female specgram
    specgram(out(i).fc.fem,2048, Fs, [], 2000); colormap('HOT'); ylim([500 6000]); caxis([-30 20]);        
    subplot(312); % Female specgram (male recording)
    specgram(out(i).fc.mal,2048, Fs, [], 2000); colormap('HOT'); ylim([500 6000]); caxis([-30 20]);        
    subplot(313); % Female-centric diff envelope plot     
            yforplot = out(i).fc.hf - out(i).fc.hm;
            out(i).fc.pfem = find(yforplot > 0.02);
            fz = zeros(1, length(out(i).fc.atim));
            fz(out(i).fc.pfem) = 1;
            area(out(i).fc.atim, fz, 'FaceColor', 'm', 'EdgeColor', 'm'); 
            hold on;
            plot(out(i).fc.atim, yforplot, 'k');
            area(out(i).fc.atim(out(i).fc.pfem), yforplot(out(i).fc.pfem), 'FaceColor', 'w', 'EdgeColor', 'm');
            tmp = axis; axis([out(i).fc.atim(1) out(i).fc.atim(end) tmp(3) tmp(4)]); 
figure(2); clf;
    subplot(311); % Male specgram
    specgram(out(i).mc.mal,2048, Fs, [], 2000); colormap('HOT'); ylim([500 6000]); caxis([-30 20]);        
    subplot(312); % Male specgram (female recording)
    specgram(out(i).mc.fem,2048, Fs, [], 2000); colormap('HOT'); ylim([500 6000]); caxis([-30 20]);        
    subplot(313); % Male-centric diff envelope plot
            yforplot = out(i).mc.hm - out(i).mc.hf;
            out(i).mc.pmal = find(yforplot > 0.02);
            mz = zeros(1, length(out(i).mc.atim));
            mz(out(i).mc.pmal) = 1;
            area(out(i).mc.atim, mz, 'FaceColor', 'b', 'EdgeColor', 'b'); 
            hold on;
            plot(out(i).mc.atim, yforplot, 'k');
            area(out(i).mc.atim(out(i).mc.pmal), yforplot(out(i).mc.pmal), 'FaceColor', 'w', 'EdgeColor', 'b');
            tmp = axis; axis([out(i).mc.atim(1) out(i).mc.atim(end) tmp(3) tmp(4)]); 
            
pause(10);    
aa = input('Next? ');


% Process the data





end;


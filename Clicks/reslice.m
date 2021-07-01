function reslice
% load the temp.mat file
fn=uigetfile('Select TEMP file','.mat');
load(fn)
savefn=strrep(fn,'TEMP','');
timestamp=fn(11:length(fn)-4);
%% Slicer time - identify syllables and who sang what

% Have the user sort syllables on the basis of the traces
hopeandpray = 1;
while hopeandpray ~=0
    
        msyls = slicer(clicked_m);
        fsyls = slicer(clicked_f);

        if (length(fsyls) ~= length(msyls))
            fprintf('Did not make the same number of syllable types, clickturd - you have to do this again. \n ');
        end

        hopeandpray = abs(length(fsyls) - length(msyls));
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
    figure(27); clf;
    subplot(211); title(FN,'Interpreter','none');specgram(out.maleMic, 1024, out.Fs); ylim([200 5200]); 
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

figure(3); clf;
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
save(savefn,'out');

postPath='/Users/daynf/Documents/WrenData/postBirds-2/';
mpathname=[postPath,'Male-',num2str(out.dist),'m/'];
fpathname=[postPath,'Female-',num2str(out.dist),'m/'];


if exist([mpathname,'clicked']) == 0 
    mkdir([mpathname,'clicked'])
end

if exist([fpathname,'clicked']) == 0 
    mkdir([fpathname,'clicked'])
end

mfilename=['postMale',num2str(out.dist),'m_August_24_2016_',num2str(timestamp),'.wav'];
ffilename=['postFemale',num2str(out.dist),'m_August_24_2016_',num2str(timestamp),'.wav'];


%move files to completed dir
movefile([mpathname,mfilename],[mpathname,'clicked/',mfilename])
movefile([fpathname,ffilename],[fpathname,'clicked/',ffilename])
disp('Files moved')

function DualClickQuick(out)

% Speed of sound (milliseconds per meter)
ss = 1 / 343;

% Make a plot
    figure(27); clf; 

%% Male Microphone    
    ax(1) = subplot(312); title(out.pairname,'Interpreter','none'); 
    specgram(out.maleMic, 1024, out.Fs); ylim([200 6000]);
    caxis([-10 40]); colormap(flipud(gray));
    hold on;
    for j=1:length(out.msyl)
        plot([out.msyl(j).syltim(1) out.msyl(j).syltim(1)], [500 5500], 'g', 'LineWidth', 2);
        plot([out.msyl(j).syltim(2) out.msyl(j).syltim(2)], [500 5500], 'r', 'LineWidth', 2);
        if out.msyl(j).sexsyltype > 50 
            text(out.msyl(j).syltim(1)+0.05, 4000, num2str(out.msyl(j).sexsyltype-50), 'Color', 'magenta');
        else
            text(out.msyl(j).syltim(1)+0.05, 4000, num2str(out.msyl(j).sexsyltype), 'Color', 'cyan');            
        end
    end
    text(0.5, 500, 'Male Microphone');

%% Female Microphone
    ax(2) = subplot(311); specgram(out.femMic, 1024, out.Fs); ylim([200 5200]);
    caxis([-10 40]); colormap(flipud(gray));
    hold on;
    for j=1:length(out.fsyl)
        plot([out.fsyl(j).syltim(1) out.fsyl(j).syltim(1)], [500 5500], 'g', 'LineWidth', 2);
        plot([out.fsyl(j).syltim(2) out.fsyl(j).syltim(2)], [500 5500], 'r', 'LineWidth', 2);
        if out.fsyl(j).sexsyltype > 50 
            text(out.fsyl(j).syltim(1)+0.05, 4000, num2str(out.fsyl(j).sexsyltype-50), 'Color', 'magenta');
        else
            text(out.fsyl(j).syltim(1)+0.05, 4000, num2str(out.fsyl(j).sexsyltype), 'Color', 'cyan');            
        end
    end
    text(0.5, 500, 'Female Microphone');

linkaxes(ax, 'xy');

%% ISI Plot
for j = length(out.fsyl):-1:2

    fISI(j-1) = out.fsyl(j).syltim(1) - out.fsyl(j-1).syltim(2);
    mISI(j-1) = out.msyl(j).syltim(1) - out.msyl(j-1).syltim(2);

end

subplot(313); hold on;

    plot(fISI*1000, 'mo-');  
    plot(mISI*1000, 'bo-');
    isidiffs = abs(mISI - fISI)*1000;
    plot(isidiffs, 'k.-', 'MarkerSize', 16);
    xlabel('Interval number');
    ylabel('ISI msec')



    ISImsg = ['Mean ISI diff: ' num2str(mean(isidiffs)) ' Distance delay: ' num2str(2 * ss * out.distance * 1000)];

    text(1.5, max([fISI*1000 mISI*1000])*0.9, ISImsg, 'Color', 'k', 'FontSize', 18);
    



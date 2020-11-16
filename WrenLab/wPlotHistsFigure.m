function wPlotHistsFigure(in)
% This script makes the "raw" physiology plots for the manuscript.
% The plots were combined and cleaned up in Illustrator / Designer for
% publication.

% We are using the examples M17 (indices 1-2) and j160815 (indices 7-8).
% These are the prettiest examples! 
% Figure 27 has the two spectrograms

% You can modify this script to plot any of the data by changing the
% idx values below.

% [msolosyls, mduetsyls, fsolosyls, fduetsyls, Cspon, Aspon] = wData;


%% PLOT THE FIRST EXAMPLE, M17

limits = [0.2 6.0];
sWIN = 50;

% repnum = min([length(in(1).Aspikes), length(in(2).Aspikes)] );
repnum = 30;

idx = 1; % Male M17  %% IDX values for males are ODD

figure('Name', 'Male M17 idx=1 AWAKE'); clf; set(gcf,'renderer','Painters');
    MakePlot(in(idx), limits, sWIN, 4, 1);
figure('Name', 'Male M17 idx=1 URETHANE' ); clf; set(gcf,'renderer','Painters');
    MakePlot(in(idx), limits, sWIN, repnum, 2);

idx = 2; % Female M17 %% IDX values for females are EVEN, follow male idx

figure('Name', 'Female M17 idx=2 AWAKE'); clf; set(gcf,'renderer','Painters');
    MakePlot(in(idx), limits, sWIN, 4, 1);
figure('Name', 'Female M17 idx=2 URETHANE'); clf; set(gcf,'renderer','Painters');
    MakePlot(in(idx), limits, sWIN, repnum, 2);

    
figure(27); clf; 
    subplot(211); 
    specgram(in(idx).duet(in(idx).tim > limits(1) & in(idx).tim < limits(2)), 512, in(idx).Fs, [], 500);
    colormap(flipud(gray));
    ylim([500 4500]); caxis([-35 30]);
    
%% PLOT THE SECOND EXAMPLE, j160815

limits = [0 4.5];
sWIN = 50;

% repnum = min([length(in(7).Aspikes), length(in(8).Aspikes)] );
repnum = 30;

idx = 7; % Male j160815  %% IDX values for males are ODD

figure('Name', 'Male j160815 idx=7 AWAKE'); clf; set(gcf,'renderer','Painters');
    MakePlot(in(idx), limits, sWIN, 4, 1);
figure('Name', 'Male j160815 idx=7 URETHANE'); clf; set(gcf,'renderer','Painters');
    MakePlot(in(idx), limits, sWIN, repnum, 2);

idx = 8; % Female j160815 %% IDX values for females are EVEN, follow male idx

figure('Name', 'Female j160815 idx=8 AWAKE'); clf; set(gcf,'renderer','Painters');
    MakePlot(in(idx), limits, sWIN, 4, 1);
figure('Name', 'Female j160815 idx=8 URETHANE'); clf; set(gcf,'renderer','Painters');
    MakePlot(in(idx), limits, sWIN, repnum, 2);


figure(27); set(gcf,'renderer','Painters');
    subplot(212); 
    specgram(in(idx).duet(in(idx).tim > limits(1) & in(idx).tim < limits(2)), 512, in(idx).Fs, [], 500);
    colormap(flipud(gray));
    ylim([500 4500]); caxis([-20 30]);
    drawnow;

%% Plotting Function

    function MakePlot(d, w, b, r, cORa) 
        
        signal = d.duet;
        tim = d.tim;
        set(gcf, 'Color', [1,1,1]); % This sets the background to white
        
        if cORa == 1; spikes = d.Cspikes; end
        if cORa == 2 
            spikes = d.Aspikes; 
        end
        
        % Raster Plot
        axs(1)=subplot(3,1,1);
        PLOTraster(spikes, r);
        xlim([w(1) w(2)]);
        h = gca; set(h, 'XTickLabel', [], 'Ydir', 'reverse', 'YAxisLocation', 'left');
        
        % Histogram Plot
        axs(2)=subplot(3,1,2);
        swpsthdata = bs_swPSTH(spikes, w, b, 0);
        hold on; 
        plot(swpsthdata.tim, swpsthdata.spers, 'Color', '[0.4940, 0.1840, 0.5560]', 'LineWidth', 3);
        xlim([w(1) w(2)]);
        h = gca; set(h, 'XTickLabel', [], 'Box', 'off');

        % Oscillogram plot
        axs(3)=subplot(3,1,3);
        tt = find(tim > w(1) & tim < w(2));
        plot(tim(tt),signal(tt), 'k');
        h = gca; set(h, 'Ycolor', [1,1,1], 'Box', 'off'); 
        f = axis;
        axis([w(1) w(2) f(3) f(4)]);

        linkaxes(axs,'x');
        

        for j=1:3 % For each subplot
            subplot(3,1,j); hold on;
            xAx = axis;

            for k=1:length(d.syl)
                if ~isempty(d.syl(k).tim)
                plot([d.syl(k).tim(1), d.syl(k).tim(1)], [xAx(3), xAx(4)], 'g-');

                if d.sylsex(k) == 1 % This is a male syllable
                    text(d.syl(k).tim(1), xAx(4)-(0.10*xAx(4)), num2str(k), 'Color', 'b', 'FontSize',14);
                end
                if d.sylsex(k) == 2 % This is a female syllable
                    text(d.syl(k).tim(1), xAx(4)-(0.10*xAx(4)), num2str(k), 'Color', 'm', 'FontSize',14);
                end

                plot([d.syl(k).tim(2), d.syl(k).tim(2)], [xAx(3), xAx(4)], 'r-');        
                end        
            end
        end

        drawnow; pause(0.5);
        
    end % End of plot function

function PLOTraster(spiketimes, reps)

hold on;

for tt = 1:reps
   for ss = 1:length( spiketimes{tt} )
      plot( ones(1,2).*spiketimes{tt}(ss), [tt-0.3 tt+0.3], 'k', 'LineWidth', 0.5);
   end
end

ax = axis;
axis( [ax(1:2) 0 reps+1] )

end


end % End of function



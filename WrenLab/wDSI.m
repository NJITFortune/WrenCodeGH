function foo = wDSI(wren, pdnum)
% This script calculates the "Duet Selectivity Index" (which is a poor name). 
% The calculation is identical to the "Direction Selectivity Index" that
% provides a ratio of responses to moving stimuli with opposite sign
% velocities. In our case, instead of opposing directions, we compare
% activity during autogenous and heterogenous duet song elements. 
%
% Usage: foo = wDSI(wren, padding)
% wren is the structure in ChronicCompleat2019a.mat - See end for a discription.
% padding is shift of syllable windows earlier (negative) or later (positive) in seconds.
% Default is 0. See description below for more information.
% THIS SCRIPT CALCULATES ON A BIRD-BY-BIRD BASIS (DSI values are calcuated from data
% obtained in single birds). 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The boundaries of each syllable is used for this analysis. This is
% inherently problematic as we expect pre-motor activity in awake animals to occur PRIOR to
% the sound and auditory activity in urethane-anesthetized animals to occur AFTER the sound.
% We are comfortable with this approach because it is a rather unbiased. Change
% the value of padding in seconds (e.g. 0.005 or -0.003) to look at the effects.

    padding = 0.000; 

    if nargin == 2; padding = pdnum; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize our holders for dsi values (male chronic and acute, female chronic and acute)

mChron = []; mAcute = []; fChron = []; fAcute = [];

    
%% The pairs of birds

% This function is only for syllables within the duet. It excludes solo
% syllables that occured before the duet. These data are restricted to
% songs that were sung under chronic and then played back successfully
% under acute. 
% These limited data provide three critical advantages.
% First, responses to syllables are only compared to syllables in that
% duet, which is a form of control for differences in firing rates between
% recording sites and variance over time. Second, this approach does not
% use a calculation of spontaneous rates. This is possible only because values
% are based on single birds. The consequence is that the differences might be 
% somewhat diluted. Third, because we are sampling the same HVC both during 
% awake and urethane anesthetized recordings, we can use a 2-sample t-test. 
% Happiness and Joy.

%% List of Chronic and Acute duet singing with syllable indices 
[~, mduetsyls, ~, fduetsyls, ~, ~] = wData;


% For each of the duets we calculate the DSI separately for the female and
% the male. 

% m17: March 2017, male index 1, female index 2

    [mChron(end+1), mAcute(end+1)] = dsi(wren(1), sort([mduetsyls{1}, fduetsyls{1}]), padding);
    [fChron(end+1), fAcute(end+1)] = dsi(wren(2), sort([mduetsyls{1}, fduetsyls{1}]), padding);
%    [fChron(end+1), fAcute(end+1)] = dsi(wren(2), [5 6 7 8 9 10 11 12 13 14], padding);

% j160806: January 2016, 08:06am, male index 3, female index 4

    [mChron(end+1), mAcute(end+1)] = dsi(wren(3), sort([mduetsyls{2}, fduetsyls{2}]), padding);
    [fChron(end+1), fAcute(end+1)] = dsi(wren(4), sort([mduetsyls{2}, fduetsyls{2}]), padding);
%    [fChron(end+1), fAcute(end+1)] = dsi(wren(4), [3 4 5 6 7 8 9 10 11 12 13], padding);

% j160807: January 2016, 08:07am, male index 5, female index 6

    [mChron(end+1), mAcute(end+1)] = dsi(wren(5), sort([mduetsyls{3}, fduetsyls{3}]), padding);
    [fChron(end+1), fAcute(end+1)] = dsi(wren(6), sort([mduetsyls{3}, fduetsyls{3}]), padding);

% j160815: January 2016, 08:15am, male index 7, female index 8

    [mChron(end+1), mAcute(end+1)] = dsi(wren(7), sort([mduetsyls{4}, fduetsyls{4}]), padding);
    [fChron(end+1), fAcute(end+1)] = dsi(wren(8), sort([mduetsyls{4}, fduetsyls{4}]), padding);

% j161009: January 2016, 10:09am, male index 9, female index 10

    [mChron(end+1), mAcute(end+1)] = dsi(wren(9),  sort([mduetsyls{5}, fduetsyls{5}]), padding);
    [fChron(end+1), fAcute(end+1)] = dsi(wren(10), sort([mduetsyls{5}, fduetsyls{5}]), padding);

% j161022: January 2016, 10:22am, male index 11, female index 12

    [mChron(end+1), mAcute(end+1)] = dsi(wren(11), sort([mduetsyls{6}, fduetsyls{6}]), padding);
    [fChron(end+1), fAcute(end+1)] = dsi(wren(12), sort([mduetsyls{6}, fduetsyls{6}]), padding);

% j17060848: 06 January 2017, 08:48am, male index 13, female index 14
    % NOTE THAT THESE DATA ARE MISSING FEMALE ACUTE (failed recording)
    
    [mChron(end+1), mAcute(end+1)] = dsi(wren(13), sort([mduetsyls{7}, fduetsyls{7}]), padding);
    [fChron(end+1), fAcute(end+1)] = dsi(wren(14), sort([mduetsyls{7}, fduetsyls{7}]), padding);
        
% j17081733: 08 January 2017, 5:33pm, male index 15, female index 16
    % Syllables 6 and 11 are both birds at the same time. Because of this,
    % we are being more restrictive and limiting the data to only 2 syllables
    % (one male and one female) of each motif.
    % This is a special case to be highlighted elsewhere.

    [mChron(end+1), mAcute(end+1)] = dsi(wren(15), sort([mduetsyls{8}, fduetsyls{8}]), padding);
    [fChron(end+1), fAcute(end+1)] = dsi(wren(16), sort([mduetsyls{8}, fduetsyls{8}]), padding);
    
    
% Transfer to the output structure. Lazy.
    foo.mChron = mChron; foo.mActue = mAcute; foo.fChron = fChron; foo.fAcute = fAcute;
    
    foo.pad = padding;

    
%% Plotting

% Plot
    figure(1); clf; hold on;
        plot(ones(length(mChron))*1, mChron, 'bo');
        plot(ones(length(fChron))*2, fChron, 'mo');
        plot(ones(length(mAcute))*3, mAcute, 'b*');
        plot(ones(length(fAcute))*4, fAcute, 'm*');

        plot(0.9, mean(mChron), 'k*');
            errorbar(0.9, mean(mChron), std(mChron), 'k' );% /sqrt(length(mChron)));
        plot(1.9, mean(fChron), 'k*');
            errorbar(1.9, mean(fChron), std(fChron), 'k' );% /sqrt(length(fChron)));
        plot(2.9, mean(mAcute), 'k*');
            errorbar(2.9, mean(mAcute), std(mAcute), 'k' );% /sqrt(length(mAcute)));
        plot(3.9, mean(fAcute(1:end-2)), 'k*');
            errorbar(3.9, mean(fAcute(1:end-2)), std(fAcute(1:end-2)), 'k' );% /sqrt(length(fAcute(1:end-1))));
        xlim([0.5 4.5]); ylim([-1 1]);
        plot([0.5 4.5], [0 0], 'k-', 'LineWidth', 0.1);
        plot([2.5 2.5], [-1 1], 'k-', 'LineWidth', 0.1);
        
        text(1.5, -0.8, 'Chronic', 'HorizontalAlignment', 'center');
        text(3.5, -0.8, 'Acute', 'HorizontalAlignment', 'center');

        ylabel('Duet Selectivity Index');
        
%% Statistics

% Statistics - t-test: DSIs different from zero?        

    [foo.mc_diffzero.sig, foo.mc_diffzero.p, foo.mc_diffzero.ci, foo.mc_diffzero.stats] = ttest(mChron);
    [foo.fc_diffzero.sig, foo.fc_diffzero.p, foo.fc_diffzero.ci, foo.fc_diffzero.stats] = ttest(fChron);
    [foo.ma_diffzero.sig, foo.ma_diffzero.p, foo.ma_diffzero.ci, foo.ma_diffzero.stats] = ttest(mAcute);
    [foo.fa_diffzero.sig, foo.fa_diffzero.p, foo.fa_diffzero.ci, foo.fa_diffzero.stats] = ttest(fAcute);

    fprintf('Male Chronic DSI diff from zero? p = %1.2e \n', foo.mc_diffzero.p);
    fprintf('Female Chronic DSI diff from zero? p = %1.2e \n', foo.mc_diffzero.p);
    fprintf('Male Acute DSI diff from zero? p = %1.2e \n', foo.mc_diffzero.p);
    fprintf('Female Acute DSI diff from zero? p = %1.2e \n \n', foo.mc_diffzero.p);

% Statistics - Two sample t-test: do Chronic and Acute differ?

    [foo.mChronVacute.sig, foo.mChronVacute.p, foo.mChronVacute.ci, foo.mChronVacute.stats] = ttest2(mChron, mAcute);
    [foo.fChronVacute.sig, foo.fChronVacute.p, foo.fChronVacute.ci, foo.fChronVacute.stats] = ttest2(fChron, fAcute);

    fprintf('Male Chronic and Acute DSIs differ? p = %1.4f \n', foo.mChronVacute.p);
    fprintf('Female Chronic and Acute DSIs differ? p = %1.4f \n \n', foo.fChronVacute.p);

% Statistics - Two sample t-test: do males differ from females? Test for each Chronic and Acute.

    [foo.cFemVmale.sig, foo.cFemVmale.p, foo.cFemVmale.ci, foo.cFemVmale.stats] = ttest2(mChron, fChron);
    [foo.aFemVmale.sig, foo.aFemVmale.p, foo.aFemVmale.ci, foo.aFemVmale.stats] = ttest2(mAcute, fAcute);

    fprintf('Chronic Female and Chronic Male DSIs differ? p = %1.2e \n', foo.cFemVmale.p);
    fprintf('Acute Female and Acute Male DSIs differ? p = %1.2e \n', foo.aFemVmale.p);
        
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% DSI nested function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
function [chron, acute] = dsi(struct, sylist, pad)
% Calculates dsi which is: 
%   (autogenous - heterogenous) / (autogenous + heterogenous)
% where autogenous and heterogenous are the spike rates during those elements.
% sylist is the idx numbers of the syllables to be used in the analysis.
% pad is the shift (in seconds) of the syllable boundaries (negative are
% earlier in the duet and positive are later in the duet).

% We use spike rates because the numbers of replay repetitions under
% urethane - 10s of playbacks - is different from the 4 electrodes of
% TOE data from awake birds.

% Initialize initialize!

    mcspikes = 0; % Male Chronic spike count (awake)
    maspikes = 0; % Male Acute spike count (urethane)
    fcspikes = 0; % Female Chronic spike count (awake)
    faspikes = 0; % Female Acute spike count (urethane)
    mdur = 0; % Cumulative duration of male syllables
    fdur = 0; % Cumulative duration of female syllables

%% Count spikes for each syllable in the list    
    
    for j=1:length(sylist) % For every duet syllable as specified by the user in rango

            ChronSpkCnt = 0; % Temporary spike count for Awake electrodes in the cell array of TOEs
            AcutSpkcnt = 0; % Temporary spike count for Chronic replay reps in the cell array of TOEs

            % Count the number of spikes in the clicked syllable range 
            for k = 1:length(struct.Cspikes)
                ChronSpkCnt = ChronSpkCnt + length(find(struct.Cspikes{k} > struct.syl(sylist(j)).tim(1)-pad & ...
                    struct.Cspikes{k} < struct.syl(sylist(j)).tim(2))-pad);
            end
            for k = 1:length(struct.Aspikes)
                AcutSpkcnt = AcutSpkcnt + length(find(struct.Aspikes{k} > struct.syl(sylist(j)).tim(1)+pad & ...
                    struct.Aspikes{k} < struct.syl(sylist(j)).tim(2))+pad);
            end

            % Sort into male and female.
            if struct.sylsex(sylist(j)) == 1 % We have a male syllable
                mcspikes = mcspikes + ChronSpkCnt;
                maspikes = maspikes + AcutSpkcnt;
                mdur = mdur + (struct.syl(sylist(j)).tim(2) - struct.syl(sylist(j)).tim(1));
            end

            if struct.sylsex(sylist(j)) == 2 % We have a female syllable
                fcspikes = fcspikes + ChronSpkCnt;
                faspikes = faspikes + AcutSpkcnt;
                fdur = fdur + (struct.syl(sylist(j)).tim(2) - struct.syl(sylist(j)).tim(1));
            end

            if struct.sylsex(sylist(j)) == 3 % Both birds sang at the same time
                %fcspikes = fcspikes + Ctmpspikecount; 
                %mcspikes = mcspikes + ChronSpkCnt;
                %fdur = fdur + (struct(i).syl(j).tim(2) - struct(i).syl(j).tim(1));
                %mdur = mdur + (struct.syl(sylist(j)).tim(2) - struct.syl(sylist(j)).tim(1));
            end
            
    end

%% Convert to spike rate (spikes per second)
    
        mcspikerate = mcspikes / mdur;
        fcspikerate = fcspikes / fdur;
        maspikerate = maspikes / mdur;
        faspikerate = faspikes / fdur;

%% Calculate DSI

        if struct.sexy == 1 % We have a male wren            
            chron = (mcspikerate - fcspikerate) / (mcspikerate + fcspikerate);            
            acute = (maspikerate - faspikerate) / (maspikerate + faspikerate); 
        end
        
        if struct.sexy == 2 % We have a female wren         
            chron = (fcspikerate - mcspikerate) / (mcspikerate + fcspikerate);            
            acute = (faspikerate - maspikerate) / (maspikerate + faspikerate); 
        end
        
end % End dsi nested function



end

% Plot with a different layout
% 
%     figure(2); clf; hold on;
%         plot(ones(length(mChron))*1, mChron, 'bo');
%         plot(ones(length(fChron))*3, fChron, 'mo');
%         plot(ones(length(mAcute))*2, mAcute, 'b*');
%         plot(ones(length(fAcute))*4, fAcute, 'm*');
% 
%         plot(0.9, mean(mChron), 'k*');
%             errorbar(0.9, mean(mChron), std(mChron), 'k' );% /sqrt(length(mChron)));
%         plot(2.9, mean(fChron), 'k*');
%             errorbar(2.9, mean(fChron), std(fChron), 'k' );% /sqrt(length(fChron)));
%         plot(1.9, mean(mAcute), 'k*');
%             errorbar(1.9, mean(mAcute), std(mAcute), 'k' );% /sqrt(length(mAcute)));
%         plot(3.9, mean(fAcute(1:end-1)), 'k*');
%             errorbar(3.9, mean(fAcute(1:end-1)), std(fAcute(1:end-1)), 'k' );% /sqrt(length(fAcute(1:end-1))));
%         
%         xlim([0.5 4.5]); ylim([-1 1]);
%         plot([0.5 4.5], [0 0], 'k-', 'LineWidth', 0.1);
%         plot([2.5 2.5], [-1 1], 'k-', 'LineWidth', 0.1);
%         
%         text(1.5, -0.8, 'Male', 'HorizontalAlignment', 'center');
%         text(3.5, -0.8, 'Female', 'HorizontalAlignment', 'center');
% 
%         ylabel('Duet Selectivity Index');

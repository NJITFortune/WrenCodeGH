%% Duet 9, 0734 birds ***********************************
% close all; clear all;
% load /Users/eric/ChronicCompleat2018c.mat % OLD DATA
% 
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0734/f_0734.mat
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0734/Duet_Female_2016-01-13T07-34-30aClicked.mat % 
%     femout = out; clear out
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0734/m_0734.mat
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0734/Duet_Male_2016-01-13T07-34-32_2aClicked.mat % 
%     malout = out; clear out
% 
%     w(18).wrensex = 'F'; 
%     w(17).wrensex = 'M'; 
%     w(18).sexy = 2; 
%     w(17).sexy = 1;
%     w(18).id = 'j160734'; 
%     w(17).id = 'm160734';
%     w(18).Aspikes = {}; %  NO ACUTE DATA
%     w(17).Aspikes = {}; %  NO ACUTE DATA
%     w(18).Fs = 1/Duet_Female_2016_01_13T07_34_30_32_bit__Ch5.interval; 
%     w(17).Fs = 1/Duet_Female_2016_01_13T07_34_30_32_bit__Ch5.interval;
%     w(18).duet = Duet_Female_2016_01_13T07_34_30_32_bit__Ch5.values;
%     w(17).duet = w(18).duet;
%     w(18).tim = 1/w(18).Fs:1/w(18).Fs:Duet_Female_2016_01_13T07_34_30_32_bit__Ch5.length * Duet_Female_2016_01_13T07_34_30_32_bit__Ch5.interval;
%     w(17).tim = w(18).tim;
%     
%     % Female chronic spikes
%     w(18).Cspikes{1} = Duet_Female_2016_01_13T07_34_30_32_bit__Ch13.times;
%     w(18).Cspikes{2} = Duet_Female_2016_01_13T07_34_30_32_bit__Ch15.times;
%     w(18).Cspikes{3} = Duet_Female_2016_01_13T07_34_30_32_bit__Ch17.times;
%     w(18).Cspikes{4} = Duet_Female_2016_01_13T07_34_30_32_bit__Ch18.times;
%     % Male chronic spikes
%     w(17).Cspikes{1} = Duet_Male_2016_01_13T07_34_32_2_32_bit__Ch15.times;
%     w(17).Cspikes{2} = Duet_Male_2016_01_13T07_34_32_2_32_bit__Ch16.times;
%     w(17).Cspikes{3} = Duet_Male_2016_01_13T07_34_32_2_32_bit__Ch17.times;
%     w(17).Cspikes{4} = Duet_Male_2016_01_13T07_34_32_2_32_bit__Ch18.times;
% 
% 
%     % Fix the syllables to 'modern' standards
%     w(18).syl(1).tim(1) = femout(1).syltim(1);
%     w(18).syl(1).tim(2) = femout(3).syltim(2);
%     w(18).syl(2).tim(1) = femout(4).syltim(1);
%     w(18).syl(2).tim(2) = femout(6).syltim(2);
%     w(18).syl(3).tim(1) = femout(7).syltim(1);
%     w(18).syl(3).tim(2) = femout(9).syltim(2);
%     w(18).syl(4).tim(1) = femout(10).syltim(1);
%     w(18).syl(4).tim(2) = femout(10).syltim(2);    
% 
%     w(17).syl = w(18).syl;
% 
%     w(18).sylsex = [2 2 2 2];
%     w(17).sylsex = w(18).sylsex;
% 
% % Confirm everything plot
%     figure(2); clf;
%         xa(1) = subplot(411); hold on;
%         bs_raster(w(18).Cspikes, 'm');
%         xa(2) = subplot(412); hold on;
%         bs_raster(w(17).Cspikes, 'b');
%         xa(3) = subplot(413); hold on;        
%         plot(w(18).tim, 2.5+(2.5*w(18).duet/(max(abs(w(17).duet)) )), 'k');
%         xa(4) = subplot(414); hold on;        
%         specgram(w(17).duet, 512, w(18).Fs, [], 500); ylim([500 4500]); colormap('HOT'); caxis([-10 40]);        
%         
%         linkaxes(xa, 'x');
%         
%  for i = 1:length(w(18).syl)
%     for j = 1:3
%         subplot(4,1,j);
%         plot([w(17).syl(i).tim(1) w(18).syl(i).tim(1)], [0 5], 'g');
%         plot([w(18).syl(i).tim(2) w(17).syl(i).tim(2)], [0 5], 'r');
%     end
%     subplot(414);         
%         plot([w(17).syl(i).tim(1) w(18).syl(i).tim(1)], [500 4500], 'g');
%         plot([w(18).syl(i).tim(2) w(17).syl(i).tim(2)], [500 4500], 'r');
% 
%  end
%         
% save /Volumes/ExFAT/Chronically/ChronicCompleat2018d.mat w; 



%% Duet 10, 0800a birds ***********************************
% close all; clear all;
% load /Users/eric/ChronicCompleat2018d.mat % OLD DATA
% 
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0800/Duet_Female_2016-01-13T08-00-14a.mat
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0800/Duet_Female_2016-01-13T08-00-14aClicked.mat % 
%     femout = out; clear out
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0800/Duet_Male_2016-01-13T08-00-16_2a.mat
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0800/Duet_Male_2016-01-13T08-00-16_2aClicked.mat % 
%     malout = out; clear out
% 
%     w(20).wrensex = 'F'; 
%     w(19).wrensex = 'M'; 
%     w(20).sexy = 2; 
%     w(19).sexy = 1;
%     w(20).id = 'j160800a'; 
%     w(19).id = 'm160800a';
%     w(20).Aspikes = {}; %  NO ACUTE DATA
%     w(19).Aspikes = {}; %  NO ACUTE DATA
%     w(20).Fs = 1/Duet_Female_2016_01_13T08_00_14_Ch5.interval; 
%     w(19).Fs = 1/Duet_Female_2016_01_13T08_00_14_Ch5.interval;
%     w(20).duet = Duet_Female_2016_01_13T08_00_14_Ch5.values;
%     w(19).duet = w(20).duet;
%     w(20).tim = 1/w(20).Fs:1/w(20).Fs:Duet_Female_2016_01_13T08_00_14_Ch5.length * Duet_Female_2016_01_13T08_00_14_Ch5.interval;
%     w(19).tim = w(20).tim;
%     
%     % Female chronic spikes
%     w(20).Cspikes{1} = Duet_Female_2016_01_13T08_00_14_Ch11.times;
%     w(20).Cspikes{2} = Duet_Female_2016_01_13T08_00_14_Ch12.times;
%     w(20).Cspikes{3} = Duet_Female_2016_01_13T08_00_14_Ch13.times;
%     w(20).Cspikes{4} = Duet_Female_2016_01_13T08_00_14_Ch14.times;
%     % Male chronic spikes
%     w(19).Cspikes{1} = Duet_Male_2016_01_13T08_00_16_2_Ch11.times;
%     w(19).Cspikes{2} = Duet_Male_2016_01_13T08_00_16_2_Ch12.times;
%     w(19).Cspikes{3} = Duet_Male_2016_01_13T08_00_16_2_Ch13.times;
%     w(19).Cspikes{4} = Duet_Male_2016_01_13T08_00_16_2_Ch14.times;
% 
% 
%     % Fix the syllables to 'modern' standards
%     w(20).syl(1).tim(1) = femout(1).syltim(1);
%     w(20).syl(1).tim(2) = femout(3).syltim(2);
% 
%     w(19).syl = w(20).syl;
% 
%     w(20).sylsex = 2;
%     w(19).sylsex = w(20).sylsex;
% 
% % Confirm everything plot
%     figure(2); clf;
%         xa(1) = subplot(411); hold on;
%         bs_raster(w(20).Cspikes, 'm');
%         xa(2) = subplot(412); hold on;
%         bs_raster(w(19).Cspikes, 'b');
%         xa(3) = subplot(413); hold on;        
%         plot(w(20).tim, 2.5+(2.5*w(20).duet/(max(abs(w(19).duet)) )), 'k');
%         xa(4) = subplot(414); hold on;        
%         specgram(w(19).duet, 512, w(20).Fs, [], 500); ylim([500 4500]); colormap('HOT'); caxis([-10 40]);        
%         
%         linkaxes(xa, 'x');
%         
%  for i = 1:length(w(20).syl)
%     for j = 1:3
%         subplot(4,1,j);
%         plot([w(19).syl(i).tim(1) w(20).syl(i).tim(1)], [0 5], 'g');
%         plot([w(20).syl(i).tim(2) w(19).syl(i).tim(2)], [0 5], 'r');
%     end
%     subplot(414);         
%         plot([w(19).syl(i).tim(1) w(19).syl(i).tim(1)], [500 4500], 'g');
%         plot([w(20).syl(i).tim(2) w(20).syl(i).tim(2)], [500 4500], 'r');
% 
%  end
        
% save /Volumes/ExFAT/Chronically/ChronicCompleat2018d.mat w; 

%% Duet 11, 0800b birds ***********************************
% close all; clear all;
% load /Users/eric/ChronicCompleat2018d.mat % OLD DATA
% 
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0800/Duet_Female_2016-01-13T08-00-14b.mat
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0800/Duet_Female_2016-01-13T08-00-14bClicked.mat % 
%     femout = out; clear out
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0800/Duet_Male_2016-01-13T08-00-16_2b.mat
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0800/Duet_Male_2016-01-13T08-00-16_2bClicked.mat % 
%     malout = out; clear out
% 
%     w(22).wrensex = 'F'; 
%     w(21).wrensex = 'M'; 
%     w(22).sexy = 2; 
%     w(21).sexy = 1;
%     w(22).id = 'j160800b'; 
%     w(21).id = 'm160800b';
%     w(22).Aspikes = {}; %  NO ACUTE DATA
%     w(21).Aspikes = {}; %  NO ACUTE DATA
%     w(22).Fs = 1/Duet_Female_2016_01_13T08_00_14_Ch5.interval; 
%     w(21).Fs = 1/Duet_Female_2016_01_13T08_00_14_Ch5.interval;
%     w(22).duet = Duet_Female_2016_01_13T08_00_14_Ch5.values;
%     w(21).duet = w(22).duet;
%     w(22).tim = 1/w(22).Fs:1/w(22).Fs:Duet_Female_2016_01_13T08_00_14_Ch5.length * Duet_Female_2016_01_13T08_00_14_Ch5.interval;
%     w(21).tim = w(22).tim;
%     
%     % Female chronic spikes
%     w(22).Cspikes{1} = Duet_Female_2016_01_13T08_00_14_Ch11.times;
%     w(22).Cspikes{2} = Duet_Female_2016_01_13T08_00_14_Ch12.times;
%     w(22).Cspikes{3} = Duet_Female_2016_01_13T08_00_14_Ch13.times;
%     w(22).Cspikes{4} = Duet_Female_2016_01_13T08_00_14_Ch14.times;
%     % Male chronic spikes
%     w(21).Cspikes{1} = Duet_Male_2016_01_13T08_00_16_2_Ch11.times;
%     w(21).Cspikes{2} = Duet_Male_2016_01_13T08_00_16_2_Ch12.times;
%     w(21).Cspikes{3} = Duet_Male_2016_01_13T08_00_16_2_Ch13.times;
%     w(21).Cspikes{4} = Duet_Male_2016_01_13T08_00_16_2_Ch14.times;
% 
% 
%     % Fix the syllables to 'modern' standards
%     w(22).syl(1).tim(1) = femout(1).syltim(1);
%     w(22).syl(1).tim(2) = femout(3).syltim(2);
% 
%     w(21).syl = w(22).syl;
% 
%     w(22).sylsex = 2;
%     w(21).sylsex = w(22).sylsex;
% 
% % Confirm everything plot
%     figure(2); clf;
%         xa(1) = subplot(411); hold on;
%         bs_raster(w(22).Cspikes, 'm');
%         xa(2) = subplot(412); hold on;
%         bs_raster(w(21).Cspikes, 'b');
%         xa(3) = subplot(413); hold on;        
%         plot(w(22).tim, 2.5+(2.5*w(22).duet/(max(abs(w(21).duet)) )), 'k');
%         xa(4) = subplot(414); hold on;        
%         specgram(w(21).duet, 512, w(22).Fs, [], 500); ylim([500 4500]); colormap('HOT'); caxis([-10 40]);        
%         
%         linkaxes(xa, 'x');
%         
%  for i = 1:length(w(22).syl)
%     for j = 1:3
%         subplot(4,1,j);
%         plot([w(21).syl(i).tim(1) w(22).syl(i).tim(1)], [0 5], 'g');
%         plot([w(22).syl(i).tim(2) w(21).syl(i).tim(2)], [0 5], 'r');
%     end
%     subplot(414);         
%         plot([w(21).syl(i).tim(1) w(21).syl(i).tim(1)], [500 4500], 'g');
%         plot([w(22).syl(i).tim(2) w(22).syl(i).tim(2)], [500 4500], 'r');
% 
%  end
%         
% save /Volumes/ExFAT/Chronically/ChronicCompleat2018d.mat w; 


%% Duet 12, 0928a birds ***********************************
% close all; clear all;
% load /Users/eric/ChronicCompleat2018d.mat % OLD DATA
% 
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0928/Duet_Female_2016-01-13T09-28-38a.mat
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0928/Duet_Female_2016-01-13T09-28-38aClicked.mat % 
%     femout = out; clear out
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0928/Duet_Male_2016-01-13T09-28-41_2a.mat
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0928/Duet_Male_2016-01-13T09-28-41_2aClicked.mat % 
%     malout = out; clear out
% 
%     f = 24; m = 23;
%     
%     w(f).wrensex = 'F'; 
%     w(m).wrensex = 'M'; 
%     w(f).sexy = 2; 
%     w(m).sexy = 1;
%     w(f).id = 'j160928a'; 
%     w(m).id = 'm160928a';
%     w(f).Aspikes = {}; %  NO ACUTE DATA
%     w(m).Aspikes = {}; %  NO ACUTE DATA
%     w(f).Fs = 1/Duet_Female_2016_01_13T09_28_38_Ch5.interval; 
%     w(m).Fs = 1/Duet_Female_2016_01_13T09_28_38_Ch5.interval;
%     w(f).duet = Duet_Female_2016_01_13T09_28_38_Ch5.values;
%     w(m).duet = w(f).duet;
%     w(f).tim = 1/w(f).Fs:1/w(f).Fs:Duet_Female_2016_01_13T09_28_38_Ch5.length * Duet_Female_2016_01_13T09_28_38_Ch5.interval;
%     w(m).tim = w(f).tim;
%     
%     % Female chronic spikes
%     w(f).Cspikes{1} = Duet_Female_2016_01_13T09_28_38_Ch11.times;
%     w(f).Cspikes{2} = Duet_Female_2016_01_13T09_28_38_Ch12.times;
%     w(f).Cspikes{3} = Duet_Female_2016_01_13T09_28_38_Ch14.times;
%     w(f).Cspikes{4} = Duet_Female_2016_01_13T09_28_38_Ch15.times;
%     % Male chronic spike
%     w(m).Cspikes{1} = Duet_Male_2016_01_13T09_28_41_2_Ch11.times;
%     w(m).Cspikes{2} = Duet_Male_2016_01_13T09_28_41_2_Ch12.times;
%     w(m).Cspikes{3} = Duet_Male_2016_01_13T09_28_41_2_Ch13.times;
%     w(m).Cspikes{4} = Duet_Male_2016_01_13T09_28_41_2_Ch13.times;
% 
% 
%     % Fix the syllables to 'modern' standards
%     w(f).syl(1).tim(1) = femout(1).syltim(1);
%     w(f).syl(1).tim(2) = femout(3).syltim(2);
% 
%     w(m).syl = w(f).syl;
% 
%     w(f).sylsex = 2;
%     w(m).sylsex = w(f).sylsex;
% 
% % Confirm everything plot
%     figure(2); clf;
%         xa(1) = subplot(411); hold on;
%         bs_raster(w(f).Cspikes, 'm');
%         xa(2) = subplot(412); hold on;
%         bs_raster(w(m).Cspikes, 'b');
%         xa(3) = subplot(413); hold on;        
%         plot(w(f).tim, 2.5+(2.5*w(f).duet/(max(abs(w(m).duet)) )), 'k');
%         xa(4) = subplot(414); hold on;        
%         specgram(w(m).duet, 512, w(f).Fs, [], 500); ylim([500 4500]); colormap('HOT'); caxis([-10 40]);        
%         
%         linkaxes(xa, 'x');
%         
%  for i = 1:length(w(f).syl)
%     for j = 1:3
%         subplot(4,1,j);
%         plot([w(f).syl(i).tim(1) w(f).syl(i).tim(1)], [0 5], 'g');
%         plot([w(m).syl(i).tim(2) w(m).syl(i).tim(2)], [0 5], 'r');
%     end
%     subplot(414);         
%         plot([w(m).syl(i).tim(1) w(f).syl(i).tim(1)], [500 4500], 'g');
%         plot([w(f).syl(i).tim(2) w(m).syl(i).tim(2)], [500 4500], 'r');
% 
%  end
%         
% save /Volumes/ExFAT/Chronically/ChronicCompleat2018d.mat w; 


%% Duet 13, 0928b birds ***********************************
% close all; clear all;
% load /Users/eric/ChronicCompleat2018d.mat % OLD DATA
% 
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0928/Duet_Female_2016-01-13T09-28-38b.mat
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0928/Duet_Female_2016-01-13T09-28-38bClicked.mat % 
%     femout = out; clear out
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0928/Duet_Male_2016-01-13T09-28-41_2b.mat
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0928/Duet_Male_2016-01-13T09-28-41_2bClicked.mat % 
%     malout = out; clear out
% 
%     f = 26; m = 25;
%     
%     w(f).wrensex = 'F'; 
%     w(m).wrensex = 'M'; 
%     w(f).sexy = 2; 
%     w(m).sexy = 1;
%     w(f).id = 'j160928b'; 
%     w(m).id = 'm160928b';
%     w(f).Aspikes = {}; %  NO ACUTE DATA
%     w(m).Aspikes = {}; %  NO ACUTE DATA
%     w(f).Fs = 1/Duet_Female_2016_01_13T09_28_38_Ch5.interval; 
%     w(m).Fs = 1/Duet_Female_2016_01_13T09_28_38_Ch5.interval;
%     w(f).duet = Duet_Female_2016_01_13T09_28_38_Ch5.values;
%     w(m).duet = w(f).duet;
%     w(f).tim = 1/w(f).Fs:1/w(f).Fs:Duet_Female_2016_01_13T09_28_38_Ch5.length * Duet_Female_2016_01_13T09_28_38_Ch5.interval;
%     w(m).tim = w(f).tim;
%     
%     % Female chronic spikes
%     w(f).Cspikes{1} = Duet_Female_2016_01_13T09_28_38_Ch11.times;
%     w(f).Cspikes{2} = Duet_Female_2016_01_13T09_28_38_Ch12.times;
%     w(f).Cspikes{3} = Duet_Female_2016_01_13T09_28_38_Ch14.times;
%     w(f).Cspikes{4} = Duet_Female_2016_01_13T09_28_38_Ch15.times;
%     % Male chronic spike
%     w(m).Cspikes{1} = Duet_Male_2016_01_13T09_28_41_2_Ch11.times;
%     w(m).Cspikes{2} = Duet_Male_2016_01_13T09_28_41_2_Ch12.times;
%     w(m).Cspikes{3} = Duet_Male_2016_01_13T09_28_41_2_Ch13.times;
%     w(m).Cspikes{4} = Duet_Male_2016_01_13T09_28_41_2_Ch13.times;
% 
% 
%     % Fix the syllables to 'modern' standards
%     w(f).syl(1).tim(1) = femout(1).syltim(1);
%     w(f).syl(1).tim(2) = femout(4).syltim(2);
% 
%     w(m).syl = w(f).syl;
% 
%     w(f).sylsex = 2;
%     w(m).sylsex = w(f).sylsex;
% 
% % Confirm everything plot
%     figure(2); clf;
%         xa(1) = subplot(411); hold on;
%         bs_raster(w(f).Cspikes, 'm');
%         xa(2) = subplot(412); hold on;
%         bs_raster(w(m).Cspikes, 'b');
%         xa(3) = subplot(413); hold on;        
%         plot(w(f).tim, 2.5+(2.5*w(f).duet/(max(abs(w(m).duet)) )), 'k');
%         xa(4) = subplot(414); hold on;        
%         specgram(w(m).duet, 512, w(f).Fs, [], 500); ylim([500 4500]); colormap('HOT'); caxis([-10 40]);        
%         
%         linkaxes(xa, 'x');
%         
%  for i = 1:length(w(f).syl)
%     for j = 1:3
%         subplot(4,1,j);
%         plot([w(f).syl(i).tim(1) w(f).syl(i).tim(1)], [0 5], 'g');
%         plot([w(m).syl(i).tim(2) w(m).syl(i).tim(2)], [0 5], 'r');
%     end
%     subplot(414);         
%         plot([w(m).syl(i).tim(1) w(f).syl(i).tim(1)], [500 4500], 'g');
%         plot([w(f).syl(i).tim(2) w(m).syl(i).tim(2)], [500 4500], 'r');
% 
%  end
%         
% save /Volumes/ExFAT/Chronically/ChronicCompleat2018d.mat w; 


%% Duet 14, 0934a birds ***********************************
% close all; clear all;
% load /Users/eric/ChronicCompleat2018d.mat % OLD DATA
% 
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0934/Duet_Female_2016-01-13T09-34-01a.mat
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0934/Duet_Female_2016-01-13T09-34-01aClicked.mat % 
%     femout = out; clear out
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0934/Duet_Male_2016-01-13T09-34-03_2a.mat
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0934/Duet_Male_2016-01-13T09-34-03_2aClicked.mat % 
%     malout = out; clear out
% 
%     f = 28; m = 27;
%     
%     w(f).wrensex = 'F'; 
%     w(m).wrensex = 'M'; 
%     w(f).sexy = 2; 
%     w(m).sexy = 1;
%     w(f).id = 'j160934a'; 
%     w(m).id = 'm160934a';
%     w(f).Aspikes = {}; %  NO ACUTE DATA
%     w(m).Aspikes = {}; %  NO ACUTE DATA
%     w(f).Fs = 1/Duet_Female_2016_01_13T09_34_01_Ch5.interval; 
%     w(m).Fs = 1/Duet_Female_2016_01_13T09_34_01_Ch5.interval;
%     w(f).duet = Duet_Female_2016_01_13T09_34_01_Ch5.values;
%     w(m).duet = w(f).duet;
%     w(f).tim = 1/w(f).Fs:1/w(f).Fs:Duet_Female_2016_01_13T09_34_01_Ch5.length * Duet_Female_2016_01_13T09_34_01_Ch5.interval;
%     w(m).tim = w(f).tim;
%     
%     % Female chronic spikes
%     w(f).Cspikes{1} = Duet_Female_2016_01_13T09_34_01_Ch12.times;
%     w(f).Cspikes{2} = Duet_Female_2016_01_13T09_34_01_Ch13.times;
%     w(f).Cspikes{3} = Duet_Female_2016_01_13T09_34_01_Ch14.times;
%     w(f).Cspikes{4} = Duet_Female_2016_01_13T09_34_01_Ch15.times;
%     % Male chronic spike
%     w(m).Cspikes{1} = Duet_Male_2016_01_13T09_34_03_2_Ch11.times;
%     w(m).Cspikes{2} = Duet_Male_2016_01_13T09_34_03_2_Ch12.times;
%     w(m).Cspikes{3} = Duet_Male_2016_01_13T09_34_03_2_Ch14.times;
%     w(m).Cspikes{4} = Duet_Male_2016_01_13T09_34_03_2_Ch15.times;
% 
% 
%     % Fix the syllables to 'modern' standards
%     w(f).syl(1).tim(1) = femout(1).syltim(1);
%     w(f).syl(1).tim(2) = femout(3).syltim(2);
% 
%     w(m).syl = w(f).syl;
% 
%     w(f).sylsex = 2;
%     w(m).sylsex = w(f).sylsex;
% 
% % Confirm everything plot
%     figure(2); clf;
%         xa(1) = subplot(411); hold on;
%         bs_raster(w(f).Cspikes, 'm');
%         xa(2) = subplot(412); hold on;
%         bs_raster(w(m).Cspikes, 'b');
%         xa(3) = subplot(413); hold on;        
%         plot(w(f).tim, 2.5+(2.5*w(f).duet/(max(abs(w(m).duet)) )), 'k');
%         xa(4) = subplot(414); hold on;        
%         specgram(w(m).duet, 512, w(f).Fs, [], 500); ylim([500 4500]); colormap('HOT'); caxis([-10 40]);        
%         
%         linkaxes(xa, 'x');
%         
%  for i = 1:length(w(f).syl)
%     for j = 1:3
%         subplot(4,1,j);
%         plot([w(f).syl(i).tim(1) w(f).syl(i).tim(1)], [0 5], 'g');
%         plot([w(m).syl(i).tim(2) w(m).syl(i).tim(2)], [0 5], 'r');
%     end
%     subplot(414);         
%         plot([w(m).syl(i).tim(1) w(f).syl(i).tim(1)], [500 4500], 'g');
%         plot([w(f).syl(i).tim(2) w(m).syl(i).tim(2)], [500 4500], 'r');
% 
%  end
%         
% save /Volumes/ExFAT/Chronically/ChronicCompleat2018d.mat w; 


%% Duet 14, 0934b birds ***********************************
% close all; clear all;
% load /Users/eric/ChronicCompleat2018d.mat % OLD DATA
% 
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0934/Duet_Female_2016-01-13T09-34-01b.mat
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0934/Duet_Female_2016-01-13T09-34-01bClicked.mat % 
%     femout = out; clear out
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0934/Duet_Male_2016-01-13T09-34-03_2b.mat
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0934/Duet_Male_2016-01-13T09-34-03_2bClicked.mat % 
%     malout = out; clear out
% 
%     f = 30; m = 29;
%     
%     w(f).wrensex = 'F'; 
%     w(m).wrensex = 'M'; 
%     w(f).sexy = 2; 
%     w(m).sexy = 1;
%     w(f).id = 'j160934b'; 
%     w(m).id = 'm160934b';
%     w(f).Aspikes = {}; %  NO ACUTE DATA
%     w(m).Aspikes = {}; %  NO ACUTE DATA
%     w(f).Fs = 1/Duet_Female_2016_01_13T09_34_01_Ch5.interval; 
%     w(m).Fs = 1/Duet_Female_2016_01_13T09_34_01_Ch5.interval;
%     w(f).duet = Duet_Female_2016_01_13T09_34_01_Ch5.values;
%     w(m).duet = w(f).duet;
%     w(f).tim = 1/w(f).Fs:1/w(f).Fs:Duet_Female_2016_01_13T09_34_01_Ch5.length * Duet_Female_2016_01_13T09_34_01_Ch5.interval;
%     w(m).tim = w(f).tim;
%     
%     % Female chronic spikes
%     w(f).Cspikes{1} = Duet_Female_2016_01_13T09_34_01_Ch12.times;
%     w(f).Cspikes{2} = Duet_Female_2016_01_13T09_34_01_Ch13.times;
%     w(f).Cspikes{3} = Duet_Female_2016_01_13T09_34_01_Ch14.times;
%     w(f).Cspikes{4} = Duet_Female_2016_01_13T09_34_01_Ch15.times;
%     % Male chronic spike
%     w(m).Cspikes{1} = Duet_Male_2016_01_13T09_34_03_2_Ch11.times;
%     w(m).Cspikes{2} = Duet_Male_2016_01_13T09_34_03_2_Ch12.times;
%     w(m).Cspikes{3} = Duet_Male_2016_01_13T09_34_03_2_Ch14.times;
%     w(m).Cspikes{4} = Duet_Male_2016_01_13T09_34_03_2_Ch15.times;
% 
% 
%     % Fix the syllables to 'modern' standards
%     w(f).syl(1).tim(1) = femout(1).syltim(1);
%     w(f).syl(1).tim(2) = femout(3).syltim(2);
% 
%     w(m).syl = w(f).syl;
% 
%     w(f).sylsex = 2;
%     w(m).sylsex = w(f).sylsex;
% 
% % Confirm everything plot
%     figure(2); clf;
%         xa(1) = subplot(411); hold on;
%         bs_raster(w(f).Cspikes, 'm');
%         xa(2) = subplot(412); hold on;
%         bs_raster(w(m).Cspikes, 'b');
%         xa(3) = subplot(413); hold on;        
%         plot(w(f).tim, 2.5+(2.5*w(f).duet/(max(abs(w(m).duet)) )), 'k');
%         xa(4) = subplot(414); hold on;        
%         specgram(w(m).duet, 512, w(f).Fs, [], 500); ylim([500 4500]); colormap('HOT'); caxis([-10 40]);        
%         
%         linkaxes(xa, 'x');
%         
%  for i = 1:length(w(f).syl)
%     for j = 1:3
%         subplot(4,1,j);
%         plot([w(f).syl(i).tim(1) w(f).syl(i).tim(1)], [0 5], 'g');
%         plot([w(m).syl(i).tim(2) w(m).syl(i).tim(2)], [0 5], 'r');
%     end
%     subplot(414);         
%         plot([w(m).syl(i).tim(1) w(f).syl(i).tim(1)], [500 4500], 'g');
%         plot([w(f).syl(i).tim(2) w(m).syl(i).tim(2)], [500 4500], 'r');
% 
%  end
%         
% save /Volumes/ExFAT/Chronically/ChronicCompleat2018d.mat w; 

%% Duet 14, 0934c birds ***********************************
% close all; clear all;
% load /Users/eric/ChronicCompleat2018d.mat % OLD DATA
% 
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0934/Duet_Female_2016-01-13T09-34-01c.mat
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0934/Duet_Female_2016-01-13T09-34-01cClicked.mat % 
%     femout = out; clear out
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0934/Duet_Male_2016-01-13T09-34-03_2c.mat
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0934/Duet_Male_2016-01-13T09-34-03_2cClicked.mat % 
%     malout = out; clear out
% 
%     f = 32; m = 31;
%     
%     w(f).wrensex = 'F'; 
%     w(m).wrensex = 'M'; 
%     w(f).sexy = 2; 
%     w(m).sexy = 1;
%     w(f).id = 'j160934c'; 
%     w(m).id = 'm160934c';
%     w(f).Aspikes = {}; %  NO ACUTE DATA
%     w(m).Aspikes = {}; %  NO ACUTE DATA
%     w(f).Fs = 1/Duet_Female_2016_01_13T09_34_01_Ch5.interval; 
%     w(m).Fs = 1/Duet_Female_2016_01_13T09_34_01_Ch5.interval;
%     w(f).duet = Duet_Female_2016_01_13T09_34_01_Ch5.values;
%     w(m).duet = w(f).duet;
%     w(f).tim = 1/w(f).Fs:1/w(f).Fs:Duet_Female_2016_01_13T09_34_01_Ch5.length * Duet_Female_2016_01_13T09_34_01_Ch5.interval;
%     w(m).tim = w(f).tim;
%     
%     % Female chronic spikes
%     w(f).Cspikes{1} = Duet_Female_2016_01_13T09_34_01_Ch12.times;
%     w(f).Cspikes{2} = Duet_Female_2016_01_13T09_34_01_Ch13.times;
%     w(f).Cspikes{3} = Duet_Female_2016_01_13T09_34_01_Ch14.times;
%     w(f).Cspikes{4} = Duet_Female_2016_01_13T09_34_01_Ch15.times;
%     % Male chronic spike
%     w(m).Cspikes{1} = Duet_Male_2016_01_13T09_34_03_2_Ch11.times;
%     w(m).Cspikes{2} = Duet_Male_2016_01_13T09_34_03_2_Ch12.times;
%     w(m).Cspikes{3} = Duet_Male_2016_01_13T09_34_03_2_Ch14.times;
%     w(m).Cspikes{4} = Duet_Male_2016_01_13T09_34_03_2_Ch15.times;
% 
% 
%     % Fix the syllables to 'modern' standards
%     w(f).syl(1).tim(1) = femout(1).syltim(1);
%     w(f).syl(1).tim(2) = femout(3).syltim(2);
% 
%     w(m).syl = w(f).syl;
% 
%     w(f).sylsex = 2;
%     w(m).sylsex = w(f).sylsex;
% 
% % Confirm everything plot
%     figure(2); clf;
%         xa(1) = subplot(411); hold on;
%         bs_raster(w(f).Cspikes, 'm');
%         xa(2) = subplot(412); hold on;
%         bs_raster(w(m).Cspikes, 'b');
%         xa(3) = subplot(413); hold on;        
%         plot(w(f).tim, 2.5+(2.5*w(f).duet/(max(abs(w(m).duet)) )), 'k');
%         xa(4) = subplot(414); hold on;        
%         specgram(w(m).duet, 512, w(f).Fs, [], 500); ylim([500 4500]); colormap('HOT'); caxis([-10 40]);        
%         
%         linkaxes(xa, 'x');
%         
%  for i = 1:length(w(f).syl)
%     for j = 1:3
%         subplot(4,1,j);
%         plot([w(f).syl(i).tim(1) w(f).syl(i).tim(1)], [0 5], 'g');
%         plot([w(m).syl(i).tim(2) w(m).syl(i).tim(2)], [0 5], 'r');
%     end
%     subplot(414);         
%         plot([w(m).syl(i).tim(1) w(f).syl(i).tim(1)], [500 4500], 'g');
%         plot([w(f).syl(i).tim(2) w(m).syl(i).tim(2)], [500 4500], 'r');
% 
%  end
%         
% save /Volumes/ExFAT/Chronically/ChronicCompleat2018d.mat w; 

%% Duet 14, 0934d birds ***********************************
% close all; clear all;
% load /Users/eric/ChronicCompleat2018d.mat % OLD DATA
% 
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0934/Duet_Female_2016-01-13T09-34-01d.mat
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0934/Duet_Female_2016-01-13T09-34-01dClicked.mat % 
%     femout = out; clear out
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0934/Duet_Male_2016-01-13T09-34-03_2d.mat
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0934/Duet_Male_2016-01-13T09-34-03_2dClicked.mat % 
%     malout = out; clear out
% 
%     f = 34; m = 33;
%     
%     w(f).wrensex = 'F'; 
%     w(m).wrensex = 'M'; 
%     w(f).sexy = 2; 
%     w(m).sexy = 1;
%     w(f).id = 'j160934d'; 
%     w(m).id = 'm160934d';
%     w(f).Aspikes = {}; %  NO ACUTE DATA
%     w(m).Aspikes = {}; %  NO ACUTE DATA
%     w(f).Fs = 1/Duet_Female_2016_01_13T09_34_01_Ch5.interval; 
%     w(m).Fs = 1/Duet_Female_2016_01_13T09_34_01_Ch5.interval;
%     w(f).duet = Duet_Female_2016_01_13T09_34_01_Ch5.values;
%     w(m).duet = w(f).duet;
%     w(f).tim = 1/w(f).Fs:1/w(f).Fs:Duet_Female_2016_01_13T09_34_01_Ch5.length * Duet_Female_2016_01_13T09_34_01_Ch5.interval;
%     w(m).tim = w(f).tim;
%     
%     % Female chronic spikes
%     w(f).Cspikes{1} = Duet_Female_2016_01_13T09_34_01_Ch12.times;
%     w(f).Cspikes{2} = Duet_Female_2016_01_13T09_34_01_Ch13.times;
%     w(f).Cspikes{3} = Duet_Female_2016_01_13T09_34_01_Ch14.times;
%     w(f).Cspikes{4} = Duet_Female_2016_01_13T09_34_01_Ch15.times;
%     % Male chronic spike
%     w(m).Cspikes{1} = Duet_Male_2016_01_13T09_34_03_2_Ch11.times;
%     w(m).Cspikes{2} = Duet_Male_2016_01_13T09_34_03_2_Ch12.times;
%     w(m).Cspikes{3} = Duet_Male_2016_01_13T09_34_03_2_Ch14.times;
%     w(m).Cspikes{4} = Duet_Male_2016_01_13T09_34_03_2_Ch15.times;
% 
% 
%     % Fix the syllables to 'modern' standards
%     w(f).syl(1).tim(1) = femout(1).syltim(1);
%     w(f).syl(1).tim(2) = femout(3).syltim(2);
% 
%     w(m).syl = w(f).syl;
% 
%     w(f).sylsex = 2;
%     w(m).sylsex = w(f).sylsex;
% 
% % Confirm everything plot
%     figure(2); clf;
%         xa(1) = subplot(411); hold on;
%         bs_raster(w(f).Cspikes, 'm');
%         xa(2) = subplot(412); hold on;
%         bs_raster(w(m).Cspikes, 'b');
%         xa(3) = subplot(413); hold on;        
%         plot(w(f).tim, 2.5+(2.5*w(f).duet/(max(abs(w(m).duet)) )), 'k');
%         xa(4) = subplot(414); hold on;        
%         specgram(w(m).duet, 512, w(f).Fs, [], 500); ylim([500 4500]); colormap('HOT'); caxis([-10 40]);        
%         
%         linkaxes(xa, 'x');
%         
%  for i = 1:length(w(f).syl)
%     for j = 1:3
%         subplot(4,1,j);
%         plot([w(f).syl(i).tim(1) w(f).syl(i).tim(1)], [0 5], 'g');
%         plot([w(m).syl(i).tim(2) w(m).syl(i).tim(2)], [0 5], 'r');
%     end
%     subplot(414);         
%         plot([w(m).syl(i).tim(1) w(f).syl(i).tim(1)], [500 4500], 'g');
%         plot([w(f).syl(i).tim(2) w(m).syl(i).tim(2)], [500 4500], 'r');
% 
%  end
%         
% save /Volumes/ExFAT/Chronically/ChronicCompleat2018d.mat w; 

%% Duet 14, 0951a birds ***********************************
% close all; clear all;
% load /Users/eric/ChronicCompleat2018d.mat % OLD DATA
% 
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0951/Duet_Female_2016-01-13T09-51-19a.mat
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0951/Duet_Female_2016-01-13T09-51-19aClicked.mat % 
%     femout = out; clear out
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0951/Duet_Male_2016-01-13T09-51-20_2a.mat
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0951/Duet_Male_2016-01-13T09-51-20_2aClicked.mat % 
%     malout = out; clear out
% 
%     f = 36; m = 35;
%     
%     w(f).wrensex = 'F'; 
%     w(m).wrensex = 'M'; 
%     w(f).sexy = 2; 
%     w(m).sexy = 1;
%     w(f).id = 'j160951a'; 
%     w(m).id = 'm160951a';
%     w(f).Aspikes = {}; %  NO ACUTE DATA
%     w(m).Aspikes = {}; %  NO ACUTE DATA
%     w(f).Fs = 1/Duet_Female_2016_01_13T09_51_19_Ch5.interval; 
%     w(m).Fs = 1/Duet_Female_2016_01_13T09_51_19_Ch5.interval;
%     w(f).duet = Duet_Female_2016_01_13T09_51_19_Ch5.values;
%     w(m).duet = w(f).duet;
%     w(f).tim = 1/w(f).Fs:1/w(f).Fs:Duet_Female_2016_01_13T09_51_19_Ch5.length * Duet_Female_2016_01_13T09_51_19_Ch5.interval;
%     w(m).tim = w(f).tim;
%     
%     % Female chronic spikes
%     w(f).Cspikes{1} = Duet_Female_2016_01_13T09_51_19_Ch11.times;
%     w(f).Cspikes{2} = Duet_Female_2016_01_13T09_51_19_Ch12.times;
%     w(f).Cspikes{3} = Duet_Female_2016_01_13T09_51_19_Ch14.times;
%     w(f).Cspikes{4} = Duet_Female_2016_01_13T09_51_19_Ch15.times;
%     % Male chronic spike
%     w(m).Cspikes{1} = Duet_Male_2016_01_13T09_51_20_2_Ch11.times;
%     w(m).Cspikes{2} = Duet_Male_2016_01_13T09_51_20_2_Ch12.times;
%     w(m).Cspikes{3} = Duet_Male_2016_01_13T09_51_20_2_Ch13.times;
%     w(m).Cspikes{4} = Duet_Male_2016_01_13T09_51_20_2_Ch15.times;
% 
% 
%     % Fix the syllables to 'modern' standards
%     w(f).syl(1).tim(1) = femout(1).syltim(1);
%     w(f).syl(1).tim(2) = femout(3).syltim(2);
% 
%     w(m).syl = w(f).syl;
% 
%     w(f).sylsex = 2;
%     w(m).sylsex = w(f).sylsex;
% 
% % Confirm everything plot
%     figure(2); clf;
%         xa(1) = subplot(411); hold on;
%         bs_raster(w(f).Cspikes, 'm');
%         xa(2) = subplot(412); hold on;
%         bs_raster(w(m).Cspikes, 'b');
%         xa(3) = subplot(413); hold on;        
%         plot(w(f).tim, 2.5+(2.5*w(f).duet/(max(abs(w(m).duet)) )), 'k');
%         xa(4) = subplot(414); hold on;        
%         specgram(w(m).duet, 512, w(f).Fs, [], 500); ylim([500 4500]); colormap('HOT'); caxis([-10 40]);        
%         
%         linkaxes(xa, 'x');
%         
%  for i = 1:length(w(f).syl)
%     for j = 1:3
%         subplot(4,1,j);
%         plot([w(f).syl(i).tim(1) w(f).syl(i).tim(1)], [0 5], 'g');
%         plot([w(m).syl(i).tim(2) w(m).syl(i).tim(2)], [0 5], 'r');
%     end
%     subplot(414);         
%         plot([w(m).syl(i).tim(1) w(f).syl(i).tim(1)], [500 4500], 'g');
%         plot([w(f).syl(i).tim(2) w(m).syl(i).tim(2)], [500 4500], 'r');
% 
%  end
%         
% save /Volumes/ExFAT/Chronically/ChronicCompleat2018d.mat w; 


%% Duet 14, 0951b birds ***********************************
% close all; clear all;
% load /Users/eric/ChronicCompleat2018d.mat % OLD DATA
% 
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0951/Duet_Female_2016-01-13T09-51-19b.mat
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0951/Duet_Female_2016-01-13T09-51-19bClicked.mat % 
%     femout = out; clear out
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0951/Duet_Male_2016-01-13T09-51-20_2b.mat
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0951/Duet_Male_2016-01-13T09-51-20_2bClicked.mat % 
%     malout = out; clear out
% 
%     f = 38; m = 37;
%     
%     w(f).wrensex = 'F'; 
%     w(m).wrensex = 'M'; 
%     w(f).sexy = 2; 
%     w(m).sexy = 1;
%     w(f).id = 'j160951b'; 
%     w(m).id = 'm160951b';
%     w(f).Aspikes = {}; %  NO ACUTE DATA
%     w(m).Aspikes = {}; %  NO ACUTE DATA
%     w(f).Fs = 1/Duet_Female_2016_01_13T09_51_19_Ch5.interval; 
%     w(m).Fs = 1/Duet_Female_2016_01_13T09_51_19_Ch5.interval;
%     w(f).duet = Duet_Female_2016_01_13T09_51_19_Ch5.values;
%     w(m).duet = w(f).duet;
%     w(f).tim = 1/w(f).Fs:1/w(f).Fs:Duet_Female_2016_01_13T09_51_19_Ch5.length * Duet_Female_2016_01_13T09_51_19_Ch5.interval;
%     w(m).tim = w(f).tim;
%     
%     % Female chronic spikes
%     w(f).Cspikes{1} = Duet_Female_2016_01_13T09_51_19_Ch11.times;
%     w(f).Cspikes{2} = Duet_Female_2016_01_13T09_51_19_Ch12.times;
%     w(f).Cspikes{3} = Duet_Female_2016_01_13T09_51_19_Ch14.times;
%     w(f).Cspikes{4} = Duet_Female_2016_01_13T09_51_19_Ch15.times;
%     % Male chronic spike
%     w(m).Cspikes{1} = Duet_Male_2016_01_13T09_51_20_2_Ch11.times;
%     w(m).Cspikes{2} = Duet_Male_2016_01_13T09_51_20_2_Ch12.times;
%     w(m).Cspikes{3} = Duet_Male_2016_01_13T09_51_20_2_Ch13.times;
%     w(m).Cspikes{4} = Duet_Male_2016_01_13T09_51_20_2_Ch15.times;
% 
% 
%     % Fix the syllables to 'modern' standards
%     w(f).syl(1).tim(1) = femout(1).syltim(1);
%     w(f).syl(1).tim(2) = femout(3).syltim(2);
% 
%     w(m).syl = w(f).syl;
% 
%     w(f).sylsex = 2;
%     w(m).sylsex = w(f).sylsex;
% 
% % Confirm everything plot
%     figure(2); clf;
%         xa(1) = subplot(411); hold on;
%         bs_raster(w(f).Cspikes, 'm');
%         xa(2) = subplot(412); hold on;
%         bs_raster(w(m).Cspikes, 'b');
%         xa(3) = subplot(413); hold on;        
%         plot(w(f).tim, 2.5+(2.5*w(f).duet/(max(abs(w(m).duet)) )), 'k');
%         xa(4) = subplot(414); hold on;        
%         specgram(w(m).duet, 512, w(f).Fs, [], 500); ylim([500 4500]); colormap('HOT'); caxis([-10 40]);        
%         
%         linkaxes(xa, 'x');
%         
%  for i = 1:length(w(f).syl)
%     for j = 1:3
%         subplot(4,1,j);
%         plot([w(f).syl(i).tim(1) w(f).syl(i).tim(1)], [0 5], 'g');
%         plot([w(m).syl(i).tim(2) w(m).syl(i).tim(2)], [0 5], 'r');
%     end
%     subplot(414);         
%         plot([w(m).syl(i).tim(1) w(f).syl(i).tim(1)], [500 4500], 'g');
%         plot([w(f).syl(i).tim(2) w(m).syl(i).tim(2)], [500 4500], 'r');
% 
%  end
%         
% save /Volumes/ExFAT/Chronically/ChronicCompleat2018d.mat w; 


%% Duet 20, 0951c birds ***********************************
% close all; clear all;
% load /Users/eric/ChronicCompleat2018d.mat % OLD DATA
% 
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0951/Duet_Female_2016-01-13T09-51-19c.mat
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0951/Duet_Female_2016-01-13T09-51-19cClicked.mat % 
%     femout = out; clear out
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0951/Duet_Male_2016-01-13T09-51-20_2c.mat
% load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0951/Duet_Male_2016-01-13T09-51-20_2cClicked.mat % 
%     malout = out; clear out
% 
%     f = 40; m = 39;
%     
%     w(f).wrensex = 'F'; 
%     w(m).wrensex = 'M'; 
%     w(f).sexy = 2; 
%     w(m).sexy = 1;
%     w(f).id = 'j160951c'; 
%     w(m).id = 'm160951c';
%     w(f).Aspikes = {}; %  NO ACUTE DATA
%     w(m).Aspikes = {}; %  NO ACUTE DATA
%     w(f).Fs = 1/Duet_Female_2016_01_13T09_51_19_Ch5.interval; 
%     w(m).Fs = 1/Duet_Female_2016_01_13T09_51_19_Ch5.interval;
%     w(f).duet = Duet_Female_2016_01_13T09_51_19_Ch5.values;
%     w(m).duet = w(f).duet;
%     w(f).tim = 1/w(f).Fs:1/w(f).Fs:Duet_Female_2016_01_13T09_51_19_Ch5.length * Duet_Female_2016_01_13T09_51_19_Ch5.interval;
%     w(m).tim = w(f).tim;
%     
%     % Female chronic spikes
%     w(f).Cspikes{1} = Duet_Female_2016_01_13T09_51_19_Ch11.times;
%     w(f).Cspikes{2} = Duet_Female_2016_01_13T09_51_19_Ch12.times;
%     w(f).Cspikes{3} = Duet_Female_2016_01_13T09_51_19_Ch14.times;
%     w(f).Cspikes{4} = Duet_Female_2016_01_13T09_51_19_Ch15.times;
%     % Male chronic spike
%     w(m).Cspikes{1} = Duet_Male_2016_01_13T09_51_20_2_Ch11.times;
%     w(m).Cspikes{2} = Duet_Male_2016_01_13T09_51_20_2_Ch12.times;
%     w(m).Cspikes{3} = Duet_Male_2016_01_13T09_51_20_2_Ch13.times;
%     w(m).Cspikes{4} = Duet_Male_2016_01_13T09_51_20_2_Ch15.times;
% 
% 
%     % Fix the syllables to 'modern' standards
%     w(f).syl(1).tim(1) = femout(1).syltim(1);
%     w(f).syl(1).tim(2) = femout(3).syltim(2);
% 
%     w(m).syl = w(f).syl;
% 
%     w(f).sylsex = 2;
%     w(m).sylsex = w(f).sylsex;
% 
% % Confirm everything plot
%     figure(2); clf;
%         xa(1) = subplot(411); hold on;
%         bs_raster(w(f).Cspikes, 'm');
%         xa(2) = subplot(412); hold on;
%         bs_raster(w(m).Cspikes, 'b');
%         xa(3) = subplot(413); hold on;        
%         plot(w(f).tim, 2.5+(2.5*w(f).duet/(max(abs(w(m).duet)) )), 'k');
%         xa(4) = subplot(414); hold on;        
%         specgram(w(m).duet, 512, w(f).Fs, [], 500); ylim([500 4500]); colormap('HOT'); caxis([-10 40]);        
%         
%         linkaxes(xa, 'x');
%         
%  for i = 1:length(w(f).syl)
%     for j = 1:3
%         subplot(4,1,j);
%         plot([w(f).syl(i).tim(1) w(f).syl(i).tim(1)], [0 5], 'g');
%         plot([w(m).syl(i).tim(2) w(m).syl(i).tim(2)], [0 5], 'r');
%     end
%     subplot(414);         
%         plot([w(m).syl(i).tim(1) w(f).syl(i).tim(1)], [500 4500], 'g');
%         plot([w(f).syl(i).tim(2) w(m).syl(i).tim(2)], [500 4500], 'r');
% 
%  end
%         
% save /Volumes/ExFAT/Chronically/ChronicCompleat2018d.mat w; 




%% Duet 21, 0951d birds ***********************************
close all; clear all;
load /Users/eric/ChronicCompleat2018d.mat % OLD DATA

load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0951/Duet_Female_2016-01-13T09-51-19d.mat
load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0951/Duet_Female_2016-01-13T09-51-19dClicked.mat % 
    femout = out; clear out
load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0951/Duet_Male_2016-01-13T09-51-20_2d.mat
load /Users/eric/WrenSourceData/WrenAdd-Aug2018/0951/Duet_Male_2016-01-13T09-51-20_2dClicked.mat % 
    malout = out; clear out

    f = 42; m = 41;
    
    w(f).wrensex = 'F'; 
    w(m).wrensex = 'M'; 
    w(f).sexy = 2; 
    w(m).sexy = 1;
    w(f).id = 'j160951d'; 
    w(m).id = 'm160951d';
    w(f).Aspikes = {}; %  NO ACUTE DATA
    w(m).Aspikes = {}; %  NO ACUTE DATA
    w(f).Fs = 1/Duet_Female_2016_01_13T09_51_19_Ch5.interval; 
    w(m).Fs = 1/Duet_Female_2016_01_13T09_51_19_Ch5.interval;
    w(f).duet = Duet_Female_2016_01_13T09_51_19_Ch5.values;
    w(m).duet = w(f).duet;
    w(f).tim = 1/w(f).Fs:1/w(f).Fs:Duet_Female_2016_01_13T09_51_19_Ch5.length * Duet_Female_2016_01_13T09_51_19_Ch5.interval;
    w(m).tim = w(f).tim;
    
    % Female chronic spikes
    w(f).Cspikes{1} = Duet_Female_2016_01_13T09_51_19_Ch11.times;
    w(f).Cspikes{2} = Duet_Female_2016_01_13T09_51_19_Ch12.times;
    w(f).Cspikes{3} = Duet_Female_2016_01_13T09_51_19_Ch14.times;
    w(f).Cspikes{4} = Duet_Female_2016_01_13T09_51_19_Ch15.times;
    % Male chronic spike
    w(m).Cspikes{1} = Duet_Male_2016_01_13T09_51_20_2_Ch11.times;
    w(m).Cspikes{2} = Duet_Male_2016_01_13T09_51_20_2_Ch12.times;
    w(m).Cspikes{3} = Duet_Male_2016_01_13T09_51_20_2_Ch13.times;
    w(m).Cspikes{4} = Duet_Male_2016_01_13T09_51_20_2_Ch15.times;


    % Fix the syllables to 'modern' standards
    w(f).syl(1).tim(1) = femout(1).syltim(1);
    w(f).syl(1).tim(2) = femout(3).syltim(2);

    w(m).syl = w(f).syl;

    w(f).sylsex = 2;
    w(m).sylsex = w(f).sylsex;

% Confirm everything plot
    figure(2); clf;
        xa(1) = subplot(411); hold on;
        bs_raster(w(f).Cspikes, 'm');
        xa(2) = subplot(412); hold on;
        bs_raster(w(m).Cspikes, 'b');
        xa(3) = subplot(413); hold on;        
        plot(w(f).tim, 2.5+(2.5*w(f).duet/(max(abs(w(m).duet)) )), 'k');
        xa(4) = subplot(414); hold on;        
        specgram(w(m).duet, 512, w(f).Fs, [], 500); ylim([500 4500]); colormap('HOT'); caxis([-10 40]);        
        
        linkaxes(xa, 'x');
        
 for i = 1:length(w(f).syl)
    for j = 1:3
        subplot(4,1,j);
        plot([w(f).syl(i).tim(1) w(f).syl(i).tim(1)], [0 5], 'g');
        plot([w(m).syl(i).tim(2) w(m).syl(i).tim(2)], [0 5], 'r');
    end
    subplot(414);         
        plot([w(m).syl(i).tim(1) w(f).syl(i).tim(1)], [500 4500], 'g');
        plot([w(f).syl(i).tim(2) w(m).syl(i).tim(2)], [500 4500], 'r');

 end
        
% save /Volumes/ExFAT/Chronically/ChronicCompleat2018d.mat w; 





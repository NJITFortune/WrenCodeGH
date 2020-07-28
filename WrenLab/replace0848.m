% Load old data
load /Users/eric/Sync/Wren/ChronicCompleat2019f.mat
% Load the new data from Spike2 for entry 0848 2017 MUTE birds
load ~/Sync/Wren/Mute_Male_2017-01-06T08-48-55_newExport.mat
load ~/Sync/Wren/Mute_Female_2017-01-06T08-48-55_newExport.mat

% Sample rate
Fs = 10000;

% Time stamps for each recording
ftim = 1/Fs:1/Fs:Mute_Female_2017_01_06T08_48_55_export_Ch5.length/Fs;
mtim = 1/Fs:1/Fs:Mute_Male_2017_01_06T08_48_55_export_Ch5.length/Fs;

% Each recording
fsig = Mute_Female_2017_01_06T08_48_55_export_Ch5.values / max(abs(Mute_Female_2017_01_06T08_48_55_export_Ch5.values));
msig = Mute_Male_2017_01_06T08_48_55_export_Ch5.values / max(abs(Mute_Male_2017_01_06T08_48_55_export_Ch5.values));

% These 'offsets' are to synchronize to the previous time base 
% The female is set to match the duet signal in our structure.
% The male was synchronized to the female by lining up the digital pulses.
    femaleoffset = 788.0457;
    maleoffset = 789.6242;
    
% And move the zero to sync with the zero as set in
% ChronicCompleat2019f.mat for w(13).tim
    ftim = ftim - femaleoffset;
    mtim = mtim - maleoffset;

% Old data for copying into updated data    
syl = w(13).syl;
sylsex = w(13).sylsex;

% New syllable boundaries as determined by clicking
    newsyl(1).tim = [-28.8602 -28.5003];
    newsyl(2).tim = [-26.8926 -26.5327];
    newsyl(3).tim = [-24.0852 -23.7253];
    
% Each of the new syllables are solo male.
    newsylsex = [1 1 1];
    


% Update the song sample    
    
w(13).tim = ftim(ftim > -30 & ftim < 10);    
w(14).tim = ftim(ftim > -30 & ftim < 10);    
w(13).duet = fsig(ftim > -30 & ftim < 10);    
w(14).duet = fsig(ftim > -30 & ftim < 10);    
    
% Update the syllable times and identities    

w(13).sylsex = [newsylsex sylsex];
w(14).sylsex = [newsylsex sylsex];

for j=1:3
    w(13).syl(j).tim = newsyl(j).tim;
    w(14).syl(j).tim = newsyl(j).tim;
end

for j=1:length(syl)
    w(13).syl(j+3).tim = syl(j).tim;
    w(14).syl(j+3).tim = syl(j).tim;
end

% Update the spike times

    fCh1 = Mute_Female_2017_01_06T08_48_55_export_Ch11.times;
    fCh2 = Mute_Female_2017_01_06T08_48_55_export_Ch12.times;
    fCh3 = Mute_Female_2017_01_06T08_48_55_export_Ch13.times;
    fCh4 = Mute_Female_2017_01_06T08_48_55_export_Ch14.times;

    mCh1 = Mute_Male_2017_01_06T08_48_55_export_Ch11.times;
    mCh2 = Mute_Male_2017_01_06T08_48_55_export_Ch12.times;
    mCh3 = Mute_Male_2017_01_06T08_48_55_export_Ch13.times;
    mCh4 = Mute_Male_2017_01_06T08_48_55_export_Ch14.times;

w(13).Cspikes{1} = mCh1(mCh1 > -30 & mCh1 < 10);
w(13).Cspikes{2} = mCh1(mCh2 > -30 & mCh2 < 10);
w(13).Cspikes{3} = mCh1(mCh3 > -30 & mCh3 < 10);
w(13).Cspikes{4} = mCh1(mCh4 > -30 & mCh4 < 10);

w(14).Cspikes{1} = fCh1(fCh1 > -30 & fCh1 < 10);
w(14).Cspikes{2} = fCh1(fCh2 > -30 & fCh2 < 10);
w(14).Cspikes{3} = fCh1(fCh3 > -30 & fCh3 < 10);
w(14).Cspikes{4} = fCh1(fCh4 > -30 & fCh4 < 10);


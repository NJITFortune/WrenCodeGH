function [mdata fdata dat tim] = wrenxc(m, f, Fs, td)
% [mdata fdata ampdiff tim] = wrenxx(male, female, Fs, td);
% male and female are the waveforms - should be exactly the same length
% td is the time difference as determined by the taps at the beginning of
% the recording (in Seconds). Assumes that the male is delayed (+ values).
% If the female is delayed, then td is negative.

%% Prepare the signals

% Make sure signals exactly the same length and subtract any offset
    len = min([length(m) length(f)]);
    mdata = m(1:len) - mean(m(1:len)); 
    fdata = f(1:len) - mean(f(1:len));

% Make an initial time base
    tim = 1/Fs:1/Fs:len/Fs;

% Usually the male recording is more retarded, so we cut the front off of it.
if td > 0;
    tt = find(tim > td); 
    mdata = mdata(tt); 
    fdata = fdata(1:length(mdata));
end;

% But sometimes the female recording is more retarded than the male.
if td < 0;
    tt = find(tim > td); 
    fdata = fdata(tt); 
    mdata = mdata(1:length(fdata));
end;

% Make a revised time base
tim = 1/Fs:1/Fs:length(m)/Fs;

%% Get the Amplitude envelope and differences between microphones

% Filter strategy 1: medfilt
%    filtwin = Fs / 50; % 50msec, 20Hz low-pass cutoff
%    mam = medfilt2(abs(m),[filtwin 1]);    
%    maf = medfilt2(abs(f),[filtwin 1]);   

% Filter strategy 2: filtfilt
        low = 20; % 20Hz low pass cutoff
        ordr = 5;
        Wnlow = low*2/Fs;
        [bb,aa] = butter(ordr,Wnlow,'low');
        mam = filtfilt(bb,aa,abs(mdata));
        maf = filtfilt(bb,aa,abs(fdata));
    
    
    dat = mam - maf;

end

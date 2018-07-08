function out = rp(signalA, signalB, Fs)
% out = rp(signalA, signalB, Fs);
% signalA and signalB should be the same length (I think)
% Fs is sample rate

%% stretch the samples to that they are the same length
if length(signalA) ~= length(signalB);
    
    % We are going to stretch to the longer of the two samples, so find
    % the length of the longer signal.
    len = max([length(signalA) length(signalB)]);

    % Make time sequences for each input and the output. I was lazy and
    % just do both signals instead of just the shorter one.  May need to
    % change this.
    timA = 1/Fs:1/Fs:length(signalA)/Fs;
    timB = 1/Fs:1/Fs:length(signalB)/Fs;
    resamp = 1/Fs:1/Fs:len/Fs;

    % Use interp1 to resample both signals.
    signalA = interp1(timA, signalA, resamp);
    signalB = interp1(timB, signalB, resamp);
end;

%% Spectrographic analysis
% Caclulate the spectrogram for the first signal
    [~, ~, ~,op] = spectrogram(signalA,1024,1000,1024,Fs);
    op = 10*log10(abs(op));

% Caclulate the spectrogram for the second signal
    [~,cf,ct,cp] = spectrogram(signalB,1024,1000,1024,Fs);
    cp = 10*log10(abs(cp));

% Make an output structure
    out.spec = cp - op;
    out.f = cf;
    out.t = ct;
    out.a = signalA;
    out.b = signalB;

%% Plotting (maybe remove this)
    mesh(ct,cf,out.spec);
    colormap('HOT');
    view(0,90); 
    ylim([700 7000]);

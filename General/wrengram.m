function wrengram(sig, Fs, cmap)
% Usage: wrengram (signal, Fs, [colormap]);
% colormap is either 1 for HOT or 2 for grayscale.
% wrengram makes a nice spectrogram for Wren songs
% with an embedded oscillogram. 

if nargin == 2; cmap=1; end;

mag = 300; % A scalar for the specgram
cx = [0 50]; % Color Axis map

        foo = sig - mean(sig); % Eliminate DC offset
        mx = max( abs(foo) ); foo = foo * mag/mx; % Normalize to +/- 300
        tim = 1/Fs:1/Fs:length(sig)/Fs; % Compute time series for data

        ax(1) = subplot(2,1,2,'position', [0 0 1 0.2]);
        plot(tim(1:5:end), foo(1:5:end), 'b');
        axis([tim(1) tim(end) -mag mag]);

	ax(2) = subplot(2,1,1,'position', [0 0.2 1 0.8]);
        specgram(sig,1024,Fs,[],1000); 	
	%spectrogram(sig,1024,1000,1024,Fs,'yaxis');

if cmap == 1; colormap('HOT'); end;

if cmap == 2; gg = abs(gray - 1); colormap(gg); end;

linkaxes(ax,'x');
ylim([500 7000]); caxis(cx);


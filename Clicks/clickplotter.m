function [ clicktimes ] = clickplotter( data, Fs, preclick )
% clicktimes = clickplotter(data, Fs);
% data is a chunk of raw wren wrecorded data
% Fs is the sample rate in Hz
% This function plots the chunk of song and has the user click all of the
% syllable starts and ends in order. 
% NOT A STANDALONE FUNCTION
% 
% requires OSCSON

    yvals = [100 4500];
    

%% Plot and get clicks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        figure(1); clf; 
        
		figprop = get(gcf,'Position'); 
        set(gcf,'Position',[figprop(1) figprop(2) 1000 400]);

        oscson(data,Fs);
%        oscsonlow(data,Fs,[-50 60]);

	% Plot prior click in the main window

		if nargin > 2 
                figure(1); 
                hold on;
                plot([preclick(1) preclick(1)], yvals, 'm'); 
                hold off; 
        end

    % get clicks for the starts and ends of the syllables
    
        figure(1); [clicktimes, ~] = ginputc('Color','w','ShowPoints',true,'ConnectPoints',false);        
        
end


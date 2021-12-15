function out = ClickFix(in)
% out = ClickFix(in) 
% This fixes some clicking issues for the distance data.

%% SETUP
out = in; % Copy the input to the output

buff = 0.100; % Buffer before and after the syllable in seconds (0.050 seems appropriate)
pltrng = 50; % This is in % for caxis.


%% Cycle through every duet in the sample (With current plan, just a single duet)
for j = 1:length(in)

% Plot the specgrams for each microphone so that we can track our progress    
    figure(1); clf;
        ax(1) = subplot(211); specgram(in(j).femMic, 2048, in(j).Fs, [], 2000); hold on;
        ax(2) = subplot(212); specgram(in(j).maleMic, 2048, in(j).Fs, [], 2000); hold on; 
        linkaxes(ax, 'xy');
        ylim([200 6000]);

% Plot the autogenous traces across both specgrams
    tt =  find([in(j).fsyl.sexsyltype] > 49); % FEMALE SYLLABLES IN TOP PLOT   
    for k = 1:length(tt)
        subplot(211); plot(in(j).fsyl(tt(k)).traceTim + in(j).fsyl(tt(k)).syltim(1), in(j).fsyl(tt(k)).traceFreq, 'm-', 'LineWidth', 2);
    end
    tt =  find([in(j).msyl.sexsyltype] < 49); % MALE SYLLABLES IN BOTTOM PLOT
    for k = 1:length(tt)
        subplot(212); plot(in(j).msyl(tt(k)).traceTim + in(j).msyl(tt(k)).syltim(1), in(j).msyl(tt(k)).traceFreq, 'b-', 'LineWidth', 2);
    end

%% Cycle through each syllable to make corrections
   for k = 1:length([in(j).fsyl.sexsyltype])
       
       figure(1); % This puts the bars along the top of the sonogram to mark user progress
       if in(j).fsyl(k).sexsyltype > 49
           subplot(211); plot([in(j).fsyl(k).syltim(1), in(j).fsyl(k).syltim(2)], [5500 5500], 'm-', 'LineWidth', 4);
           subplot(212); plot([in(j).msyl(k).syltim(1), in(j).msyl(k).syltim(2)], [5500 5500], 'm-', 'LineWidth', 4);
       else
           subplot(211); plot([in(j).fsyl(k).syltim(1), in(j).fsyl(k).syltim(2)], [5500 5500], 'b-', 'LineWidth', 4);
           subplot(212); plot([in(j).msyl(k).syltim(1), in(j).msyl(k).syltim(2)], [5500 5500], 'b-', 'LineWidth', 4);
       end
       if k > 1
           subplot(211); plot([in(j).fsyl(k-1).syltim(1), in(j).fsyl(k-1).syltim(2)], [5500 5500], 'k-', 'LineWidth', 5);
           subplot(212); plot([in(j).msyl(k-1).syltim(1), in(j).msyl(k-1).syltim(2)], [5500 5500], 'k-', 'LineWidth', 5);
       end



aaa = 9; % Just a dummy to get things started for our while loop...

while aaa > 1 
       
       figure(2); clf;

       longerDur = max([out(j).fsyl(k).sylen, out(j).msyl(k).sylen]);
       
       % Set which2fix: 1 for fixing the heterogenous syllable, 2 for both
%        if in(j).fsyl(k).sexsyltype > 49
%             if in(j).fsyl(k).sylen > in(j).msyl(k).sylen
%                 which2fix = 1;
%             else
%                 which2fix = 2;
%             end           
%        end
%        if in(j).fsyl(k).sexsyltype < 49
%             if in(j).msyl(k).sylen > in(j).fsyl(k).sylen
%                 which2fix = 1;
%             else
%                 which2fix = 2;
%             end           
%        end
       
       figure(2); axx(1) = subplot(211); % FEMALE MICROPHONE
       
       tim = 1/out(j).Fs:1/out(j).Fs:length(out(j).femMic)/out(j).Fs;
       
       specgram(out(j).femMic(tim > out(j).fsyl(k).syltim(1)-buff & tim < out(j).fsyl(k).syltim(1)+longerDur+buff), 2048, in(j).Fs, [], 2000); 
       cx{1} = caxis; caxis([cx{1}(1)*(pltrng/100), floor(cx{1}(2))]); hold on;
       text(0.05, 5500, 'Female Microphone', 'Color', 'w');

            if (out(j).fsyl(k).sexsyltype > 49) % Female syllable
                text(0.05, 5000, 'This is a female syllable', 'Color', 'm');
                plot(out(j).fsyl(k).traceTim+buff, out(j).fsyl(k).traceFreq, 'm-', 'LineWidth', 3); % HEAVY TRACE
                plot([buff, buff], [200 6000], 'g-', 'LineWidth', 1); % START LINE
                plot([buff+out(j).fsyl(k).sylen, buff+out(j).fsyl(k).sylen], [200 6000], 'r-', 'LineWidth', 1); % END LINE
            end
            if (out(j).fsyl(k).sexsyltype < 49) % Male syllable
                text(0.05, 5000, 'Blue is other microphone, Yellow is this microphone', 'Color', 'y');
                plot(out(j).msyl(k).traceTim+buff, out(j).msyl(k).traceFreq, 'b-', 'LineWidth', 2); % Autogenous trace
                plot(out(j).fsyl(k).traceTim+buff, out(j).fsyl(k).traceFreq, 'y-', 'LineWidth', 3); % Heterogenous trace
                plot([buff, buff], [200 6000], 'g-', 'LineWidth', 1); % START LINE
                plot([buff+out(j).fsyl(k).sylen, buff+out(j).fsyl(k).sylen], [200 6000], 'r-', 'LineWidth', 1); % END LINE
            end

       figure(2); axx(2) = subplot(212); % Male Microphone

       tim = 1/out(j).Fs:1/out(j).Fs:length(out(j).maleMic)/out(j).Fs;
            
       specgram(out(j).maleMic(tim > out(j).msyl(k).syltim(1)-buff & tim < out(j).msyl(k).syltim(1)+longerDur+buff), 2048, out(j).Fs, [], 2000); 
       cx{2} = caxis; caxis([cx{2}(1)*(pltrng/100), floor(cx{2}(2))]); hold on;
       text(0.05, 5500, 'Male Microphone', 'Color', 'w');
       
            if (out(j).msyl(k).sexsyltype < 49) % Male syllable
                text(0.05, 5000, 'This is a male syllable', 'Color', 'b');
                plot(out(j).msyl(k).traceTim+buff, out(j).msyl(k).traceFreq, 'b-', 'LineWidth', 2);
                plot([buff, buff], [200 6000], 'g-', 'LineWidth', 1); % START LINE 
                plot([buff+out(j).msyl(k).sylen, buff+out(j).msyl(k).sylen], [200 6000], 'r-', 'LineWidth', 1); % END LINE
            end
            if (out(j).msyl(k).sexsyltype > 49) % Female syllable
                text(0.05, 5000, 'Magenta is other microphone, Yellow is this microphone', 'Color', 'y');
                plot(out(j).fsyl(k).traceTim+buff, out(j).fsyl(k).traceFreq, 'm-', 'LineWidth', 2); % Autogenous trace
                plot(out(j).msyl(k).traceTim+buff, out(j).msyl(k).traceFreq, 'y-', 'LineWidth', 3); % Heterogenous trace
                plot([buff, buff], [200 6000], 'g-', 'LineWidth', 1); % START LINE
                plot([buff+out(j).msyl(k).sylen, buff+out(j).msyl(k).sylen], [200 6000], 'r-', 'LineWidth', 1); % END LINE
            end
 
        linkaxes(axx, 'xy');
        
        figure(2); subplot(211); colormap('GRAY'); ylim([200 6000]); 
        pause(0.1);
        grid(axx, "on")
        axx(1).GridColor = [0 0.9 1];
        axx(1).GridAlpha = 0.5;
        axx(2).GridColor = [0 0.9 1];
        axx(2).GridAlpha = 0.5;

    
    fprintf('1 FIX end only and move to next syllable. \n');
    fprintf('2 get one click for start of heterogenous and review. \n');
    fprintf('3 reclick both (3 clicks, start and end of autogenous, start of heterogenous) and review. \n');
    fprintf('10 to 90 change contrast (lower is higher). \n');
    aaa = input('ACTION: ');

    if aaa >= 10
        pltrng = aaa;
    end
    
    if aaa == 0 % Extend heterogenous evenly in both directions
        if (out(j).fsyl(k).sexsyltype > 49) % Female is autogenous
            amt = (out(j).fsyl(k).sylen - out(j).msyl(k).sylen)/2;
            out(j).msyl(k).sylen = out(j).fsyl(k).sylen;
            out(j).msyl(k).syltim(1) = out(j).msyl(k).syltim(1)-amt;
            out(j).msyl(k).syltim(2) = out(j).msyl(k).syltim(1) + out(j).fsyl(k).sylen;
            out(j).msyl(k).sylidx(1) = round(out(j).msyl(k).sylidx(1) - (amt * out(j).Fs));
            out(j).msyl(k).sylidx(2) = out(j).msyl(k).sylidx(1) + (out(j).fsyl(k).sylidx(2) - out(j).fsyl(k).sylidx(1));

            tmp = syldat(out(j).maleMic(out(j).msyl(k).sylidx(1):out(j).msyl(k).sylidx(2)), out(j).Fs);
            out(j).msyl(k).traceTim = tmp.trace_tim;
            out(j).msyl(k).traceFreq = tmp.trace_freq;
            out(j).msyl(k).trace.peakfreq = tmp.trace_peakf;
            out(j).msyl(k).trace.slopemean = tmp.trace_slopemean;
            out(j).msyl(k).trace.slopestd = tmp.trace_slopestd;
            out(j).msyl(k).trace.slopevar = tmp.trace_slopevar;
            
        end
        if (out(j).fsyl(k).sexsyltype < 49) % Male is autogenous
            amt = (out(j).msyl(k).sylen - out(j).fsyl(k).sylen)/2;
            out(j).fsyl(k).sylen = out(j).msyl(k).sylen;
            out(j).fsyl(k).syltim(1) = out(j).fsyl(k).syltim(1)-amt;
            out(j).fsyl(k).syltim(2) = out(j).fsyl(k).syltim(1) + out(j).msyl(k).sylen;
            out(j).fsyl(k).sylidx(1) = round(out(j).fsyl(k).sylidx(1) - (amt * out(j).Fs));
            out(j).fsyl(k).sylidx(2) = out(j).fsyl(k).sylidx(1) + (out(j).msyl(k).sylidx(2) - out(j).msyl(k).sylidx(1));

            tmp = syldat(out(j).femMic(out(j).fsyl(k).sylidx(1):out(j).fsyl(k).sylidx(2)), out(j).Fs);
            out(j).fsyl(k).traceTim = tmp.trace_tim;
            out(j).fsyl(k).traceFreq = tmp.trace_freq;
            out(j).fsyl(k).trace.peakfreq = tmp.trace_peakf;
            out(j).fsyl(k).trace.slopemean = tmp.trace_slopemean;
            out(j).fsyl(k).trace.slopestd = tmp.trace_slopestd;
            out(j).fsyl(k).trace.slopevar = tmp.trace_slopevar;        
        end
    end
    
    if aaa == 1 % Extend the end of heterogenous
        if (out(j).fsyl(k).sexsyltype > 49) % Female is autogenous
            out(j).msyl(k).sylen = out(j).fsyl(k).sylen;
            out(j).msyl(k).syltim(2) = out(j).msyl(k).syltim(1) + out(j).fsyl(k).sylen;
            out(j).msyl(k).sylidx(2) = out(j).msyl(k).sylidx(1) + (out(j).fsyl(k).sylidx(2) - out(j).fsyl(k).sylidx(1));

            tmp = syldat(out(j).maleMic(out(j).msyl(k).sylidx(1):out(j).msyl(k).sylidx(2)), out(j).Fs);
            out(j).msyl(k).traceTim = tmp.trace_tim;
            out(j).msyl(k).traceFreq = tmp.trace_freq;
            out(j).msyl(k).trace.peakfreq = tmp.trace_peakf;
            out(j).msyl(k).trace.slopemean = tmp.trace_slopemean;
            out(j).msyl(k).trace.slopestd = tmp.trace_slopestd;
            out(j).msyl(k).trace.slopevar = tmp.trace_slopevar;
        end
        if (out(j).fsyl(k).sexsyltype < 49) % Male is autogenous
            out(j).fsyl(k).sylen = out(j).msyl(k).sylen;
            out(j).fsyl(k).syltim(2) = out(j).fsyl(k).syltim(1) + out(j).msyl(k).sylen;
            out(j).fsyl(k).sylidx(2) = out(j).fsyl(k).sylidx(1) + (out(j).msyl(k).sylidx(2) - out(j).msyl(k).sylidx(1));

            tmp = syldat(out(j).femMic(out(j).fsyl(k).sylidx(1):out(j).fsyl(k).sylidx(2)), out(j).Fs);
            out(j).fsyl(k).traceTim = tmp.trace_tim;
            out(j).fsyl(k).traceFreq = tmp.trace_freq;
            out(j).fsyl(k).trace.peakfreq = tmp.trace_peakf;
            out(j).fsyl(k).trace.slopemean = tmp.trace_slopemean;
            out(j).fsyl(k).trace.slopestd = tmp.trace_slopestd;
            out(j).fsyl(k).trace.slopevar = tmp.trace_slopevar;        
        end
     end
    
    if aaa == 2 % Get single new click for heterogenous at start
        if (out(j).fsyl(k).sexsyltype > 49) % Female is autogenous
            fprintf('Click on start syllable at the male microphone.\n');
            out(j).msyl(k).sylen = out(j).fsyl(k).sylen;
            figure(2); [newclk, ~] = ginput(1); 
            out(j).msyl(k).syltim(1) = out(j).msyl(k).syltim(1) + (newclk - buff);
            out(j).msyl(k).syltim(2) = out(j).msyl(k).syltim(1) + out(j).fsyl(k).sylen;
            out(j).msyl(k).sylidx(1) = round(out(j).msyl(k).sylidx(1) + ((newclk - buff) * out(j).Fs));
            out(j).msyl(k).sylidx(2) = out(j).msyl(k).sylidx(1) + (out(j).fsyl(k).sylidx(2) - out(j).fsyl(k).sylidx(1));

            tmp = syldat(out(j).maleMic(out(j).msyl(k).sylidx(1):out(j).msyl(k).sylidx(2)), out(j).Fs);
            out(j).msyl(k).traceTim = tmp.trace_tim;
            out(j).msyl(k).traceFreq = tmp.trace_freq;
            out(j).msyl(k).trace.peakfreq = tmp.trace_peakf;
            out(j).msyl(k).trace.slopemean = tmp.trace_slopemean;
            out(j).msyl(k).trace.slopestd = tmp.trace_slopestd;
            out(j).msyl(k).trace.slopevar = tmp.trace_slopevar;            
        end
        if (out(j).fsyl(k).sexsyltype < 49) % Male is autogenous
            fprintf('Click on start syllable at the female microphone.\n');
            out(j).fsyl(k).sylen = out(j).msyl(k).sylen;
            figure(2); [newclk, ~] = ginput(1); 
            out(j).fsyl(k).syltim(1) = out(j).fsyl(k).syltim(1) + (newclk - buff);
            out(j).fsyl(k).syltim(2) = out(j).fsyl(k).syltim(1) + out(j).msyl(k).sylen;
            out(j).fsyl(k).sylidx(1) = round(out(j).fsyl(k).sylidx(1) + ((newclk - buff) * out(j).Fs));
            out(j).fsyl(k).sylidx(2) = out(j).fsyl(k).sylidx(1) + (out(j).msyl(k).sylidx(2) - out(j).msyl(k).sylidx(1));

            tmp = syldat(out(j).femMic(out(j).fsyl(k).sylidx(1):out(j).fsyl(k).sylidx(2)), out(j).Fs);
            out(j).fsyl(k).traceTim = tmp.trace_tim;
            out(j).fsyl(k).traceFreq = tmp.trace_freq;
            out(j).fsyl(k).trace.peakfreq = tmp.trace_peakf;
            out(j).fsyl(k).trace.slopemean = tmp.trace_slopemean;
            out(j).fsyl(k).trace.slopestd = tmp.trace_slopestd;
            out(j).fsyl(k).trace.slopevar = tmp.trace_slopevar;        
        end
    end
    
    if aaa == 3 % Get new clicks for autogenous and start click for heterogenous

        if (out(j).fsyl(k).sexsyltype > 49) % Female is autogenous
            fprintf('Click on start and end of female syllable, female microphone.\n');
            figure(2); [newclx, ~] = ginput(2); newclx = sort(newclx);
                out(j).fsyl(k).sylen = newclx(2) - newclx(1);
                out(j).fsyl(k).syltim(1) = out(j).fsyl(k).syltim(1) + (newclx(1) - buff);
                out(j).fsyl(k).syltim(2) = out(j).fsyl(k).syltim(1) + out(j).fsyl(k).sylen;
                out(j).fsyl(k).sylidx(1) = out(j).fsyl(k).sylidx(1) + round(((newclx(1) - buff) * out(j).Fs));
                out(j).fsyl(k).sylidx(2) = out(j).fsyl(k).sylidx(1) + round((out(j).fsyl(k).sylen * out(j).Fs));
            
            tmp = syldat(out(j).femMic(out(j).fsyl(k).sylidx(1):out(j).fsyl(k).sylidx(2)), out(j).Fs);
                out(j).fsyl(k).traceTim = tmp.trace_tim;
                out(j).fsyl(k).traceFreq = tmp.trace_freq;
                out(j).fsyl(k).trace.peakfreq = tmp.trace_peakf;
                out(j).fsyl(k).trace.slopemean = tmp.trace_slopemean;
                out(j).fsyl(k).trace.slopestd = tmp.trace_slopestd;
                out(j).fsyl(k).trace.slopevar = tmp.trace_slopevar;        

            fprintf('Click on start syllable at the male microphone.\n');
            figure(2); [newclk, ~] = ginput(1);
                out(j).msyl(k).sylen = out(j).fsyl(k).sylen;
                out(j).msyl(k).syltim(1) = out(j).msyl(k).syltim(1) + (newclk - buff);
                out(j).msyl(k).syltim(2) = out(j).msyl(k).syltim(1) + out(j).fsyl(k).sylen;
                out(j).msyl(k).sylidx(1) = out(j).msyl(k).sylidx(1) + round(((newclk - buff) * out(j).Fs));
                out(j).msyl(k).sylidx(2) = out(j).msyl(k).sylidx(1) + round((out(j).fsyl(k).sylen * out(j).Fs));

            tmp = syldat(out(j).maleMic(out(j).msyl(k).sylidx(1):out(j).msyl(k).sylidx(2)), out(j).Fs);
                out(j).msyl(k).traceTim = tmp.trace_tim;
                out(j).msyl(k).traceFreq = tmp.trace_freq;
                out(j).msyl(k).trace.peakfreq = tmp.trace_peakf;
                out(j).msyl(k).trace.slopemean = tmp.trace_slopemean;
                out(j).msyl(k).trace.slopestd = tmp.trace_slopestd;
                out(j).msyl(k).trace.slopevar = tmp.trace_slopevar;            
        end
        
        if (out(j).fsyl(k).sexsyltype < 49) % Male is autogenous
            fprintf('Click on start and end of male syllable, male microphone.\n');
            figure(2); [newclx, ~] = ginput(2); newclx = sort(newclx);
                out(j).msyl(k).sylen = newclx(2) - newclx(1);
                out(j).msyl(k).syltim(1) = out(j).msyl(k).syltim(1) + (newclx(1) - buff);
                out(j).msyl(k).syltim(2) = out(j).msyl(k).syltim(1) + out(j).msyl(k).sylen;
                out(j).msyl(k).sylidx(1) = out(j).msyl(k).sylidx(1) + round(((newclx(1) - buff) * out(j).Fs));
                out(j).msyl(k).sylidx(2) = out(j).msyl(k).sylidx(1) + round((out(j).msyl(k).sylen * out(j).Fs));
            
            tmp = syldat(out(j).maleMic(out(j).msyl(k).sylidx(1):out(j).msyl(k).sylidx(2)), out(j).Fs);
                out(j).msyl(k).traceTim = tmp.trace_tim;
                out(j).msyl(k).traceFreq = tmp.trace_freq;
                out(j).msyl(k).trace.peakfreq = tmp.trace_peakf;
                out(j).msyl(k).trace.slopemean = tmp.trace_slopemean;
                out(j).msyl(k).trace.slopestd = tmp.trace_slopestd;
                out(j).msyl(k).trace.slopevar = tmp.trace_slopevar;            

            fprintf('Click on start syllable at the female microphone.\n');
            figure(2); [newclk, ~] = ginput(1);
                out(j).fsyl(k).sylen = out(j).msyl(k).sylen;
                out(j).fsyl(k).syltim(1) = out(j).fsyl(k).syltim(1) + (newclk - buff);
                out(j).fsyl(k).syltim(2) = out(j).fsyl(k).syltim(1) + out(j).msyl(k).sylen;
                out(j).fsyl(k).sylidx(1) = out(j).fsyl(k).sylidx(1) + round(((newclk - buff) * out(j).Fs));
                out(j).fsyl(k).sylidx(2) = out(j).fsyl(k).sylidx(1) + round((out(j).msyl(k).sylen * out(j).Fs));

            tmp = syldat(out(j).femMic(out(j).fsyl(k).sylidx(1):out(j).fsyl(k).sylidx(2)), out(j).Fs);
                out(j).fsyl(k).traceTim = tmp.trace_tim;
                out(j).fsyl(k).traceFreq = tmp.trace_freq;
                out(j).fsyl(k).trace.peakfreq = tmp.trace_peakf;
                out(j).fsyl(k).trace.slopemean = tmp.trace_slopemean;
                out(j).fsyl(k).trace.slopestd = tmp.trace_slopestd;
                out(j).fsyl(k).trace.slopevar = tmp.trace_slopevar;        
        end
           
    end
   
end % replot and re-ask if clicked.
    
   end % End of syllable loop

end % End of loop by duet

end % End of main function



function foo = syldat(sng, Fs)

% Trace the frequency of the syllable
	siglen = length(sng);
	tinytim = 1/Fs:1/Fs:siglen/Fs;
    stp = 0.005; % Step size
	startwin = 0;
	p = 1;
    nfft = 512; % FFT size
    
	while startwin < max(tinytim) + stp
            % Matlab removed the find command here
        	ttdata = sng(tinytim > startwin & tinytim < startwin + stp);
        	L = length(ttdata);

        	t_fftdata = fft(ttdata,nfft)/L;
        	t_fftdata = 2 * abs(t_fftdata(1:nfft/2+1));
        	t_freqs = Fs/2*linspace(0,1,nfft/2+1);

        	[maxval, idx] = max(t_fftdata);
        	peakamp(p) = maxval;
        	peakfreq(p) = t_freqs(idx);
        	stw(p) = startwin + stp/2;

        	p = p + 1; startwin = startwin + stp/2;
	end

	foo.trace_freq = medfilt1(peakfreq(1:end-2),5);
    foo.trace_amp = medfilt1(peakamp(1:end-2),5);
    foo.trace_tim = stw(1:end-2);
    
	foo.trace_slopemean = mean(diff(foo.trace_freq));
	foo.trace_slopestd = std(diff(foo.trace_freq));
	foo.trace_slopevar = var(diff(foo.trace_freq));
	[foo.trace_peakf, indf] = max(foo.trace_freq);
	foo.trace_peakt = stw(indf);
	foo.trace_peakpt = foo.trace_peakt / stw(end);
end
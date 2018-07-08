function out = groupsyllables(struct);

%band                  freqtracepeakpercent  loud                  slopemean             spectrum              
%fftfreqs              freqtracepeaktim      maxfreq               slopestd              syl                   
%freqtim               Fs                    meanloud              slopevar              sylen                 
%freqtrace             ind                   minfreq               song                  time                  
%freqtracepeak         ISI                   peakfreq              specfilt              

close all;

figure(1); plot(struct.syl(:,1), struct.sylen, '*-');



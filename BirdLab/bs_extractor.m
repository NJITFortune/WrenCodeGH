function [burst_spikes,iso_spikes] = bs_extractor( spike_times, isi_thresh )
% [s_times, b_times, b_spikes] = burst_calc( spike_times, isi_thresh )
% VERSION 18 August 2008
% finds bursting events in a list of SPIKE_TIMES,
% given an interspike interval threshold ISI_THRESH
%
% returns a list of single-spike times, burst-event times,
% and a cell array in which each element is an array with
% the times of the spikes within a single burst
%
% JAB 7/14/07
% Modded by ESF for the cell array data

reps = length(spike_times);
isi_thresh = isi_thresh/1000;

% loop for each repetition in the cell array of spike times

for i=1:reps
    
    sptimes = ( spike_times{i} );
    isis = diff( sptimes ); % ISIs
    good_diffs = find ( isis <= isi_thresh );
    len=length(good_diffs);
    
    burstspikes = [];
    keeper=1;
    
    for getit=1:(len-2)
        burstspikes( keeper ) = sptimes( good_diffs(getit) );
        keeper=keeper+1;
        
        if ( (sptimes( (good_diffs(getit)+1) ) - sptimes( good_diffs(getit) ) ) < isi_thresh )
            if ( sptimes( (good_diffs(getit)+2) ) - sptimes( (good_diffs(getit)+1) ) > isi_thresh )
                
                burstspikes( keeper ) = sptimes( (good_diffs(getit)+1) );
                keeper=keeper+1;
            end
        end
    end
    
    
    burst_spikes{i} = burstspikes';
    isos = setdiff( sptimes, burstspikes );
    iso_spikes{i} = isos(:);
    
    clear sptimes burstspikes;
    
end


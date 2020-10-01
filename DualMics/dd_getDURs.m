function out = dd_getDURs(in)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fsyldurs = [];
msyldurs = [];

for d = 1:length(idx)
    
    numsyls = length(in(idx(d)).fsyl);
    
    % Cycle through every syllable except the last
    for s=1:numsyls

        if in(idx(d)).fsyl(s).sexsyltype < 49 
            fsyldurs(end+1) = in(idx(d)).fsyl(s).sylen;
        end
        if in(idx(d)).msyl(s).sexsyltype > 49 
            msyldurs(end+1) = in(idx(d)).msyl(s).sylen;
        end
        
    end
    
    
    
end


out.fsyldurs = fsyldurs;
out.msyldurs = msyldurs;


end


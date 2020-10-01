function out = dd_getDURs(in)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fsyldurs = [];
msyldurs = [];

for d = 1:length(in)
    
    numsyls = length(in(d).fsyl);
    
    % Cycle through every syllable except the last
    for s=1:numsyls

        if in(d).msyl(s).sexsyltype < 49 
            fsyldurs(end+1) = in(d).msyl(s).sylen;
        end
        if in(d).fsyl(s).sexsyltype > 49 
            msyldurs(end+1) = in(d).msyl(s).sylen;
        end
        
    end
    
    
    
end


out.fsyldurs = fsyldurs;
out.msyldurs = msyldurs;


end


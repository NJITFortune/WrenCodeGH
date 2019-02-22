function [vector_strength, phasespikes] = wPhasor(spiketimes, tims)
    
    for j = length(tims):-1:1
        
    prespikes = pi * (spiketimes(spiketimes > tims(j,1) & spiketimes < tims(j,2)) / (tims(j,2) - tims(j,1)));
    postspikes = pi + (pi * (spiketimes(spiketimes > tims(j,2) & spiketimes < tims(j,3)) / (tims(j,3) - tims(j,3))) );

    phasespikes{j} = sort([prespikes postspikes]);
    vector_strength(j) = sqrt(mean(cos(phasespikes{j})).^2 + mean(sin(phasespikes{j})).^2);

    end

figure(27); clf; bs_raster(phasespikes);
    
    
end

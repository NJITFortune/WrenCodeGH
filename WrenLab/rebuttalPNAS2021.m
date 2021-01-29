function rebuttalPNAS2021

[ww, Fs] = audioread('C:\Users\eric_fortune\Downloads\wrenagade.wav');

cmap = [0 33]; 
nfft = [1024 512 1024 512];

plotit(ww, Fs, cmap, nfft);

maxtim = length(ww)/Fs;
    newFs = 10000;    
    newtim = 1/newFs:1/newFs:maxtim;
    tim = 1/Fs:1/Fs:maxtim;
neww = resample(ww, tim, newFs); 

plotit(neww, newFs, cmap, nfft);



function plotit(data, plotFs, colrmap, fftsize)

figure;
clrmappping = flipud(gray);
overlap = 0.95;
fntsz = 12;

subplot(411);
    specgram(data, fftsize(1), plotFs, [], ceil(fftsize(1)*overlap)); 
    colormap(clrmappping);
    ylim([500 5000]);
    caxis(colrmap);
    txtxt = ['Fs=' num2str(plotFs) ', nFFT=' num2str(fftsize(1)) ', coloraxis=[' num2str(colrmap) ']'];    
    text(0.2, 4400, txtxt, 'FontSize', fntsz, 'Color', 'r'); 
subplot(412);
    specgram(data, fftsize(2), plotFs, [], ceil(fftsize(2)*overlap)); 
    colormap(clrmappping);
    ylim([500 5000]);
    caxis(colrmap);
    txtxt = ['Fs=' num2str(plotFs) ', nFFT=' num2str(fftsize(2)) ', coloraxis=[' num2str(colrmap) ']'];    
    text(0.2, 4400, txtxt, 'FontSize', fntsz, 'Color', 'r');
subplot(413);
    specgram(data, fftsize(3), plotFs, [], ceil(fftsize(3)*overlap)); 
    colormap(clrmappping);
    ylim([500 5000]);
    caxis(colrmap);
    txtxt = ['Fs=' num2str(plotFs) ', nFFT=' num2str(fftsize(3)) ', coloraxis=[' num2str(colrmap) ']'];    
    text(0.2, 4400, txtxt, 'FontSize', fntsz, 'Color', 'r');
subplot(414);
    specgram(data, fftsize(4), plotFs, [], ceil(fftsize(4)*overlap)); 
    colormap(clrmappping);
    ylim([500 5000]);
    caxis(colrmap);
    txtxt = ['Fs=' num2str(plotFs) ', nFFT=' num2str(fftsize(4)) ', coloraxis=[' num2str(colrmap) ']'];    
    text(0.2, 4400, txtxt, 'FontSize', fntsz, 'Color', 'r');
        

end

end

% load /Users/eric/Sync/Wren/DistanceData-V10.mat


j = 0;


while j < length(dd)
    
fprintf('Last j was %i.  ', j); j = input('What next? ');

    figure(3);clf; 
        xx(1) = subplot(211); specgram(dd(j).femMic, 1024, dd(j).Fs); 
        xx(2) = subplot(212); specgram(dd(j).maleMic, 1024, dd(j).Fs); 
        linkaxes(xx,'xy'); ylim([0 8000]);
    drawnow;

dat = [num2str(dd(j).day) '_' num2str(dd(j).month) '_' num2str(dd(j).year)];

fprintf('Date is: %s \n', dat);

k = input('Enter a number if bad: ');
    if isempty(k)
    loca = input('Enter Location. ');
    nam = input('Enter the pair unique four character name: ');
    uniquestr = input('Enter time in ms from midnight! ');
    
malefilename =   [nam 'Male' num2str(dd(j).distance)   'm_' dat '_' loca '_' num2str(uniquestr) '.wav'];
femalefilename = [nam 'Female' num2str(dd(j).distance) 'm_' dat '_' loca '_' num2str(uniquestr) '.wav'];


fprintf('Writing %s and %s ...\n', malefilename, femalefilename);

audiowrite(malefilename, dd(j).maleMic, dd(j).Fs);
audiowrite(femalefilename, dd(j).femMic, dd(j).Fs);

fprintf(' \n');
    end
end



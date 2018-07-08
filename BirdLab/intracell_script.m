%%
load dual_hvc1

keys = dual_hvc1_Ch31;

stim = dual_hvc1_Ch3;

spk_lft = dual_hvc1_Ch401;
spk_rht = dual_hvc1_Ch402;

labels = {'BOSD','BOSU1','BOSU2','CON','REV'};
raster = bs_converter(spk_lft,keys,stim,[10 8], labels);
%raster = bs_converter(spk_lft,keys,stim,[0.2 15]);



Fs = raster.SampleRate;

%%

win = [-raster(1).pretrig raster(1).posttrig];

icellBOS = bs


outBOS = bs_calculator(raster(1).Spikes, raster(1).Stim, Fs, win);


%%
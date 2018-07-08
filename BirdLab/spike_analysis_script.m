%%
load dual_hvc1

keys = dual_hvc1_Ch31;

stim = dual_hvc1_Ch3;

spk_lft = dual_hvc1_Ch401;
spk_rht = dual_hvc1_Ch402;

labels = {'BOSD','BOSU1','BOSU2','CON','REV'};

raster = bs_converter(spikes,keys,stim,[6 6], labels, psp);

Fs = raster.SampleRate;%raster = bs_converter(spk_lft,keys,stim,[0.2 15]);



%%
win = [-raster(1).pretrig raster(1).posttrig];



out = bs_calculator(raster(1).Spikes, raster(1).Stim, Fs, win);



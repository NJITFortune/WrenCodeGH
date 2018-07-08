function [song] = songeditrev(filename);
%form: [song] = wavedit(filename,savedir);
%
%example: bos1 = wavedit('3_August_12_2006_15_05_21.wav','c:/songs');
%
%This function loads a wave file (filename), asks the user to select the 'song' out of
%the wave file, gets the user to name the wave file, and then filters and saves the
%file in savedir.  Default savedir is current directory.
%made/modified TAN 8/17/06 from filtdatc,wavedit

fs=44100;

format short g
format compact

%savedir='c:/editedsongs';

f=findstr(filename,'wav');
g=findstr(filename,'dat');
if ~isempty(f)
    [x,fs,nbits]=wavread(filename);
elseif ~isempty(g)
    [a]=loaddatvi(filename);
    x=a(5,:);
else
    error('The filename must contain the extension wav or dat.')
end;
specgram(x,[],fs);
v=axis;
figure(1)
clf
subplot(2,1,1)
plot(x)%plot the unfiltered waveform in blue
hold on

hold on
%make filter to filter song 500 to 10000 Hz (Hamming fixed impulse
%response)
order=256;
b=fir1(order,[500/(fs/2) 10000/(fs/2)]);
xf=conv(x,b);
xf2=xf(order/2:length(xf)-order/2-1);


figure(1)
plot(xf2,'r');%plot the filtered waveform in red
axis tight
subplot(2,1,2)
specgram(xf2,[],fs)%plot the sonogram of the filtered waveform
v=axis;
axis([v(1) v(2) 0 10000])

%normalize
song=xf2./max(abs(xf2));%normalize to max
song=song*0.999;

%get user to choose song and save
input('Using the magnifier tool, select the song on the red trace and hit enter.\n Note only x axis matters.\n')

v=axis;
b=song(round(v(1)):round(v(2)));
sound(b,44100);
nm=input('Enter the name of the song.\n Example: bu70b3\n','s');
rnm=[nm 'r.wav'];
nm=[nm '.wav'];

wavwrite(b,44100,16,nm);

r=reverse(b);
wavwrite(r,44100,16,rnm);


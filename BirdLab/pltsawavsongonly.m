function pltsawavsongonly
%form: pltsawavs(firstletterofmonth)
%
%To begin, combine all files of interest (must be 4 at same time point) into one 
%directory. YOu can use combinesafolders.m
%first letter of month should be capitalized in form 'J'

format compact;
format short g;



%find all trials
ext=['_']

dfiles = dir;
files = [];

%  get files of type ext....
for i = 1:length(dfiles)
   if ~isempty(findstr(dfiles(i).name, ext))
      files = [files i];
   end;
end;

files = dfiles(files);

%****************************************************************************************************************************
numsongs=length(files);
for i = 1:length(files)

   sprintf('%d %% done',round(100*i/numsongs))
    f1=files(i).name
   %f2=['2' f1(2:length(f1))];
   % f3=['3' f1(2:length(f1))];
   % f4=['4' f1(2:length(f1))]

    x1 = wavread(f1);
   % x2=wavread(f2);
    %x3=wavread(f3);
    %x4=wavread(f4);

        %make song recording prettier
        x1=diff(x1);

    
    figure(1)
    clf
    
    %subplot(4,1,1)
    specgram(x1,[],44100)
    axis tight
    
   % subplot(4,1,2)
   % plot(x2,'k')
   % axis tight
    
    %subplot(4,1,3)
   % plot(x3,'k')
   % axis tight

   %subplot(4,1,4)
   % plot(x4,'k')
   % axis tight
input('Hit enter to move to the next song.\n');

end;

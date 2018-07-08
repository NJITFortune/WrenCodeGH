function pltsawavsongonly_sid_noah(first_file)
%form: pltsawavs(firstletterofmonth)
%
%To begin, combine all files of interest (must be 4 at same time point) into one
%directory. YOu can use combinesafolders.m
%first letter of month should be capitalized in form 'J'

format compact;
format short g;




%find all trials
ext='_';

files = dir(['*' ext '*']);

if nargin == 0
    first_file = files(1).name;
end

reached_first = 0;

stop_code = 0;
%****************************************************************************************************************************
numsongs=length(files);

% for i = 1:length(files)
i = 1;

show_menu();

while ~stop_code
    
    sprintf('%d %% done',round(100*i/numsongs))
    f1=files(i).name;
    disp(f1)
    
    if strcmp(f1,first_file)
        reached_first = 1;
    elseif reached_first==0
        i = i+1;
    end
    if reached_first
        
        [x1, Fs] = wavread(f1);
        x1 = x1/max(abs(x1));

        
        %make song recording prettier
        x1=diff(x1);
        
        
        figure(1)
        clf
        
        specgram(x1,[],44100)
        axis tight
        
        disp('Hit n to move to the next song.\n');
        [~,~,user_inp] = ginput(1);
        if isempty(user_inp)
            user_inp = 'n';
        end
        if user_inp == 'p'
            user_inp = 'b';
        end
        switch user_inp
            case 's'
                disp('Playing sound....')
                sound(x1,Fs);
                disp('... done!');
            case 'q'
                stop_code = 1;
            case 'b'
                i = i-1;
            case 'e'
                songeditrev(f1);
            case 'j'
                i = input('Enter song number -- numerical value only');
                max([i 0]);
                i = round(min([max([i 0]); numsongs]));
            otherwise
                i = i+1;
        end
    end
    if i > numsongs
        disp('That was the last song!\n')
        stop_code=1;
    end
end;


function  show_menu();
disp('n: next song\n')
disp('b/p: back/previous song\n')
disp('e: edit song\n')
disp('s: sound play (listen to) song\n')
disp('q: quit\n')


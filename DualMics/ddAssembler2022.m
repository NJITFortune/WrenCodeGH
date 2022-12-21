% This script assembles the NEW and IMPROVED Distance Data structure.
%
% 
%% OLDE (pairnum = 1)  
% This is a Bellavista pair from November 2012 with both vision and no vision

cd olde

% ZERO M
idx = 1;
    load olde0m14232153A.mat
    dd(idx) = foo;
    dd(idx).vision = 1; % Vision
    dd(idx).pairnum = 1;
    clear foo;

idx = 2;
    load olde0m14232156A.mat
        foo.pairnum = 1;
    dd(idx) = foo;
    dd(idx).vision = 1; % Vision
    dd(idx).pairnum = 1;
    clear foo;

idx = 3;
    load olde0m14232157A.mat
        foo.pairnum = 1;
    dd(idx) = foo;
    dd(idx).vision = 1; % Vision
    clear foo;
    
% ONE M
% 8am - 9:05am Vision
% 9:05am and later No Vision

idx = 4;
    load olde1m14151234A.mat
        foo.pairnum = 1;
    dd(idx) = foo;
    dd(idx).vision = 9; % Not yet known
    clear foo;
    
idx = 5;
    load olde1m14232170A.mat
        out.pairnum = 1;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;

idx = 6;
    load olde1m14152103A.mat
        foo.pairnum = 1;
    dd(idx) = foo;
    dd(idx).vision = 9; % Not yet known
    clear foo;

idx = 7;
    load olde1m14232171A.mat
        out.pairnum = 1;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;
    
idx = 8;
    load olde1m14232172A.mat
        out.pairnum = 1;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;    
    
idx = 9;
    load olde1m14232173A.mat
        out.pairnum = 1;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;

idx = 10;
    load olde1m14232174A.mat
        out.pairnum = 1;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;
    
% THREE M
% 9:40am Vision
% 12:00pm No Vision

idx = 11;
    load olde3m14152117A.mat
        foo.pairnum = 1;
    dd(idx) = foo;
    dd(idx).vision = 9; % Not yet known
    clear foo;
    
idx = 12;
    load olde3m14152121A.mat
        foo.pairnum = 1;
    dd(idx) = foo;
    dd(idx).vision = 9; % Not yet known
    clear foo;

idx = 13;
    load olde3m14232181A.mat
        out.pairnum = 1;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;

% FIVE M
% 12:17pm Vision
% 2:47pm No Vision

idx = 14;
    load olde5m14152125A.mat
        foo.pairnum = 1;
    dd(idx) = foo;
    dd(idx).vision = 9; % Not yet known
    clear foo;
    
idx = 15;
    load olde5m14152126A.mat
        foo.pairnum = 1;
    dd(idx) = foo;
    dd(idx).vision = 9; % Not yet known
    clear foo;

idx = 16;
    load olde5m14232192A.mat
        foo.pairnum = 1;
    dd(idx) = foo;
    dd(idx).vision = 9; % Not yet known
    clear foo;    
    
idx = 17;
    load olde5m32222205A.mat
        out.pairnum = 1;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;

% SEVEN M
% All No Vision

idx = 18;
    load olde7m14152146A.mat
        foo.pairnum = 1;
    dd(idx) = foo;
    dd(idx).vision = 0; % No Vision
    clear foo;
   
idx = 19;
    load olde7m14152148A.mat 
        foo.pairnum = 1;
    dd(idx) = foo;
    dd(idx).vision = 0; % No Vision
    clear foo;
    
idx = 20;
    load olde7m14232199A.mat
        foo.pairnum = 1;
    dd(idx) = foo;
    dd(idx).vision = 0; % No Vision
    clear foo;
    
idx = 21;
    load olde7m32222230A.mat
        out.pairnum = 1;
    dd(idx) = out;
    dd(idx).vision = 0; % No Vision
    clear out;

idx = 22;
    load olde7m32222231A.mat
        out.pairnum = 1;
    dd(idx) = out;
    dd(idx).vision = 0; % No Vision
    clear out;
    
idx = 23;
    load olde7m32222232A.mat
        out.pairnum = 1;
    dd(idx) = out;
    dd(idx).vision = 0; % No Vision
    clear out;    
    
%% PULL (pairnum = 2)  
% This is a 2016 Yanayacu pair 5-January-2016 - BEND BIRDS
% All vision
cd ../pull

% ZERO M

idx = 24;
    load pull0m23452345A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;    

idx = 25;
    load pull0m23452355A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;    
        
idx = 26;
    load pull0m23452365A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;    

idx = 27;
    load pull0m23452377A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;    

% TWO M

idx = 28;
    load pull2m23452383A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;    

idx = 29;
    load pull2m23452389A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;    

idx = 30;
    load pull2m23452391A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;    

% pull2m23452391AX.mat

idx = 31;
    load pull2m23452393A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;   

idx = 32;
    load pull2m23452395A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;    

idx = 33;
    load pull2m23452398A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;    

idx = 34;
    load pull2m23452399A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;    
    
idx = 35;
    load pull2m23452532A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;    
    
idx = 36;
    load pull2m23452535A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;    
    
idx = 37;
    load pull2m23452537A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;    

idx = 38;
    load pull2m23452538A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;    

idx = 39;
    load pull2m23452612A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;  

idx = 40;
    load pull2m23452622A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;  

idx = 41;
    load pull2m23452624A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;  

idx = 42;
    load pull2m23452625A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;  

idx = 42;
    load pull2m23452627A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out; 

idx = 43;
    load pull2m23452628A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out; 

idx = 44;
    load pull2m23452630A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out; 

idx = 45;
    load pull2m23452633A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out; 

idx = 46;
    load pull2m23452634A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out; 

% THREE M

idx = 47;
    load pull3m23452648A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out; 

idx = 48;
    load pull3m23452650A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out; 

idx = 49;
    load pull3m23452651A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out; 

idx = 50;
    load pull3m23452655A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out; 

idx = 51;
    load pull3m23452656A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out; 

% NINE M

idx = 52;
    load pull9m23452657A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out; 

idx = 53;
    load pull9m23452658A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out; 

idx = 54;
    load pull9m23452659A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out; 

idx = 55;
    load pull9m23452661A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out; 

idx = 56;
    load pull9m23452666A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out; 

idx = 57;
    load pull9m23452667A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision 
    clear out; 

idx = 58;
    load pull9m23452669A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision 
    clear out; 

idx = 59;
    load pull9m23452677A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision 
    clear out;

% pull9m23452678A.mat Error that needs fixing as of Sept. 2022

idx = 60;
    load pull9m23452680A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision 
    clear out;

% TEN M

idx = 61;
    load pull10m23453333A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision 
    clear out;

%% SIXT (pairnum = 3)  
% This is Bellavista 26-November-2012
cd ../sixt
% THREE M
% Vision unknown from notes

idx = 62;
    load sixt3m32222240A.mat
        out.pairnum = 3;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;

% FIVE M
% 11am No Vision

idx = 63;
    load sixt5m32222233A.mat
        out.pairnum = 3;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;

idx = 64;
    load sixt5m32222234A.mat
        out.pairnum = 3;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;

% SEVEN M
% 2:15pm No Vision
% 4:15pm Vision
 
idx = 65;
    load sixt7m32222235A.mat
        out.pairnum = 3;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;

idx = 66;
    load sixt7m32222236A.mat
        out.pairnum = 3;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;

idx = 67;
    load sixt7m32222237A.mat
        out.pairnum = 3;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;

idx = 68;
    load sixt7m32222238A.mat
        out.pairnum = 3;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;

idx = 69;
    load sixt7m32222241A.mat
        out.pairnum = 3;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;

idx = 70;
    load sixt7m32222242A.mat
        out.pairnum = 3;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;


%% YEGE (pairnum = 4)  
% Siempre Verde 8-January-2014 9-January-2014
% 9-January Vision up to 9:35am 3,5,7,0
% 9-January No Vision 2,5 10am, 11am
% 9-January Vision 2m 5pm
cd ../yege
% ZERO M

idx = 71;
    load yege0m43221001A.mat
        foo.pairnum = 4;
    dd(idx) = foo;
    dd(idx).vision = 9; % Not yet known
    clear foo;

% TWO M

idx = 72;
    load yege2m43211219A.mat
        foo.pairnum = 4;
    dd(idx) = foo;
    dd(idx).vision = 9; % Not yet known
    clear foo;


% FIVE M

idx = 73;
    load yege5m44421658A.mat
        foo.pairnum = 4;
    dd(idx) = foo;
    dd(idx).vision = 9; % Not yet known
    clear foo;



%% POST (pairnum = 5)  
% This is I HAVE NO IDEA
% 24-August-2016

cd ../post
% THREE M

idx = 74;
    load post3m41137938A.mat
        out.pairnum = 5;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out;

idx = 75;
    load post3m41165559A.mat
        out.pairnum = 5;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out;

% idx = 74;
%     load post3m42523714A.mat
%         out.pairnum = 5;
%     dd(idx) = out;
%     dd(idx).vision = 1; % Vision
%     clear out;

% FIVE M

idx = 76;
    load post5m25997390A.mat
        out.pairnum = 5;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out;

idx = 77;
    load post5m26553835A.mat
        out.pairnum = 5;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out;

idx = 78;
    load post5m28848247A.mat
        out.pairnum = 5;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out;

%% End of script

cd ..


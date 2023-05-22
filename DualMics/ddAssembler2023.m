% This script assembles the NEW and IMPROVED Distance Data structure.
%
% 
%% OLDE (pairnum = 1)  
% This is a Bellavista pair from November 2012 with both vision and no vision
% Red Male / Blue Female

cd olde

% ZERO M (Vision Approved)
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
    load olde1m14151234A.mat % ZOOM0002_1m_male_vision.mp3
        foo.pairnum = 1;
    dd(idx) = foo;
    dd(idx).vision = 1; % Update May 19 2023
    clear foo;
    
idx = 5;
    load olde1m14232170A.mat % Assuming ZOOM0002_1m_male_vision.mp3
        out.pairnum = 1;
    dd(idx) = out;
    dd(idx).vision = 1; % Update May 19 2023
    clear out;

idx = 6;
    load olde1m14152103A.mat % ZOOM0002_1m_male_vision.mp3
        foo.pairnum = 1;
    dd(idx) = foo;
    dd(idx).vision = 1; % Update May 19 2023
    clear foo;

idx = 7;
    load olde1m14232171A.mat % Assuming ZOOM0002_1m_male_vision.mp3
        out.pairnum = 1;
    dd(idx) = out;
    dd(idx).vision = 1; % Not yet known
    clear out;
    
idx = 8;
    load olde1m14232172A.mat % zoom0003_1m_male_noVision.mp3 (23rd!)
        out.pairnum = 1;
    dd(idx) = out;
    dd(idx).vision = 0; % Update May 19 2023
    clear out;    
    
idx = 9; 
    load olde1m14232173A.mat % zoom0003_1m_male_noVision.mp3 (23rd!)
        out.pairnum = 1;
    dd(idx) = out;
    dd(idx).vision = 0; % Update May 19 2023
    clear out;

idx = 10;
    load olde1m14232174A.mat % zoom0003_1m_male_noVision.mp3 (23rd!)
        out.pairnum = 1;
    dd(idx) = out;
    dd(idx).vision = 0; % Update May 19 2023
    clear out;
    
% THREE M
% 9:40am Vision
% 12:00pm No Vision

idx = 11;
    load olde3m14152117A.mat % zoom0006_3m_male_vision.mp3 (22nd)
        foo.pairnum = 1;
    dd(idx) = foo;
    dd(idx).vision = 1; % Update May 19 2023
    clear foo;
    
idx = 12;
    load olde3m14152121A.mat % zoom0006_3m_male_vision.mp3 (22nd)
        foo.pairnum = 1;
    dd(idx) = foo;
    dd(idx).vision = 1; % Update May 19 2023
    clear foo;

idx = 13;
    load olde3m14232181A.mat % Assuming zoom0006_3m_male_vision.mp3 (22nd)
        out.pairnum = 1;
    dd(idx) = out;
    dd(idx).vision = 1; % Update May 19 2023
    clear out;

%%%%%%%% NEED TO ADD NOVISION CONDITION ZOOM0006_3m_male_noVision.mp3
%%%%%%%% Added below Indices 108, 109, 110


% FIVE M
% 12:17pm Vision
% 2:47pm No Vision

idx = 14;
    load olde5m14152125A.mat % ZOOM0007_5m_male_vision
        foo.pairnum = 1;
    dd(idx) = foo;
    dd(idx).vision = 1; % Update May 19 2023
    clear foo;
    
idx = 15;
    load olde5m14152126A.mat % ZOOM0007_5m_male_vision
        foo.pairnum = 1;
    dd(idx) = foo;
    dd(idx).vision = 1; % Update May 19 2023
    clear foo;

idx = 16;
    load olde5m14232192A.mat % ZOOM0008_5m_male_noVision
        foo.pairnum = 1;
    dd(idx) = foo;
    dd(idx).vision = 0; % Update May 19 2023
    clear foo;    
    
idx = 17;
    load olde5m32222205A.mat % Looks like from ZOOM0008_5m_male_noVision
        out.pairnum = 1;
    dd(idx) = out;
    dd(idx).vision = 0; % Assumed - looks consistant but didn't find
    clear out;

%%%%%%%% May NEED TO ADD NOVISION CONDITION ZOOM0008_5m_male_noVision.mp3
%%%%%%%% Added another 5m NoVision for entry 111

% SEVEN M

idx = 18;
    load olde7m14152146A.mat % ZOOM0008_7m_male_vision  22nd!!!
        foo.pairnum = 1;
    dd(idx) = foo;
    dd(idx).vision = 1; % Vision
    clear foo;
   
%%%%%%%% Entry 18 is the only full usable duet from ZOOM0008_7m_male_vision.mp3  

idx = 19;
    load olde7m14152148A.mat % ZOOM0009_7m_male_noVision 22nd!!!
        foo.pairnum = 1;
    dd(idx) = foo;
    dd(idx).vision = 0; % No Vision
    clear foo;
    
idx = 20;
    load olde7m14232199A.mat % ZOOM0009_7m_male_noVision 23rd!!!
        foo.pairnum = 1;
    dd(idx) = foo;
    dd(idx).vision = 0; % No Vision
    clear foo;
    
idx = 21;
    load olde7m32222230A.mat % ZOOM0007_7m_male_noVision 24th!!!
        out.pairnum = 1;
    dd(idx) = out;
    dd(idx).vision = 0; % No Vision
    clear out;

idx = 22;
    load olde7m32222231A.mat % ZOOM0007_7m_male_noVision 24th!!!
        out.pairnum = 1;
    dd(idx) = out;
    dd(idx).vision = 0; % No Vision
    clear out;
    
idx = 23;
    load olde7m32222232A.mat % ZOOM0007_7m_male_noVision 24th!!!
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
    dd(idx).vision = 1; % Vision
    clear out;    

idx = 25;
    load pull0m23452355A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out;    
        
idx = 26;
    load pull0m23452365A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out;    

idx = 27;
    load pull0m23452377A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out;    

% TWO M

idx = 28;
    load pull2m23452383A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out;    

idx = 29;
    load pull2m23452389A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out;    

idx = 30;
    load pull2m23452391A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out;    

% pull2m23452391AX.mat

idx = 31;
    load pull2m23452393A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out;   

idx = 32;
    load pull2m23452395A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out;    

idx = 33;
    load pull2m23452398A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out;    

idx = 34;
    load pull2m23452399A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out;    
    
idx = 35;
    load pull2m23452532A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out;    
    
idx = 36;
    load pull2m23452535A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out;    
    
idx = 37;
    load pull2m23452537A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out;    

idx = 38;
    load pull2m23452538A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out;    

idx = 39;
    load pull2m23452612A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out;  

idx = 40;
    load pull2m23452622A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out;  

idx = 41;
    load pull2m23452624A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out;  

idx = 42;
    load pull2m23452625A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out;  

idx = 43;
    load pull2m23452627A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out; 

idx = 44;
    load pull2m23452628A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out; 

idx = 45;
    load pull2m23452630A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out; 

idx = 46;
    load pull2m23452633A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out; 

idx = 47;
    load pull2m23452634A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out; 

% THREE M

idx = 48;
    load pull3m23452648A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out; 

idx = 49;
    load pull3m23452650A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out; 

idx = 50;
    load pull3m23452651A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out; 

idx = 51;
    load pull3m23452655A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out; 

idx = 52;
    load pull3m23452656A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out; 

% NINE M

idx = 53;
    load pull9m23452657A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out; 

idx = 54;
    load pull9m23452658A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out; 

idx = 55;
    load pull9m23452659A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out; 

idx = 56;
    load pull9m23452661A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out; 

idx = 57;
    load pull9m23452666A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision
    clear out; 

idx = 58;
    load pull9m23452667A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision 
    clear out; 

idx = 59;
    load pull9m23452669A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision 
    clear out; 

idx = 60;
    load pull9m23452677A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision 
    clear out;

idx = 61;
    load pull9m23452680A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision 
    clear out;

% TEN M

idx = 62;
    load pull10m23453333A.mat
        out.pairnum = 2;
    dd(idx) = out;
    dd(idx).vision = 1; % Vision 
    clear out;

%% SIXT (pairnum = 3)  
% This is Bellavista 26-November-2012
% Green Male / Blue Female 

cd ../sixt
% THREE M
% Vision unknown from notes [27th No Vision] [22nd with Vision]



% FIVE M
% 11am No Vision [26th No vision] [22nd vision]

idx = 62; %% This is properly labeled in dataset as 5m
    load sixt3m32222240A.mat % ZOOM0001_5m_maleGreen_noVision.mp3
        out.pairnum = 3;
    dd(idx) = out;
    dd(idx).vision = 0; % See below
    clear out;

%%%% AS of 19-May-2023 entry 62 and 63 are duplicates!!!!

idx = 63;
    load sixt5m32222233A.mat % ZOOM0001_5m_maleGreen_noVision.mp3
        out.pairnum = 3;
    dd(idx) = out;
    dd(idx).vision = 0; % See above
    clear out;

idx = 64;
    load sixt5m32222234A.mat % ZOOM0001_5m_maleGreen_noVision.mp3
        out.pairnum = 3;
    dd(idx) = out;
    dd(idx).vision = 0; % Updated 19 May 2023 - No Vision
    clear out;

% SEVEN M
% 2:15pm No Vision 4:15pm Vision [26th]
% Vision 10:45am, NoVision 1pm [27th]
% Vision 4pm, NoVision 5:15pm [22th]

idx = 65;
    load sixt7m32222235A.mat % ZOOM0003_7m_maleGreen_vision
        out.pairnum = 3;
    dd(idx) = out;
    dd(idx).vision = 1; % Updated 19 May 2023
    clear out;

idx = 66;
    load sixt7m32222236A.mat % ZOOM0003_7m_maleGreen_vision
        out.pairnum = 3;
    dd(idx) = out;
    dd(idx).vision = 1; % Updated 19 May 2023
    clear out;

idx = 67;
    load sixt7m32222237A.mat % ZOOM0002_7m_maleGreen_noVision.mp3
        out.pairnum = 3;
    dd(idx) = out;
    dd(idx).vision = 0; % Updated 19 May 2023
    clear out;

idx = 68;
    load sixt7m32222238A.mat % ZOOM0002_7m_maleGreen_noVision.mp3
        out.pairnum = 3;
    dd(idx) = out;
    dd(idx).vision = 0; % Updated 19 May 2023
    clear out;

idx = 69;
    load sixt7m32222241A.mat % CAN'T seem to find this entry
        out.pairnum = 3;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;

idx = 70;
    load sixt7m32222242A.mat % CAN'T seem to find this entry
        out.pairnum = 3;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;


%% YEGE (pairnum = 4)  
% Siempre Verde 8-January-2014 9-January-2014 (may be listed in places as 2016)
% 7-Jan 2m 0m all vision
% 8-Jan no notes assume vision
% 9-January Vision up to 9:35am 3,5,7,2,0
% 9-January No Vision 2, 5 10am, 11am
% 9-January Vision 2m 5pm

%%%%%%%%% 19 May 2023... ADD DATA  Drive Simpre_Verde_1.8.14_distance

cd ../yege
% ZERO M

idx = 71;
    load yege0m43221001A.mat % ZOOM0002_greenMale_0m.mp3
        foo.pairnum = 4;
    dd(idx) = foo;
    dd(idx).vision = 1; % Doesn't say in audio or notes, but probably vision
    clear foo;

% TWO M

idx = 72;
    load yege2m43211219A.mat % ZOOM0001_greenMale_2m.mp3 8-January
        foo.pairnum = 4;
    dd(idx) = foo;
    dd(idx).vision = 1; % Doesn't say in audio or notes, but probably vision
    clear foo;

% FIVE M

idx = 73;
    load yege5m44421658A.mat % ZOOM0008_greenMale_5m_noVision.mp3
        foo.pairnum = 4;
    dd(idx) = foo;
    dd(idx).vision = 0; % Updated 19 May 2023
    clear foo;

%%%%%%%%%%%%%%% NEED TO ADD 2M and 5M both vision and no vision
%%%%%%%%%%%%%%% Added 2M below as entries 112-114
%%%%%%%%%%%%%%% Added 5M below as entries 115-117


%% POST (pairnum = 5)  
% 24-August-2016
% ALL WITH VISION

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

%% CARL (pairnum = 6)  
% 27-July-2012
% ALL WITH VISION

cd ../carl
% TWO M

idx = 79;
    load carl2m_A.mat
        foo.pairnum = 6;
    dd(idx) = foo;
    dd(idx).vision = 1; % Vision
    clear foo;

% FOUR M

idx = 80;
    load carl4m_A10_44.mat
        foo.pairnum = 6;
    dd(idx) = foo;
    dd(idx).vision = 1; % Vision
    clear foo;

% FIVE M

idx = 81;
    load carl5m_A11_26.mat
        foo.pairnum = 6;
    dd(idx) = foo;
    dd(idx).vision = 1; %  Vision
    clear foo;

% SIX M

idx = 82;
    load carl6m_A11_32.mat
        foo.pairnum = 6;
    dd(idx) = foo;
    dd(idx).vision = 1; %  Vision
    clear foo;

idx = 83;
    load carl6m_A11_46.mat
        foo.pairnum = 6;
    dd(idx) = foo;
    dd(idx).vision = 1; %  Vision
    clear foo;

% SEVEN M

idx = 84;
    load carl7m_A12_12.mat
        foo.pairnum = 6;
    dd(idx) = foo;
    dd(idx).vision = 1; %  Vision
    clear foo;

%% BELA (pairnum = 7)  
% 25-January-2013
% ALL WITH VISION

cd ../bela

% ZERO M
idx = 85;
    load bela0m_1a_vision-Fixed.mat
        foo.pairnum = 7;
        foo.timestamp = 101;
    dd(idx) = foo;
    dd(idx).vision = 1; % Vision
    clear foo;

idx = 86;
    load bela0m_1b_vision-Fixed.mat
        foo.pairnum = 7;
        foo.timestamp = 102;
    dd(idx) = foo;
    dd(idx).vision = 1; % Vision
    clear foo;

idx = 87;
    load bela0m_1c_vision-Fixed.mat
        foo.pairnum = 7;
        foo.timestamp = 103;
    dd(idx) = foo;
    dd(idx).vision = 1; % Vision
    clear foo;

idx = 88;
    load bela0m_1d_vision-Fixed.mat
        foo.pairnum = 7;
        foo.timestamp = 104;
    dd(idx) = foo;
    dd(idx).vision = 1; % Vision
    clear foo;

idx = 89;
    load bela0m_2a_vision-Fixed.mat
        foo.pairnum = 7;
        foo.timestamp = 105;
    dd(idx) = foo;
    dd(idx).vision = 1; % Vision
    clear foo;

idx = 90;
    load bela5m_3a_vision-Fixed.mat
        foo.pairnum = 7;
        foo.timestamp = 121;
    dd(idx) = foo;
    dd(idx).vision = 1; % Vision
    clear foo;

idx = 91;
    load bela5m_3b_vision-Fixed.mat
        foo.pairnum = 7;
        foo.timestamp = 122;
    dd(idx) = foo;
    dd(idx).vision = 1; % Vision
    clear foo;

idx = 92;
    load bela5m_3c_vision-Fixed.mat
        foo.pairnum = 7;
        foo.timestamp = 123;
    dd(idx) = foo;
    dd(idx).vision = 1; % Vision
    clear foo;

idx = 93;
    load bela5m_3d_vision-Fixed.mat
        foo.pairnum = 7;
        foo.timestamp = 124;
    dd(idx) = foo;
    dd(idx).vision = 1; % Vision
    clear foo;

idx = 94;
    load bela5m_4_vision-Fixed.mat
        foo.pairnum = 7;
        foo.timestamp = 125;
    dd(idx) = foo;
    dd(idx).vision = 1; % Vision
    clear foo;

idx = 95;
    load bela5m_5a_vision-Fixed.mat
        foo.pairnum = 7;
        foo.timestamp = 126;
    dd(idx) = foo;
    dd(idx).vision = 1; % Vision
    clear foo;

idx = 96;
    load bela5m_5b_vision-Fixed.mat
        foo.pairnum = 7;
        foo.timestamp = 127;
    dd(idx) = foo;
    dd(idx).vision = 1; % Vision
    clear foo;

idx = 97;
    load bela5m_5c_vision-Fixed.mat
        foo.pairnum = 7;
        foo.timestamp = 129;
    dd(idx) = foo;
    dd(idx).vision = 1; % Vision
    clear foo;

idx = 98;
    load bela1m_3_vision_Bellavista-esf.mat
        foo.pairnum = 7;
        foo.timestamp = 130;
    dd(idx) = foo;
    dd(idx).vision = 1; % Vision
    clear foo;    

idx = 99;
    load bela1m_4_vision_Bellavista-esf.mat
        foo.pairnum = 7;
        foo.timestamp = 131;
    dd(idx) = foo;
    dd(idx).vision = 1; % Vision
    clear foo;    
    
idx = 100;
    load bela1m_8a_vision_Bellavista-esf.mat
        foo.pairnum = 7;
        foo.timestamp = 132;
    dd(idx) = foo;
    dd(idx).vision = 1; % Vision
    clear foo;    

idx = 101;
    load bela7m_1_vision_Bellavista2013esf.mat
        foo.pairnum = 7;
        foo.timestamp = 133;
    dd(idx) = foo;
    dd(idx).vision = 1; % Vision
    clear foo;    

idx = 102;
    load bela7m_2_vision_Bellavista2013esf.mat
        foo.pairnum = 7;
        foo.timestamp = 134;
    dd(idx) = foo;
    dd(idx).vision = 1; % Vision
    clear foo;    

idx = 103;
    load bela7m_3_vision_Bellavista2013esf.mat
        foo.pairnum = 7;
        foo.timestamp = 135;
    dd(idx) = foo;
    dd(idx).vision = 1; % Vision
    clear foo;    

idx = 104;
    load bela7m_8_vision_Bellavista2013esf.mat
        foo.pairnum = 7;
        foo.timestamp = 136;
    dd(idx) = foo;
    dd(idx).vision = 1; % Vision
    clear foo;    

idx = 105;
    load bela7m_25_1_2013_Bellavista_9aVision.mat
        foo.pairnum = 7;
        foo.timestamp = 137;
    dd(idx) = foo;
    dd(idx).vision = 1; % Vision
    clear foo;    

idx = 106;
    load bela7m_25_1_2013_Bellavista_9bVision.mat
        foo.pairnum = 7;
        foo.timestamp = 138;
    dd(idx) = foo;
    dd(idx).vision = 1; % Vision
    clear foo;    

idx = 107;
    load bela7m_25_1_2013_Bellavista_9cVision.mat
        foo.pairnum = 7;
        foo.timestamp = 139;
    dd(idx) = foo;
    dd(idx).vision = 1; % Vision
    clear foo;    




%% Additions

% OLDE ZOOM0006_3m_male_noVision.mp3

cd ../olde

idx = 108;
    load olde3m743-fixed.mat
        foo.pairnum = 1;
        foo.timestamp = 743;
    dd(idx) = foo;
    dd(idx).vision = 0; % No Vision
    clear foo;  

idx = 109;
    load olde3m14710-fixed.mat
        foo.pairnum = 1;
        foo.timestamp = 14710;
    dd(idx) = foo;
    dd(idx).vision = 0; % No Vision
    clear foo;   

idx = 110;
    load olde3m15000-fixed.mat
        foo.pairnum = 1;
        foo.timestamp = 15000;
    dd(idx) = foo;
    dd(idx).vision = 0; % No Vision
    clear foo;  

idx = 111;
    load olde5m5749-fixed.mat
        foo.pairnum = 1;
        foo.timestamp = 5749;
    dd(idx) = foo;
    dd(idx).vision = 0; % No Vision
    clear foo;  

cd ../yege

idx = 112;
    load yege2m173050vision-fixed.mat
        foo.pairnum = 1;
        foo.timestamp = 5749;
    dd(idx) = foo;
    dd(idx).vision = 1; % Vision
    clear foo;

idx = 113;
    load yege2m322625novision-fixed.mat
        foo.pairnum = 1;
        foo.timestamp = 5749;
    dd(idx) = foo;
    dd(idx).vision = 0; % No Vision
    clear foo;

idx = 114;
    load yege2m431020novision-fixed.mat
        foo.pairnum = 1;
        foo.timestamp = 5749;
    dd(idx) = foo;
    dd(idx).vision = 0; % No Vision
    clear foo;

idx = 115;
    load yege5m1925vision-fixedfixed.mat
        foo.pairnum = 1;
        foo.timestamp = 1925;
    dd(idx) = foo;
    dd(idx).vision = 1; % Vision
    clear foo;

idx = 116;
    load yege5m1201vision-fixed.mat
        foo.pairnum = 1;
        foo.timestamp = 1201;
    dd(idx) = foo;
    dd(idx).vision = 1; % Vision
    clear foo;

idx = 117;
    load yege5m13900novision-fixed.mat
        foo.pairnum = 1;
        foo.timestamp = 1201;
    dd(idx) = foo;
    dd(idx).vision = 0; % No Vision
    clear foo;
    
%% End of script

cd ..


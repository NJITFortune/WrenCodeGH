% This script assembles the NEW and IMPROVED Distance Data structure.
%
% 
%% OLDE (pairnum = 1)  
% This is a Bellavista pair from November 2012 with both vision and no vision

idx = 1;
    load olde0m14232153A.mat
    dd(idx) = foo;
    dd(idx).vision = 9; % Not yet known
    dd(idx).pairnum = 1;
    clear foo;

idx = 2;
    load olde0m14232157A.mat
        foo.pairnum = 1;
    dd(idx) = foo;
    dd(idx).vision = 9; % Not yet known
    clear foo;
    
idx = 3;
    load olde1m14151234A.mat
        foo.pairnum = 1;
    dd(idx) = foo;
    dd(idx).vision = 9; % Not yet known
    clear foo;
    
idx = 4;
    load olde1m14232170A.mat
        out.pairnum = 1;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;

idx = 5;
    load olde1m14232171A.mat
        out.pairnum = 1;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;
    
idx = 6;
    load olde1m14232172A.mat
        out.pairnum = 1;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;    
    
idx = 7;
    load olde1m14232173A.mat
        out.pairnum = 1;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;

idx = 8;
    load olde1m14232174A.mat
        out.pairnum = 1;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;
    
idx = 9;
    load olde3m14152117A.mat
        foo.pairnum = 1;
    dd(idx) = foo;
    dd(idx).vision = 9; % Not yet known
    clear foo;
    
idx = 10;
    load olde3m14232181A.mat
        out.pairnum = 1;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;

idx = 11;
    load olde5m14152125A.mat
        foo.pairnum = 1;
    dd(idx) = foo;
    dd(idx).vision = 9; % Not yet known
    clear foo;
    
idx = 12;
    load olde5m14232192A.mat
        foo.pairnum = 1;
    dd(idx) = foo;
    dd(idx).vision = 9; % Not yet known
    clear foo;    
    
idx = 13;
    load olde5m32222205A.mat
        out.pairnum = 1;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;

idx = 14;
    load olde7m14152146A.mat
        foo.pairnum = 1;
    dd(idx) = foo;
    dd(idx).vision = 9; % Not yet known
    clear foo;
    
idx = 15;
    load olde7m14232199A.mat
        foo.pairnum = 1;
    dd(idx) = foo;
    dd(idx).vision = 9; % Not yet known
    clear foo;
    
idx = 16;
    load olde7m32222230A.mat
        out.pairnum = 1;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;

idx = 17;
    load olde7m32222231A.mat
        out.pairnum = 1;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;
    
idx = 18;
    load olde7m32222232A.mat
        out.pairnum = 1;
    dd(idx) = out;
    dd(idx).vision = 9; % Not yet known
    clear out;    
    

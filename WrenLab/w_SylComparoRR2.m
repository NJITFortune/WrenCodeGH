%function w_SylComparoRR2
%% m17 Duet 1 (w(1) male and w(2) female) Male solo (1)

% wrengram(w(1).duet, w(1).Fs, 1);
% caxis([-50, 50]);

% Male
    solo{1} = [1 2 3 4];
      Ams(1) =  getMeanSpikes(w, 1, solo{1});    
      Hms(1) =  getMeanSpikes(w, 2, solo{1});    
    solo{2} = 5;
      Ams(2) =  getMeanSpikes(w, 1, solo{2});    
      Hms(2) =  getMeanSpikes(w, 2, solo{2});    

    duet{1} = [7 11];
      Amd(1) =  getMeanSpikes(w, 1, duet{1});    
      Hmd(1) =  getMeanSpikes(w, 2, duet{1});    
    duet{2} = [9 13];
      Amd(2) =  getMeanSpikes(w, 1, duet{2});    
      Hmd(2) =  getMeanSpikes(w, 2, duet{2});    

 fprintf('Male Autogenous solo1 %2.4f duet1 %2.4f \n', Ams(1), Amd(1));
 fprintf('Female Heterogenous solo1 %2.4f duet1 %2.4f \n', Hms(1), Hmd(1));
 fprintf('Male Autogenous solo2 %2.4f duet2 %2.4f \n', Ams(2), Amd(2));
 fprintf('Female Heterogenous solo2 %2.4f duet2 %2.4f \n', Hms(2), Hmd(2));

 
%% j160806 Duet 2 (w(3) male and w(4) female) Female solo (1)

% wrengram(w(3).duet, w(3).Fs, 1);
% caxis([-50, 50]);
clear solo duet

% Female
    solo = 1;
      Afs(1) =  getMeanSpikes(w, 4, solo);    
      Hfs(1) =  getMeanSpikes(w, 3, solo);    

    duet = [5 9 13];
      Afd(1) =  getMeanSpikes(w, 4, duet);    
      Hfd(1) =  getMeanSpikes(w, 3, duet);    
    
 fprintf('Female Autogenous solo1 %2.4f duet1 %2.4f \n', Afs(1), Afd(1));
 fprintf('Male Heterogenous solo1 %2.4f duet1 %2.4f \n', Hfs(1), Hfd(1));
    
%% j160807 Duet 3 (w(5) male and w(6) female) Male and Female solo (1) and (1)

% wrengram(w(5).duet, w(5).Fs, 1);
% caxis([-50, 50]);
clear solo duet

% Male
    solo = 2;
      Ams(3) =  getMeanSpikes(w, 5, solo);    
      Hms(3) =  getMeanSpikes(w, 6, solo);    
    duet = [6 10 14 18];
      Amd(3) =  getMeanSpikes(w, 5, duet);    
      Hmd(3) =  getMeanSpikes(w, 6, duet);        
        
% Female
    solo = 1;
      Afs(2) =  getMeanSpikes(w, 6, solo);    
      Hfs(2) =  getMeanSpikes(w, 5, solo);    
    duet = [3 7 11 15];
      Afd(2) =  getMeanSpikes(w, 6, duet);    
      Hfd(2) =  getMeanSpikes(w, 5, duet);    
    
fprintf('Male Autogenous solo3 %2.4f duet1 %2.4f \n', Ams(3), Amd(3));
fprintf('Female Heterogenous solo3 %2.4f duet1 %2.4f \n', Hms(3), Hmd(3));
fprintf('Female Autogenous solo2 %2.4f duet2 %2.4f \n', Afs(2), Afd(2));
fprintf('Male Heterogenous solo2 %2.4f duet2 %2.4f \n', Hfs(2), Hfd(2));

 
 %% j160815 Duet 4 (w(7) male and w(8) female) Female solo (1)

% wrengram(w(7).duet, w(7).Fs, 2);
% caxis([-50, 50]);
clear solo duet

% Female
    solo = 1;
      Afs(3) =  getMeanSpikes(w, 8, solo);    
      Hfs(3) =  getMeanSpikes(w, 7, solo);    
    duet = [4 8];
      Afd(3) =  getMeanSpikes(w, 8, duet);    
      Hfd(3) =  getMeanSpikes(w, 7, duet);    
      
 fprintf('Female Autogenous solo3 %2.4f duet3 %2.4f \n', Afs(3), Afd(3));
 fprintf('Male Heterogenous solo3 %2.4f duet3 %2.4f \n', Hfs(3), Hfd(3));
    
    
%% j161009 Duet 5 (w(9) male and w(10) female) Female solo (2)

% wrengram(w(9).duet, w(9).Fs, 2);
% caxis([-50, 50]);
clear solo duet
% Female
    solo{1} = 1;
      Afs(4) =  getMeanSpikes(w, 10, solo{1});    
      Hfs(4) =  getMeanSpikes(w, 9, solo{1});    
    solo{2} = 2;
      Afs(5) =  getMeanSpikes(w, 10, solo{2});    
      Hfs(5) =  getMeanSpikes(w, 9, solo{2});    
    duet{1} = [4 8 12 16];
      Afd(4) =  getMeanSpikes(w, 10, duet{1});    
      Hfd(4) =  getMeanSpikes(w, 9, duet{1});    
    duet{2} = [6 10];
      Afd(5) =  getMeanSpikes(w, 10, duet{2});    
      Hfd(5) =  getMeanSpikes(w, 9, duet{2});    
      
 fprintf('Female Autogenous solo4 %2.4f duet4 %2.4f \n', Afs(4), Afd(4));
 fprintf('Male Heterogenous solo4 %2.4f duet4 %2.4f \n', Hfs(4), Hfd(4));
 fprintf('Female Autogenous solo5 %2.4f duet5 %2.4f \n', Afs(5), Afd(5));
 fprintf('Male Heterogenous solo5 %2.4f duet5 %2.4f \n', Hfs(5), Hfd(5));
    

%% j161022 Duet 6 (w(11) male and w(12) female) Female solo (2)

% wrengram(w(11).duet, w(11).Fs, 2);
% caxis([-50, 50]);
clear solo duet

% Female
    solo{1} = 1;
      Afs(6) =  getMeanSpikes(w, 12, solo{1});    
      Hfs(6) =  getMeanSpikes(w, 11, solo{1});    
    solo{2} = 2;
      Afs(7) =  getMeanSpikes(w, 12, solo{2});    
      Hfs(7) =  getMeanSpikes(w, 11, solo{2});    
    duet{1} = [3 7 11];
      Afd(6) =  getMeanSpikes(w, 12, duet{1});    
      Hfd(6) =  getMeanSpikes(w, 11, duet{1});    
    duet{2} = [5 9 13];
      Afd(7) =  getMeanSpikes(w, 12, duet{2});    
      Hfd(7) =  getMeanSpikes(w, 11, duet{2});    

 fprintf('Female Autogenous solo6 %2.4f duet6 %2.4f \n', Afs(6), Afd(6));
 fprintf('Male Heterogenous solo6 %2.4f duet6 %2.4f \n', Hfs(6), Hfd(6));
 fprintf('Female Autogenous solo7 %2.4f duet7 %2.4f \n', Afs(7), Afd(7));
 fprintf('Male Heterogenous solo7 %2.4f duet7 %2.4f \n', Hfs(7), Hfd(7));
      
      
%% j17060848 Duet 7 (w(13) male and w(14) female) Male solo (2)

% wrengram(w(13).duet, w(13).Fs, 2);
% caxis([-50, 50]);
clear solo duet

    solo{1} = [2 3 4 5 7];
      Ams(4) =  getMeanSpikes(w, 13, solo{1});    
      Hms(4) =  getMeanSpikes(w, 14, solo{1});    
    solo{2} = 6;
      Ams(5) =  getMeanSpikes(w, 13, solo{2});    
      Hms(5) =  getMeanSpikes(w, 14, solo{2});    
    duet{1} = [8 12 16 20];
      Amd(4) =  getMeanSpikes(w, 13, duet{1});    
      Hmd(4) =  getMeanSpikes(w, 14, duet{1});    
    duet{2} = [10 14 18];
      Amd(5) =  getMeanSpikes(w, 13, duet{2});    
      Hmd(5) =  getMeanSpikes(w, 14, duet{2});    
    
 fprintf('Male Autogenous solo4 %2.4f duet4 %2.4f \n', Ams(4), Amd(4));
 fprintf('Female Heterogenous solo4 %2.4f duet4 %2.4f \n', Hms(4), Hmd(4));
 fprintf('Male Autogenous solo5 %2.4f duet5 %2.4f \n', Ams(5), Amd(5));
 fprintf('Female Heterogenous solo5 %2.4f duet5 %2.4f \n', Hms(5), Hmd(5));
    
 
 %%
function out = getMeanSpikes(ww, wrenidx, in)

out = 0;

for j=1:length(in)
    
    for k=1:4
       out = out + length(find(ww(wrenidx).Cspikes{k} > ww(wrenidx).syl(in(j)).tim(1) &  ww(wrenidx).Cspikes{k} < ww(wrenidx).syl(in(j)).tim(2)));
    end
    
end

out = out/length(in);

end

%end
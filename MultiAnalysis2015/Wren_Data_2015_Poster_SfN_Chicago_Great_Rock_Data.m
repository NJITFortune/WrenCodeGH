%use this to make plots showing how the female and male change due to the
%male swtiching the syllables 

%good for chicago and for future papers
%they are more synchronized after the male changes his song

load RockFemL1bDat;
load RockFemL2bDat;
load RockFemR2bDat;
load RockFemR1bDat;
load RockMalL1bDat;
load RockMalL2bDat;
load RockMalR2bDat;
load RockMalR1bDat;
fem1 = RockFemL1b;
fem2 = RockFemL2b;
fem3 = RockFemR1b;
fem4 = RockFemR2b;
mal1 = RockMalL1b;
mal2 = RockMalL2b;
mal3 = RockMalR1b;
mal4 = RockMalR2b;

for j=1:1;
    figure(j); clf;
        ax(1) = subplot(3,1,1); 
            plot(fem1(j).srast.tim, fem1(j).srast.spers, 'b'); 
            hold on;
            plot(fem2(j).srast.tim, fem2(j).srast.spers, 'm');
            plot(fem3(j).srast.tim, fem3(j).srast.spers, 'r');
            plot(fem4(j).srast.tim, fem4(j).srast.spers, 'k');
            hold off;
         ax(2) = subplot(3,1,2); 
            plot(mal1(j).srast.tim, mal1(j).srast.spers, 'b'); 
            hold on;
            plot(mal2(j).srast.tim, mal2(j).srast.spers, 'm');
            plot(mal3(j).srast.tim, mal3(j).srast.spers, 'r');
            plot(mal4(j).srast.tim, mal4(j).srast.spers, 'k');
            hold off;       
        ax(3) = subplot(313); plot(mal1(j).tim, mal1(j).stim); 
        text(-5, 1, mal1(j).StimName);
        linkaxes(ax, 'x');
end;
        
        totalfem = fem1(1).srast.spers + fem2(1).srast.spers + fem3(1).srast.spers + fem4(1).srast.spers;
        totalmal = mal1(1).srast.spers + mal2(1).srast.spers + mal3(1).srast.spers + mal4(1).srast.spers;

 figure(2);
 xa(1) = subplot(211);
        plot(mal4(1).srast.tim, totalfem, 'm');
        hold on;
        plot(mal4(1).srast.tim, totalmal, 'b');
        hold off;

 xa(2) = subplot(212);
        plot(mal4(1).srast.tim, totalfem / max(totalfem), 'm');
        hold on;
        plot(mal4(1).srast.tim, totalmal / max(totalmal), 'b');
        hold off;
linkaxes(xa, 'x');
    
% Cross correlation of whole response

histograminterval = (mal4(1).srast.tim(2) - mal4(1).srast.tim(1));

tt = find(mal4(1).srast.tim > 0 & mal4(1).srast.tim < 34.5);
xmf = xcorr(totalmal(tt) / max(totalmal(tt)), totalfem(tt) / max(totalfem(tt)));

xmftime = histograminterval:histograminterval:histograminterval*length(xmf);
xmftime = xmftime - xmftime(end)/2;

figure(3); plot(xmftime, xmf);

% Seems like the male response is leading the female by about 35 msec (via
% xcorr) which eyeball-o-metrically matches what we see in the
% sliding-window histogram.


% Cross correlation of beginning and end responses

earlywindow = [0.6747 9.2034];
latewindow = [19.8295 28.9136];

tte = find(mal4(1).srast.tim > earlywindow(1) & mal4(1).srast.tim < earlywindow(2));
ttl = find(mal4(1).srast.tim > latewindow(1) & mal4(1).srast.tim < latewindow(2));
earlyXmf = xcorr(totalmal(tte) / max(totalmal(tte)), totalfem(tte) / max(totalfem(tte)));
lateXmf = xcorr(totalmal(ttl) / max(totalmal(ttl)), totalfem(ttl) / max(totalfem(ttl)));
    exmftime = histograminterval:histograminterval:histograminterval*length(earlyXmf);
    exmftime = exmftime - exmftime(end)/2;
    lxmftime = histograminterval:histograminterval:histograminterval*length(lateXmf);
    lxmftime = lxmftime - lxmftime(end)/2;


figure(4); 
    plot(exmftime, earlyXmf, 'r');
    hold on;
    plot(lxmftime, lateXmf, 'k');
    hold off;


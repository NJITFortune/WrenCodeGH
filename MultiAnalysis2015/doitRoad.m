load RoadFemL1aDat;
load RoadFemL2aDat;
load RoadFemL3aDat;
load RoadFemR2aDat;
load RoadFemR1aDat;
load RoadMalL1aDat;
load RoadMalR2aDat;
load RoadMalR1aDat;
fem1 = RoadFemL1a;
fem2 = RoadFemL2a;
fem3 = RoadFemR1a;
fem4 = RoadFemR2a;
fem5 = RoadFemL3a;
mal1 = RoadMalL1a;
mal2 = RoadMalR1a;
mal3 = RoadMalR2a;


%COMMENTS:

%female is leading a little bit

%main thing found: when male drops out, he seems to respond almost double
%the amount that she does (have two instances of him dropping out) -- look
%at the sonogram to see



for j=2:2;
    figure(1); clf;
        ax(1) = subplot(3,1,1); 
            plot(fem1(j).srast.tim, fem1(j).srast.spers, 'b'); 
            hold on;
            plot(fem2(j).srast.tim, fem2(j).srast.spers, 'm');
            plot(fem3(j).srast.tim, fem3(j).srast.spers, 'r');
            plot(fem4(j).srast.tim, fem4(j).srast.spers, 'k');
            plot(fem5(j).srast.tim, fem5(j).srast.spers, 'c');
            hold off;
         ax(2) = subplot(3,1,2); 
            plot(mal1(j).srast.tim, mal1(j).srast.spers, 'b'); 
            hold on;
            plot(mal2(j).srast.tim, mal2(j).srast.spers, 'm');
            plot(mal3(j).srast.tim, mal3(j).srast.spers, 'r');
            hold off;       
        ax(3) = subplot(313); plot(mal1(j).tim, mal1(j).stim); 
        text(-5, 1, mal1(j).StimName);
        linkaxes(ax, 'x');
end;

%because some of them are 15sec and some are 20sec long. OOPS
endlength = min([length(fem1(2).srast.spers) length(fem2(2).srast.spers) length(fem3(2).srast.spers) length(fem4(2).srast.spers) length(fem5(2).srast.spers)]);

totalfem = fem1(2).srast.spers(1:endlength) + fem2(2).srast.spers(1:endlength) + fem3(2).srast.spers(1:endlength) + fem4(2).srast.spers(1:endlength) + fem5(2).srast.spers(1:endlength);
totalmal = mal1(2).srast.spers(1:endlength) + mal2(2).srast.spers(1:endlength) + mal3(2).srast.spers(1:endlength);

 figure(2);
 xa(1) = subplot(211);
        plot(mal3(2).srast.tim, totalfem, 'm');
        hold on;
        plot(mal3(2).srast.tim, totalmal, 'b');
        hold off;

 xa(2) = subplot(212);
        plot(mal3(2).srast.tim, totalfem / max(totalfem), 'm');
        hold on;
        plot(mal3(2).srast.tim, totalmal / max(totalmal), 'b');
        hold off;
linkaxes(xa, 'x');



% Cross correlation of whole response

histograminterval = (mal1(2).srast.tim(2) - mal1(2).srast.tim(1));

tt = find(mal1(2).srast.tim > 0 & mal1(2).srast.tim < 6.45);
xmf = xcorr(totalmal(tt) / max(totalmal(tt)), totalfem(tt) / max(totalfem(tt)));

xmftime = histograminterval:histograminterval:histograminterval*length(xmf);
xmftime = xmftime - xmftime(end)/2;

figure(56); plot(xmftime, xmf);

% Seems like the male response is leading the female by about 35 msec (via
% xcorr) which eyeball-o-metrically matches what we see in the
% sliding-window histogram.


% Cross correlation of beginning and end responses

% earlywindow = [0.6747 9.2034];
% latewindow = [19.8295 28.9136];
% 
% tte = find(mal4(1).srast.tim > earlywindow(1) & mal4(1).srast.tim < earlywindow(2));
% ttl = find(mal4(1).srast.tim > latewindow(1) & mal4(1).srast.tim < latewindow(2));
% earlyXmf = xcorr(totalmal(tte) / max(totalmal(tte)), totalfem(tte) / max(totalfem(tte)));
% lateXmf = xcorr(totalmal(ttl) / max(totalmal(ttl)), totalfem(ttl) / max(totalfem(ttl)));
%     exmftime = histograminterval:histograminterval:histograminterval*length(earlyXmf);
%     exmftime = exmftime - exmftime(end)/2;
%     lxmftime = histograminterval:histograminterval:histograminterval*length(lateXmf);
%     lxmftime = lxmftime - lxmftime(end)/2;
% 
% 
% figure(4); 
%     plot(exmftime, earlyXmf, 'r');
%     hold on;
%     plot(lxmftime, lateXmf, 'k');
%     hold off;
    
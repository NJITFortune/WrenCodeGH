

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
        
    figure(49);
        totalfem = fem1(1).srast.spers + fem2(1).srast.spers + fem3(1).srast.spers + fem4(1).srast.spers;
        totalmal = mal1(1).srast.spers + mal2(1).srast.spers + mal3(1).srast.spers + mal4(1).srast.spers;
        plot(mal4(1).srast.tim, totalfem, 'm');
        hold on;
        plot(mal4(1).srast.tim, totalmal, 'b');
        hold off;
    
end;

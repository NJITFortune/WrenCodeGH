%code for looking at all of the stimuli in Lucky. Summing all the responses

load RoadFemL1bDat;
load RoadFemR2bDat;
load RoadFemR1bDat;

load RoadMalL1bDat;
load RoadMalR2bDat;

fem1 = RoadFemL1b;
fem3 = RoadFemR1b;
fem4 = RoadFemR2b;

mal1 = RoadMalL1b;
mal2 = RoadMalR2b;


for j=1:length(fem1);
    figure(j); clf;
        ax(1) = subplot(3,1,1); 
            plot(fem1(j).srast.tim, fem1(j).srast.spers, 'b'); 
            hold on;
            plot(fem3(j).srast.tim, fem3(j).srast.spers, 'r');
            plot(fem4(j).srast.tim, fem4(j).srast.spers, 'k');
            text(-5, 50, 'Female');
            hold off;
         ax(2) = subplot(3,1,2); 
            plot(mal1(j).srast.tim, mal1(j).srast.spers, 'b'); 
            hold on;
            plot(mal2(j).srast.tim, mal2(j).srast.spers, 'm');
            text(-5, 50, 'Male');
            hold off; 
     
        ax(3) = subplot(313); plot(mal1(j).tim, mal1(j).stim); 
            text(-5, 1, mal1(j).StimName);
            
        linkaxes(ax, 'x');
end;

endlength = min([length(fem1(2).srast.spers) length(fem3(2).srast.spers) length(fem4(2).srast.spers)]);

for k = 1:length(fem1);
    totalfem{k} = fem1(k).srast.spers(1:endlength) + fem3(k).srast.spers(1:endlength) + fem4(k).srast.spers(1:endlength);
    totalmal{k} = mal1(k).srast.spers(1:endlength) + mal2(k).srast.spers(1:endlength);
end

uu = length(fem1);

    for j=1:length(fem1);
        figure(j+uu); clf;
        xa(1) = subplot(211);
            plot(mal1(j).srast.tim(1:endlength), totalfem{j}, 'm');
            hold on;
            plot(mal1(j).srast.tim(1:endlength), totalmal{j}, 'b');
            text(15, 100, mal1(j).StimName);
            hold off;

        xa(2) = subplot(212);
            plot(mal1(j).srast.tim(1:endlength), totalfem{j} / max(totalfem{j}), 'm');
            hold on;
            plot(mal1(j).srast.tim(1:endlength), totalmal{j} / max(totalmal{j}), 'b');
            text(15, .5, mal1(j).StimName);
            hold off;
        linkaxes(xa, 'x');
    end


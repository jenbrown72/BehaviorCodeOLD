
function [] = JB_plotMovingPerformance(basicPropertiesToPlot)

percentCorrectChance = [0.5 0.5];

for j = 1:length(basicPropertiesToPlot);
    if strcmp('S1auto', basicPropertiesToPlot{j,1}.sessionType)
        sessionCode(j,1) = 1;
    elseif strcmp('S1',basicPropertiesToPlot{j,1}.sessionType)
        sessionCode(j,1) = 2;
    elseif strcmp('S2',basicPropertiesToPlot{j,1}.sessionType)
        sessionCode(j,1) = 3;
    elseif strcmp('S6',basicPropertiesToPlot{j,1}.sessionType)
        sessionCode(j,1) = 4;
    elseif strcmp('S8',basicPropertiesToPlot{j,1}.sessionType)
        sessionCode(j,1) = 5;
    elseif strcmp('S10',basicPropertiesToPlot{j,1}.sessionType)
        sessionCode(j,1) = 6;
    elseif strcmp('S12',basicPropertiesToPlot{j,1}.sessionType)
        sessionCode(j,1) = 7;
    end
    
    blockPerformance = basicPropertiesToPlot{j,1}.temptrialTypes;
    
    blockPerformance(blockPerformance==4)=1; %CR
    blockPerformance(blockPerformance==1)=1; %Hit
    blockPerformance(blockPerformance==2)=0; %FA
    blockPerformance(blockPerformance==3)=0; %Miss
    
    tally = 0;
    basicPropertiesToPlot{j,1}.rolloingperformance = [];
    for kk = 1:length(blockPerformance)
        basicPropertiesToPlot{j,1}.rolloingperformance(kk,1) = (blockPerformance(kk)+tally)/kk;
        tally = tally+blockPerformance(kk);
    end
end

ColorCodeColors = [1;0.9;0.8;0.6;0.4;0.2;0];
spaceBetweenSessions = 40;
figure(88);clf;
cummX = 1;
for jj = 1:length(basicPropertiesToPlot)
    if ~isempty(basicPropertiesToPlot{jj,1}.rolloingperformance) && (length(basicPropertiesToPlot{jj,1}.rolloingperformance)>100);
        movingAvg=movingmean( basicPropertiesToPlot{jj,1}.rolloingperformance,61);
        plot(cummX:cummX+length(movingAvg)-1,movingAvg,'-', 'Color',[0.4,ColorCodeColors(sessionCode(jj)),0.8],'LineWidth',4)
        hold on
        cummX = cummX + length(movingAvg) + spaceBetweenSessions;
    end
end
plot([min(xlim) max(xlim)],percentCorrectChance,'k--','LineWidth',2);
ylim([0 1])
ylabel('Moving fraction Correct')
xlabel('Trial No.')
end






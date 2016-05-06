function [] = JB_plotPerformanceFig(basicPropertiesToPlot,plotON)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

positionGraph1 = [14   83   503   869];
dprimeThreshold = [1 1];
percentCorrectChance = [0.5 0.5];
markerSizePlots=8;

if (plotON==1)
    f = figure;clf
    set(f,'Position',positionGraph1);
else
    figure('Visible','off');clf;
end

tempTitle = basicPropertiesToPlot{1,1}.mouseID;
tempTitle(findstr(tempTitle,'_'))=[];
set(gcf,'name',tempTitle,'numbertitle','off');

noSubPlots = 2;
subplot(noSubPlots,1,2)
numPoints = 1:1:length(basicPropertiesToPlot);

SessionTypes = {'S1auto' ; 'S1'; 'S2'; 'S6'; 'S8'; 'S10'; 'S12'};

for j = 1:length(basicPropertiesToPlot);
    if strcmp('S1auto', basicPropertiesToPlot{j,1}.sessionType)
        colorCode(j,1) = 1;
    elseif strcmp('S1',basicPropertiesToPlot{j,1}.sessionType)
        colorCode(j,1) = 2;
    elseif strcmp('S2',basicPropertiesToPlot{j,1}.sessionType)
        colorCode(j,1) = 3;
    elseif strcmp('S6',basicPropertiesToPlot{j,1}.sessionType)
        colorCode(j,1) = 4;
    elseif strcmp('S8',basicPropertiesToPlot{j,1}.sessionType)
        colorCode(j,1) = 5;
    elseif strcmp('S10',basicPropertiesToPlot{j,1}.sessionType)
        colorCode(j,1) = 6;
    elseif strcmp('S12',basicPropertiesToPlot{j,1}.sessionType)
        colorCode(j,1) = 7;
    end
end
ColorCodeColors = [0.9;0.8;0.6;0.4;0.2;0];

for jj = 1:length(basicPropertiesToPlot);
    dataTemp(jj,:)  = [basicPropertiesToPlot{jj,1}.dprime colorCode(jj,1)];
end

sessionNoStop = (find((dataTemp(:,2)==5),1))+1;

%plot performance
diffData = diff(dataTemp(:,2));
subplot(noSubPlots,1,1)
for jj = 1:sessionNoStop-1;
    if diffData(jj)==1;
        plot(numPoints(jj),basicPropertiesToPlot{jj,1}.sessionperformance,'o','MarkerSize', markerSizePlots, 'linewidth',2,'MarkerEdgeColor',[0.4,ColorCodeColors(colorCode(jj)),0.8], 'MarkerFaceColor',[0.4,ColorCodeColors(colorCode(jj)),0.8],'Color', [0.4,ColorCodeColors(colorCode(jj)),0.8]);
        hold on
    else
        plot([numPoints(jj) numPoints(jj+1)],[basicPropertiesToPlot{jj,1}.sessionperformance basicPropertiesToPlot{jj+1,1}.sessionperformance],'o-','MarkerSize', markerSizePlots, 'linewidth',2,'MarkerEdgeColor',[0.4,ColorCodeColors(colorCode(jj)),0.8], 'MarkerFaceColor',[0.4,ColorCodeColors(colorCode(jj)),0.8],'Color', [0.4,ColorCodeColors(colorCode(jj)),0.8]);
        hold on
    end
end
plot([0 24],percentCorrectChance,'k--')
ylim([0 1])
xlim([0 24])
ylabel('Performance');
xlabel('Session Number');

%plot Dprime
subplot(noSubPlots,1,2)
for jj = 1:sessionNoStop-1;
    if diffData(jj)==1;
        plot(numPoints(jj),basicPropertiesToPlot{jj,1}.dprime,'o','MarkerSize', markerSizePlots, 'linewidth',2,'MarkerEdgeColor',[0.4,ColorCodeColors(colorCode(jj)),0.8], 'MarkerFaceColor',[0.4,ColorCodeColors(colorCode(jj)),0.8],'Color', [0.4,ColorCodeColors(colorCode(jj)),0.8]);
        hold on
    else
        plot([numPoints(jj) numPoints(jj+1)],[basicPropertiesToPlot{jj,1}.dprime basicPropertiesToPlot{jj+1,1}.dprime],'o-','MarkerSize', 10, 'linewidth',2,'MarkerEdgeColor',[0.4,ColorCodeColors(colorCode(jj)),0.8], 'MarkerFaceColor',[0.4,ColorCodeColors(colorCode(jj)),0.8],'Color', [0.4,ColorCodeColors(colorCode(jj)),0.8]);
        hold on
    end
end
hold on
plot([0 24],dprimeThreshold,'k--')
xlim([0 24])
ylabel('dprime');
xlabel('Session Number');

end




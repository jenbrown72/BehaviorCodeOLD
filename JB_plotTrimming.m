function [DATAtrim] = JB_plotTrimming(basicPropertiesToPlot, plotON, None)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
%basicPropertiesToPlot
%plotON = 1, plot fig, 0 = no figure
%None = 1, include no whiskers in plots, 0 = exclude no whiskers

positionGraph1 = [897   382   422   614];
positionGraph2 = [1321 383 404 613];
positionGraph3 = [7   385   888   611];

percentCorrectChance = [0.5 0.5];
dprimeThreshold = [1 1];
trimtally=1;
numPoints = 1:1:length(basicPropertiesToPlot);
indTrim = nan(length(numPoints),1);


for h = 1:length(numPoints)
    if ~isempty(basicPropertiesToPlot{h,1}.trimming)
        indTrim(h,1)=1;
    else
        indTrim(h,1)=0;
    end
end

%take the two last sessions before trimming for baseline
firstIndTrim = find(indTrim,1);

DATAtrim.tempTrimPerformance = nan(length(numPoints),1);
DATAtrim.tempTrimdPrime = nan(length(numPoints),1);
DATAtrim.tempTrimWhiskers = cell(length(numPoints),1);
clear DATAtrim.tempTrimAnglesdPrime DATAtrim.tempTrimAnglesDiff
if ~isempty(firstIndTrim)
    disp(' ');
    if (strcmp(basicPropertiesToPlot{firstIndTrim,1}.trimType,'Row'))
        disp('Analysing a Row data Set');
        DATAtrim.trimType = 'Row';
    else
        disp('Analysing a Column data Set');
        DATAtrim.trimType = 'Column';
    end
    disp(' ');
    
    for h=1:length(numPoints)
        if h==firstIndTrim-1 || h==firstIndTrim-2
            DATAtrim.tempTrimPerformance(trimtally,1) = basicPropertiesToPlot{h,1}.sessionperformance;
            DATAtrim.tempTrimdPrime(trimtally,1) = basicPropertiesToPlot{h,1}.dprime;
            DATAtrim.tempTrimWhiskers{trimtally,1} = 'Full';
            DATAtrim.tempTrimAnglesdPrime(trimtally,:) = basicPropertiesToPlot{h,1}.pairsDprime;
            DATAtrim.tempTrimAnglesDiff(trimtally,:) = basicPropertiesToPlot{h,1}.pairsDiff;
            trimtally = trimtally+1;
            
        elseif indTrim(h,1)==1
            DATAtrim.tempTrimPerformance(trimtally,1) = basicPropertiesToPlot{h,1}.sessionperformance;
            DATAtrim.tempTrimdPrime(trimtally,1) = basicPropertiesToPlot{h,1}.dprime;
            DATAtrim.tempTrimWhiskers{trimtally,1} = basicPropertiesToPlot{h,1}.trimming;
            DATAtrim.tempTrimAnglesdPrime(trimtally,:) = basicPropertiesToPlot{h,1}.pairsDprime;
            DATAtrim.tempTrimAnglesDiff(trimtally,:) = basicPropertiesToPlot{h,1}.pairsDiff;
            trimtally = trimtally+1;
        end
    end
end


%delete empty cells
DATAtrim.tempTrimWhiskers(any(cellfun(@isempty,DATAtrim.tempTrimWhiskers),2),:)=[];
cell_A = strrep(DATAtrim.tempTrimWhiskers,' ','');

list=1;
DATAtrim.PairType = cell(length(numPoints),1);
DATAtrim.Performance = nan(length(numPoints),2);
DATAtrim.dPrime = nan(length(numPoints),2);

%Find pairs of stim types
for k = 1:length(cell_A)-1;
    if (None==0)
        if(strcmp(cell_A{k},'None'));
            continue
        end
    end
    
    if (strcmp(cell_A{k},cell_A{k+1}))
        DATAtrim.PairType{list} = cell_A{k};
        DATAtrim.Performance(list,:) = [DATAtrim.tempTrimPerformance(k), DATAtrim.tempTrimPerformance(k+1)];
        DATAtrim.dPrime(list,:) = [DATAtrim.tempTrimdPrime(k), DATAtrim.tempTrimdPrime(k+1)];
        list = list+1;
    end
end

DATAtrim.PairType(any(cellfun(@isempty,DATAtrim.PairType),2),:)=[];
DATAtrim.Performance = DATAtrim.Performance(~any(isnan(DATAtrim.Performance),2),:);
DATAtrim.dPrime = DATAtrim.dPrime(~any(isnan(DATAtrim.dPrime),2),:);

% sigPairs = cell;
clear sigPairs sigPairsPval h p tessting
testing = DATAtrim.Performance';
k=1;
%h = nan(length(testing)

for kk = 1:length(testing)-1
    [h(kk),p(kk)] = ttest(testing(:,1),testing(:,kk+1));
    if h(kk)==1
        sigPairs{k} = [1,kk];
        sigPairsPval(k) = [p(kk)];
        k = k+1;
    end
end


plotRows = 2;
plotCols = 2;

if (plotON==1)
    ffff=figure;clf
    set(ffff,'Position',positionGraph2);
else
    figure('Visible','off');clf;
end

tempTitle = basicPropertiesToPlot{1,1}.mouseID;
tempTitle(findstr(tempTitle,'_'))=[];
tempTitle = strcat(tempTitle,'_PerformanceTrimming');
set(gcf,'name',tempTitle,'numbertitle','off')
currPlot = 1;

%plot avergae of each pair.
subplot(plotRows,plotCols,currPlot);
bar(mean(DATAtrim.Performance,2), 0.5);
hold on
errorbar(mean(DATAtrim.Performance,2),std(DATAtrim.Performance,0,2),'.');
%sigstar(sigPairs,sigPairsPval);
set(gca,'XTickLabel',DATAtrim.PairType);
ylabel('average Performance')
currPlot=currPlot+1;

%plot just the first occurance of the pair
subplot(plotRows,plotCols,currPlot);
bar(DATAtrim.Performance(:,1), 0.5);
hold on
set(gca,'XTickLabel',DATAtrim.PairType);
ylabel('First Session Performance')
currPlot=currPlot+1;


testing = DATAtrim.dPrime';
k=1;
clear h p sigPairs sigPairsPval

for kk = 1+1:length(testing)
    [h(kk),p(kk)] = ttest(testing(:,1),testing(:,kk));
    if h(kk)==1
        sigPairs{k} = [1,kk];
        sigPairsPval(k) = [p(kk)];
        k = k+1;
    end
end

subplot(plotRows,plotCols,currPlot);
bar(mean(DATAtrim.dPrime,2), 0.5);
hold on
errorbar(mean(DATAtrim.dPrime,2),std(DATAtrim.dPrime,0,2),'.');
%sigstar(sigPairs,sigPairsPval);
set(gca,'XTickLabel',DATAtrim.PairType);
ylabel('average dPrime')
currPlot=currPlot+1;

subplot(plotRows,plotCols,currPlot);
bar(DATAtrim.Performance(:,1), 0.5);
set(gca,'XTickLabel',DATAtrim.PairType);
ylabel('First Session dPrime')
hold on

%plot performance and d prime over session

plotRows = 2;
plotCols = 1;
plotTally = (plotRows*plotCols);
numFigs = 1;

if (plotON==1)
    ffff=figure;clf
    set(ffff,'Position',positionGraph1);
else
    figure('Visible','off');clf;
end

tempTitle = basicPropertiesToPlot{1,1}.mouseID;
tempTitle(findstr(tempTitle,'_'))=[];
set(gcf,'name',tempTitle,'numbertitle','off')
currPlot = 1;

subplot(plotRows,plotCols,currPlot);
plot(DATAtrim.tempTrimPerformance,'o-k','LineWidth',5);
hold on
set(gca,'XTick',[1:length(DATAtrim.tempTrimPerformance)]);
set(gca,'XTickLabel',DATAtrim.tempTrimWhiskers);
ylim([0 1]);
set(gca, 'Ylim',[0 1])
ylabel('Performance');
xlabel('whisker trimming');
currPlot=currPlot+1;
plot([1 max(xlim)],percentCorrectChance,'k--','LineWidth',2)

subplot(plotRows,plotCols,currPlot);
plot(DATAtrim.tempTrimdPrime,'o-k','LineWidth',5);
hold on
set(gca,'XTick',[1:length(DATAtrim.tempTrimdPrime)]);
set(gca,'XTickLabel',DATAtrim.tempTrimWhiskers);
plot([1 max(xlim)],dprimeThreshold,'k--','LineWidth',2)
ylimit = ylim;

if ylimit(1)<0;
    yLimitSet = ylimit(1);
else
    yLimitSet = 0;
end
ylim([yLimitSet ylimit(2)]);
set(gca, 'Ylim',[yLimitSet max(DATAtrim.tempTrimdPrime)])
ylabel('d prime');
xlabel('whisker trimming');
currPlot=currPlot+1;

saveas(gca,fullfile('C:\Users\adesniklab\Documents\BehaviorRawData\currFigs\trimming',tempTitle),'jpeg');


%Plot dprime over angles
%delete empty cells

plotRows = 3;
plotCols = size((DATAtrim.tempTrimAnglesdPrime),2);
plotTally = (plotRows*plotCols);
numFigs = 1;

if (plotON==1)
    ffff=figure;clf
    set(ffff,'Position',positionGraph3);
else
    figure('Visible','off');clf;
end

tempTitle = basicPropertiesToPlot{1,1}.mouseID;
tempTitle(findstr(tempTitle,'_'))=[];
tempTitle = strcat(tempTitle,'_Angles');
set(gcf,'name',tempTitle,'numbertitle','off')
currPlot = 1;

subplot(plotRows,plotCols,currPlot:(plotCols*2));

currPlot = (plotCols*2)+1;
for jj = 1:size((DATAtrim.tempTrimAnglesdPrime),2);
    plot(DATAtrim.tempTrimAnglesdPrime(:,jj),'o-','LineWidth',5,'Color',[0,1/(jj),1/(jj)]);
    templegend{jj} = ['diffAngle ',num2str(DATAtrim.tempTrimAnglesDiff(1,jj))];
    hold on
end
legend(templegend)

set(gca,'XTick',(1:length(DATAtrim.tempTrimAnglesdPrime)));
set(gca,'XTickLabel',DATAtrim.tempTrimWhiskers);

maxY = max(DATAtrim.tempTrimAnglesdPrime(:));
minY = min(DATAtrim.tempTrimAnglesdPrime(:));
for jj = 1:size((DATAtrim.tempTrimAnglesdPrime),2);
    subplot(plotRows,plotCols,currPlot);
    plot(DATAtrim.tempTrimAnglesdPrime(:,jj),'o-','LineWidth',5,'Color',[0,1/(jj),1/(jj)]);
    hold on
    set(gca,'XTick',[1:length(DATAtrim.tempTrimAnglesdPrime)]);
    set(gca,'XTickLabel',DATAtrim.tempTrimWhiskers);
    currPlot = currPlot+1;
    set(gca,'Ylim',[minY maxY]);
    xlimit = xlim;
    ylimit = ylim;
    str = strcat('diffAngle ',templegend(jj));
    text(xlimit(1),(maxY+0.5),str)
    
end
saveas(gca,fullfile('C:\Users\adesniklab\Documents\BehaviorRawData\currFigs\trimming',tempTitle),'jpeg');
end


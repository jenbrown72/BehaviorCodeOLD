function [DATAtrim] = JB_plotTrimming(basicPropertiesToPlot, plotON, None)

 % JB_plotTrimming  plots individual trimming data 

 %   [DATAtrim] = JB_plotTrimming(basicPropertiesToPlot) - data
 %   inputed - returned from JB_basicBehaviorProperties
 %   plotON = 1 (plots output), = 0 (no plot generated) 
 %   None = 1 (included non whisker days in plot, = 0, excluded no
 %   whiskers

 % Examples:
 %   [DATAtrim] = JB_plotTrimming(basicPropertiesToPlot, 1, 0);
 %
 %   [DATAtrim] = JB_plotTrimming(basicPropertiesToPlot, 1, 1);
 %   [DATATRIM] = returns a matrix of average data to be used in group analysis
 %   


positionGraph1 = [897   382   422   614];
positionGraph2 = [1097         383         628         613];
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
%firstIndTrim = find(diff(indTrim)>0);
firstIndTrim = find(indTrim,1);

DATAtrim.performanceAll = nan(length(numPoints),1);
DATAtrim.dPrimeAll = nan(length(numPoints),1);
DATAtrim.whiskerIDAll = cell(length(numPoints),1);
clear DATAtrim.pairedAnglesdPrime DATAtrim.pairedAnglesDiff

if ~isempty(firstIndTrim)
    for kk = 1:length(firstIndTrim);
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
            DATAtrim.performanceAll(trimtally,1) = basicPropertiesToPlot{h,1}.sessionperformance;
            DATAtrim.dPrimeAll(trimtally,1) = basicPropertiesToPlot{h,1}.dprime;
            DATAtrim.whiskerIDAll{trimtally,1} = 'Full';
            DATAtrim.pairedAnglesdPrime(trimtally,:) = basicPropertiesToPlot{h,1}.pairsDprime;
            DATAtrim.pairedAnglesDiff(trimtally,:) = basicPropertiesToPlot{h,1}.pairsDiff;
            trimtally = trimtally+1;
            
        elseif indTrim(h,1)==1
            DATAtrim.performanceAll(trimtally,1) = basicPropertiesToPlot{h,1}.sessionperformance;
            DATAtrim.dPrimeAll(trimtally,1) = basicPropertiesToPlot{h,1}.dprime;
            DATAtrim.whiskerIDAll{trimtally,1} = basicPropertiesToPlot{h,1}.trimming;
            DATAtrim.pairedAnglesdPrime(trimtally,:) = basicPropertiesToPlot{h,1}.pairsDprime;
            DATAtrim.pairedAnglesDiff(trimtally,:) = basicPropertiesToPlot{h,1}.pairsDiff;
            trimtally = trimtally+1;
        end
    end
end


%delete empty cells
DATAtrim.performanceAll(isnan(DATAtrim.performanceAll))=[];
DATAtrim.dPrimeAll(isnan(DATAtrim.dPrimeAll))=[];
DATAtrim.whiskerIDAll(any(cellfun(@isempty,DATAtrim.whiskerIDAll),2),:)=[];

cell_A = strrep(DATAtrim.whiskerIDAll,' ','');
list=1;
DATAtrim.whiskerID = cell(length(numPoints),1);
DATAtrim.performance = nan(length(numPoints),5);
DATAtrim.dPrime = nan(length(numPoints),5);

for jj = 1:size((DATAtrim.pairedAnglesdPrime),2);
    tempName = ['anglesdPrime',num2str(DATAtrim.pairedAnglesDiff(1,jj))];
    angleName{jj} = tempName;
    DATAtrim.(tempName) = nan(length(numPoints),5);
end

%Find repeats of stim types

uniqueSessionID = unique(cell_A,'stable');

for kk = 1:length(uniqueSessionID);
    
    [matchID] = strcmp(cell_A,uniqueSessionID(kk));
    if any(matchID) %if there was a hit
        DATAtrim.whiskerID{list} = uniqueSessionID{kk};
        idx = find(matchID); %find the indics of the matches
        DATAtrim.performance(list,1:length(idx)) = DATAtrim.performanceAll(idx)';
        DATAtrim.dPrime(list,1:length(idx)) = DATAtrim.dPrimeAll(idx)';
    end
    
    for jj = 1:size((DATAtrim.pairedAnglesdPrime),2);
        DATAtrim.(angleName{jj})(list,1:length(idx)) = DATAtrim.pairedAnglesdPrime(idx,jj)';
    end
    list = list+1;
end

DATAtrim.performance(~any(~isnan(DATAtrim.performance),2),:) = [];%remove rows with nan
DATAtrim.dPrime(~any(~isnan(DATAtrim.dPrime),2),:) = []; %remove rows with nan

        if (None==1)
            idx = find(strcmp(DATAtrim.whiskerID,'None'));
            DATAtrim.performance(idx,:) = [];
            DATAtrim.dPrime(idx,:) = [];
        end
        

% DATAtrim.performance(:,any(isnan(DATAtrim.performance),1))=[]; %remove columns with nan
% DATAtrim.dPrime(:,any(isnan(DATAtrim.dPrime),1))=[]; %remove columns with nan
DATAtrim.whiskerID(any(cellfun(@isempty,DATAtrim.whiskerID),2),:)=[];

for jj = 1:size((DATAtrim.pairedAnglesdPrime),2);
    DATAtrim.(angleName{jj})(~any(~isnan(DATAtrim.(angleName{jj})),2),:) = [];
    %     DATAtrim.(angleName{jj})(:,any(isnan(DATAtrim.(angleName{jj})),1))=[]; %remove columns with nan
end

% sigPairs = cell;
clear sigPairs sigPairsPval h p tessting
testing = DATAtrim.performance';
k=1;
%h = nan(length(testing)

for kk = 1:size((testing),2)-1
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

%plot just the first occurance of the pair
subplot(plotRows,plotCols,currPlot);
bar(DATAtrim.performance(:,1), 0.5);
hold on
set(gca,'XTickLabel',DATAtrim.whiskerID);
ylabel('First Session Performance')
currPlot=currPlot+1;


%plot avergae of each pair.
%calculate the min number of observations per data point
minNoDATApoints = min(sum(~isnan(DATAtrim.dPrime),2));

subplot(plotRows,plotCols,currPlot);
bar(nanmean(DATAtrim.performance(:,1:minNoDATApoints),2), 0.5);
hold on
errorbar(nanmean(DATAtrim.performance,2),nanstd(DATAtrim.performance(:,1:minNoDATApoints),0,2),'.');
%sigstar(sigPairs,sigPairsPval);
set(gca,'XTickLabel',DATAtrim.whiskerID);
ylabel('average Performance')
currPlot=currPlot+1;

testing = DATAtrim.dPrime';
k=1;
clear h p sigPairs sigPairsPval

for kk = 1:size((testing),2)-1
    [h(kk),p(kk)] = ttest(testing(:,1),testing(:,kk+1));
    if h(kk)==1
        sigPairs{k} = [1,kk+1];
        sigPairsPval(k) = [p(kk)];
        k = k+1;
    end
end


subplot(plotRows,plotCols,currPlot);
bar(DATAtrim.performance(:,1), 0.5);
set(gca,'XTickLabel',DATAtrim.whiskerID);
ylabel('First Session dPrime')
hold on
currPlot=currPlot+1;


subplot(plotRows,plotCols,currPlot);
bar(nanmean(DATAtrim.dPrime(:,1:minNoDATApoints),2), 0.5);
hold on
errorbar(nanmean(DATAtrim.dPrime(:,1:minNoDATApoints),2),nanstd(DATAtrim.dPrime(:,1:minNoDATApoints),0,2),'.');
%sigstar(sigPairs,sigPairsPval);
set(gca,'XTickLabel',DATAtrim.whiskerID);
ylabel('average dPrime')
currPlot=currPlot+1;


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
plot(DATAtrim.performanceAll,'o-k','LineWidth',5);
hold on
set(gca,'XTick',[1:length(DATAtrim.performanceAll)]);
set(gca,'XTickLabel',DATAtrim.whiskerIDAll);
ylim([0 1]);
set(gca, 'Ylim',[0 1])
ylabel('Performance');
xlabel('whisker trimming');
currPlot=currPlot+1;
plot([1 max(xlim)],percentCorrectChance,'k--','LineWidth',2)

subplot(plotRows,plotCols,currPlot);
plot(DATAtrim.dPrimeAll,'o-k','LineWidth',5);
hold on
set(gca,'XTick',[1:length(DATAtrim.dPrimeAll)]);
set(gca,'XTickLabel',DATAtrim.whiskerIDAll);
plot([1 max(xlim)],dprimeThreshold,'k--','LineWidth',2)
ylimit = ylim;

if ylimit(1)<0;
    yLimitSet = ylimit(1);
else
    yLimitSet = 0;
end
ylim([yLimitSet ylimit(2)]);
set(gca, 'Ylim',[yLimitSet max(DATAtrim.dPrimeAll)])
ylabel('d prime');
xlabel('whisker trimming');
currPlot=currPlot+1;

saveas(gca,fullfile('C:\Users\adesniklab\Documents\BehaviorRawData\currFigs\trimming',tempTitle),'jpeg');
%Plot dprime over angles
%delete empty cells

plotRows = 3;
plotCols = size((DATAtrim.pairedAnglesdPrime),2);
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
for jj = 1:size((DATAtrim.pairedAnglesdPrime),2);
    plot(DATAtrim.pairedAnglesdPrime(:,jj),'o-','LineWidth',5,'Color',[0,1/(jj),1/(jj)]);
    templegend{jj} = ['diffAngle ',num2str(DATAtrim.pairedAnglesDiff(1,jj))];
    hold on
end
legend(templegend)

set(gca,'XTick',(1:length(DATAtrim.pairedAnglesdPrime)));
set(gca,'XTickLabel',DATAtrim.whiskerIDAll);

maxY = max(DATAtrim.pairedAnglesdPrime(:));
minY = min(DATAtrim.pairedAnglesdPrime(:));
for jj = 1:size((DATAtrim.pairedAnglesdPrime),2);
    subplot(plotRows,plotCols,currPlot);
    plot(DATAtrim.pairedAnglesdPrime(:,jj),'o-','LineWidth',5,'Color',[0,1/(jj),1/(jj)]);
    hold on
    set(gca,'XTick',[1:length(DATAtrim.pairedAnglesdPrime)]);
    set(gca,'XTickLabel',DATAtrim.whiskerIDAll);
    currPlot = currPlot+1;
    set(gca,'Ylim',[minY maxY]);
    xlimit = xlim;
    ylimit = ylim;
    str = strcat('diffAngle ',templegend(jj));
    text(xlimit(1),(maxY+0.5),str)
    
end
saveas(gca,fullfile('C:\Users\adesniklab\Documents\BehaviorRawData\currFigs\trimming',tempTitle),'jpeg');
end


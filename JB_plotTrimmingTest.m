function [DATAtrim] = JB_plotTrimmingTest(basicPropertiesToPlot, plotON, None)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
%basicPropertiesToPlot
%plotON = 1, plot fig, 0 = no figure
%None = 1, include no whiskers in plots, 0 = exclude no whiskers

clear DATAtrim;
positionGraph1 = [897   382   422   614];
positionGraph2 = [1097         383         628         613];
positionGraph3 = [7   385   888   611];

percentCorrectChance = [0.5 0.5];
dprimeThreshold = [1 1];
numPoints = 1:1:length(basicPropertiesToPlot);
indTrim = nan(length(numPoints),1);

for h = 1:length(numPoints)
    if ~isempty(basicPropertiesToPlot{h,1}.trimming)
        indTrim(h,1)=1;
    else
        indTrim(h,1)=0;
    end
end

%identify first session of trimming
firstIndTrim = find(diff(indTrim)>0)+1;

for ii = 1:length(firstIndTrim)
    DATAtrim.performanceAll{ii} = NaN(length(numPoints),1);
    DATAtrim.dPrimeAll{ii} = NaN(length(numPoints),1);
    DATAtrim.performance{ii} = NaN(length(numPoints),5);
    DATAtrim.dPrime{ii} = NaN(length(numPoints),5);
end

% DATAtrim.performanceAll = cell((length(firstIndTrim)),1);
% DATAtrim.dPrimeAll = cell((length(firstIndTrim)),1);
% DATAtrim.whiskerIDAll = cell((length(firstIndTrim)),1);
% DATAtrim.whiskerID = cell((length(firstIndTrim)),1);
% DATAtrim.performance = cell((length(firstIndTrim)),1);
% DATAtrim.dPrime = cell((length(firstIndTrim)),1);

if ~isempty(firstIndTrim)
    for kk = 1:length(firstIndTrim);
        % for kk = 2;
        trimtally=1;
        disp(' ');
        if (strcmp(basicPropertiesToPlot{firstIndTrim(kk),1}.trimType,'Row'))
            disp('Analysing a Row data Set');
            DATAtrim.trimType = 'Row';
        else
            disp('Analysing a Column data Set');
            DATAtrim.trimType = 'Column';
        end
        disp(' ');
        
        for h=1:length(numPoints)
            if h==firstIndTrim(kk)-1 || h==firstIndTrim(kk)-2
                DATAtrim.performanceAll{kk}(trimtally,1) = basicPropertiesToPlot{h,1}.sessionperformance;
                DATAtrim.dPrimeAll{kk}(trimtally,1) = basicPropertiesToPlot{h,1}.dprime;
                DATAtrim.whiskerIDAll{kk}{trimtally,1} = 'Full';
                DATAtrim.pairedAnglesdPrime{kk}(trimtally,:) = basicPropertiesToPlot{h,1}.pairsDprime;
                DATAtrim.pairedAnglesDiff{kk}(trimtally,:) = basicPropertiesToPlot{h,1}.pairsDiff;
                trimtally = trimtally+1;
                
            elseif indTrim(h,1)==1 && h>=firstIndTrim(kk);
                
                if kk==1 && length(firstIndTrim)>1 && h>=firstIndTrim(kk+1) ;
                    continue;
                else
                    DATAtrim.performanceAll{kk}(trimtally,1) = basicPropertiesToPlot{h,1}.sessionperformance;
                    DATAtrim.dPrimeAll{kk}(trimtally,1) = basicPropertiesToPlot{h,1}.dprime;
                    DATAtrim.whiskerIDAll{kk}{trimtally,1} = basicPropertiesToPlot{h,1}.trimming;
                    DATAtrim.pairedAnglesdPrime{kk}(trimtally,:) = basicPropertiesToPlot{h,1}.pairsDprime;
                    DATAtrim.pairedAnglesDiff{kk}(trimtally,:) = basicPropertiesToPlot{h,1}.pairsDiff;
                    trimtally = trimtally+1;
                end
            end
        end
        
        for jj = 1:size((DATAtrim.pairedAnglesdPrime{kk}),2);
            tempName = ['anglesdPrime',num2str(DATAtrim.pairedAnglesDiff{kk}(1,jj))];
            angleName{jj} = tempName;
            DATAtrim.(tempName){kk} = nan(length(numPoints),5);
        end
        
        cell_A = strrep(DATAtrim.whiskerIDAll{kk},' ','');
        list=1;
        
        %Find repeats of stim types
        uniqueSessionID = unique(cell_A,'stable');
        for kkk = 1:length(uniqueSessionID);
            [matchID] = strcmp(cell_A,uniqueSessionID(kkk));
            if any(matchID) %if there was a hit
                DATAtrim.whiskerID{list,kk} = uniqueSessionID{kkk};
                idx = find(matchID); %find the indics of the matches
                DATAtrim.performance{kk}(list,1:length(idx)) = DATAtrim.performanceAll{kk}(idx)';
                DATAtrim.dPrime{kk}(list,1:length(idx)) = DATAtrim.dPrimeAll{kk}(idx)';
            end
            
            for jj = 1:size((DATAtrim.pairedAnglesdPrime{kk}),2);
                DATAtrim.(angleName{jj}){kk}(list,1:length(idx)) = DATAtrim.pairedAnglesdPrime{kk}(idx,jj)';
            end
            list = list+1;
        end
        
        DATAtrim.performance{kk}(~any(~isnan(DATAtrim.performance{kk}),2),:) = [];%remove rows with nan
        DATAtrim.dPrime{kk}(~any(~isnan(DATAtrim.dPrime{kk}),2),:) = []; %remove rows with nan
        
        DATAtrim.performance{kk}(:,any(isnan( DATAtrim.performance{kk}),1))=[];%remove cols with nan
        DATAtrim.dPrime{kk}(:,any(isnan( DATAtrim.dPrime{kk}),1))=[];%remove cols with nan
        
        if (None==1)
            idx = find(strcmp(DATAtrim.whiskerID(:,kk),'None'));
            DATAtrim.performance{kk}(idx,:) = [];
            DATAtrim.dPrime{kk}(idx,:) = [];
        end
        
        for jj = 1:size((DATAtrim.pairedAnglesdPrime{kk}),2);
            DATAtrim.(angleName{jj}){kk}(~any(~isnan(DATAtrim.(angleName{jj}){kk}),2),:) = [];
            DATAtrim.(angleName{jj}){kk}(:,any(isnan( DATAtrim.(angleName{jj}){kk}),1))=[];
        end
        
        clear sigPairs sigPairsPval h p tessting
        testing = DATAtrim.performance{kk}';
        ttestRun=0;
        
        if size((testing),1)>2;
            ttestRun=1;
            k=1;
            for kkk = 1:size((testing),2)-1
                [h(kkk),p(kkk)] = ttest(testing(:,1),testing(:,kkk+1));
                if h(kkk)==1
                    sigPairs{k} = [1,kkk];
                    sigPairsPval(k) = [p(kkk)];
                    k = k+1;
                end
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
        bar(DATAtrim.performance{kk}(:,1), 0.5);
        hold on
        set(gca,'XTickLabel',DATAtrim.whiskerID(:,kk));
        ylabel('First Session Performance')
        currPlot=currPlot+1;
        
        %plot avergae of each pair.
        %calculate the min number of observations per data point
        subplot(plotRows,plotCols,currPlot);
        bar(nanmean(DATAtrim.performance{kk},2), 0.5);
        hold on
        errorbar(nanmean(DATAtrim.performance{kk},2),nanstd(DATAtrim.performance{kk},0,2),'.');
        if (ttestRun==1);
            sigstar(sigPairs,sigPairsPval);
        end
        set(gca,'XTickLabel',DATAtrim.whiskerID(:,kk));
        ylabel('average Performance')
        currPlot=currPlot+1;
        
        testing = DATAtrim.dPrime{kk}';
        ttestRun=0;
        
        if size((testing),1)>2;
            ttestRun=1;
            k=1;
            clear h p sigPairs sigPairsPval
            
            for kkk = 1:size((testing),2)-1
                [h(kkk),p(kkk)] = ttest(testing(:,1),testing(:,kkk+1));
                if h(kkk)==1
                    sigPairs{k} = [1,kkk+1];
                    sigPairsPval(k) = [p(kkk)];
                    k = k+1;
                end
            end
        end
        
        
        subplot(plotRows,plotCols,currPlot);
        bar(DATAtrim.dPrime{kk}(:,1), 0.5);
        set(gca,'XTickLabel',DATAtrim.whiskerID(:,kk));
        ylabel('First Session dPrime')
        hold on
        currPlot=currPlot+1;
        
        
        subplot(plotRows,plotCols,currPlot);
        bar(nanmean(DATAtrim.dPrime{kk},2), 0.5);
        hold on
        errorbar(nanmean(DATAtrim.dPrime{kk},2),nanstd(DATAtrim.dPrime{kk},0,2),'.');
        if (ttestRun==1);
            sigstar(sigPairs,sigPairsPval);
        end
        set(gca,'XTickLabel',DATAtrim.whiskerID(:,kk));
        ylabel('average dPrime')
        currPlot=currPlot+1;
        
        %plot performance and d prime over session
        plotRows = 2;
        plotCols = 1;
        
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
        plot(DATAtrim.performanceAll{kk},'o-k','LineWidth',5);
        hold on
        set(gca,'XTick',[1:length(DATAtrim.performanceAll{kk})]);
        set(gca,'XTickLabel',DATAtrim.whiskerIDAll{kk});
        ylim([0 1]);
        set(gca, 'Ylim',[0 1])
        ylabel('Performance');
        xlabel('whisker trimming');
        currPlot=currPlot+1;
        plot([1 max(xlim)],percentCorrectChance,'k--','LineWidth',2)
        
        subplot(plotRows,plotCols,currPlot);
        plot(DATAtrim.dPrimeAll{kk},'o-k','LineWidth',5);
        hold on
        set(gca,'XTick',[1:length(DATAtrim.dPrimeAll{kk})]);
        set(gca,'XTickLabel',DATAtrim.whiskerIDAll{kk});
        plot([1 max(xlim)],dprimeThreshold,'k--','LineWidth',2)
        ylimit = ylim;
        
        if ylimit(1)<0;
            yLimitSet = ylimit(1);
        else
            yLimitSet = 0;
        end
        ylim([yLimitSet ylimit(2)]);
        set(gca, 'Ylim',[yLimitSet max(DATAtrim.dPrimeAll{kk})])
        ylabel('d prime');
        xlabel('whisker trimming');
        currPlot=currPlot+1;
        
        saveas(gca,fullfile('C:\Users\adesniklab\Documents\BehaviorRawData\currFigs\trimming',tempTitle),'jpeg');
        %Plot dprime over angles
        %delete empty cells
        
        plotRows = 3;
        plotCols = size((DATAtrim.pairedAnglesdPrime{kk}),2);
        
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
        for jj = 1:size((DATAtrim.pairedAnglesdPrime{kk}),2);
            plot(DATAtrim.pairedAnglesdPrime{kk}(:,jj),'o-','LineWidth',5,'Color',[0,1/(jj),1/(jj)]);
            templegend{jj} = ['diffAngle ',num2str(DATAtrim.pairedAnglesDiff{kk}(1,jj))];
            hold on
        end
        legend(templegend)
        
        set(gca,'XTick',(1:length(DATAtrim.pairedAnglesdPrime{kk})));
        set(gca,'XTickLabel',DATAtrim.whiskerIDAll{kk});
        
        maxY = max(DATAtrim.pairedAnglesdPrime{kk}(:));
        minY = min(DATAtrim.pairedAnglesdPrime{kk}(:));
        for jj = 1:size((DATAtrim.pairedAnglesdPrime{kk}),2);
            subplot(plotRows,plotCols,currPlot);
            plot(DATAtrim.pairedAnglesdPrime{kk}(:,jj),'o-','LineWidth',5,'Color',[0,1/(jj),1/(jj)]);
            hold on
            set(gca,'XTick',[1:length(DATAtrim.pairedAnglesdPrime{kk})]);
            set(gca,'XTickLabel',DATAtrim.whiskerIDAll{kk});
            currPlot = currPlot+1;
            set(gca,'Ylim',[minY maxY]);
            xlimit = xlim;
            ylimit = ylim;
            str = strcat('diffAngle ',templegend(jj));
            text(xlimit(1),(maxY+0.5),str)
        end
        saveas(gca,fullfile('C:\Users\adesniklab\Documents\BehaviorRawData\currFigs\trimming',tempTitle),'jpeg');
    end
end
end



function [DATAavg] = JB_plotSelectionPerformance(basicPropertiesToPlot,possibleAngles,plotON,selection,average,individualTraces,stim)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%eg
%[DATAavg] = JB_plotSelectionPerformance(basicPropertiesToPlot,possibleAngles,1,[27 29 31],1,0,1)
%stim=0 all sessions, stim=1 ( stimsession on opto days), stim =2 (control sessions on opto days
%
%plot sessionType performance to get an idea of how inclinded to lick the
%mouse is
DATAavg = [];
positionGraph2 = [602   83   644   869];
positionGraph3 = [995   676   875   309];
percentCorrectChance = [0.5 0.5];
dprimeThreshold = [1 1];
% selection = [23 24];

tempTitle = basicPropertiesToPlot{1,1}.mouseID;
tempTitle(findstr(tempTitle,'_'))=[];

plotRows = 4;
plotCols = 4;
plotTally = (plotRows*plotCols);
numFigs = 1;

if exist('selection','var')
    numPoints = selection;
else  
    numPoints = 1:1:length(basicPropertiesToPlot);
end

if (exist('plotON','var') && plotON==1)
    ff=figure;clf
    set(ff,'Position',positionGraph2);
    set(gcf,'name',tempTitle,'numbertitle','off')
    if (average==1)
        fff=figure;clf
        set(fff,'Position',positionGraph3);
        set(gcf,'name',tempTitle,'numbertitle','off')
    end
else
    figure('Visible','off');clf;
end

currPlot = 1;
addtoLegend = 1;
addedWater = 0;
addedNeg = 0;
addedOpto = 0;
legendAdd = [];

for h = 1:length(numPoints)
    
    if ((basicPropertiesToPlot{numPoints(h),1}.optogenetics==1) && (stim==1))
        saved = 0;
        activeAngles = cell2mat(basicPropertiesToPlot{numPoints(h),1}.performanceSTIM);
        probLick = cell2mat(basicPropertiesToPlot{numPoints(h),1}.probLickingSTIM);
        trialTypecombo = [basicPropertiesToPlot{numPoints(h),1}.HITStim basicPropertiesToPlot{numPoints(h),1}.MISSStim;basicPropertiesToPlot{numPoints(h),1}.FAStim basicPropertiesToPlot{numPoints(h),1}.CRStim];
    elseif ((basicPropertiesToPlot{numPoints(h),1}.optogenetics==1) && (stim==2))
        saved = 0;
        activeAngles = cell2mat(basicPropertiesToPlot{numPoints(h),1}.performanceNoSTIM);
        probLick = cell2mat(basicPropertiesToPlot{numPoints(h),1}.probLickingNoSTIM);
        trialTypecombo = [basicPropertiesToPlot{numPoints(h),1}.HITNoStim basicPropertiesToPlot{numPoints(h),1}.MISSNoStim;basicPropertiesToPlot{numPoints(h),1}.FANoStim basicPropertiesToPlot{numPoints(h),1}.CRNoStim];  
    else
        saved = 0;
        activeAngles = cell2mat(basicPropertiesToPlot{numPoints(h),1}.performance);
        probLick = cell2mat(basicPropertiesToPlot{numPoints(h),1}.probLicking);
        trialTypecombo = [basicPropertiesToPlot{numPoints(h),1}.HIT basicPropertiesToPlot{numPoints(h),1}.MISS;basicPropertiesToPlot{numPoints(h),1}.FA basicPropertiesToPlot{numPoints(h),1}.CR];
    end
    
    plotAngles = possibleAngles-270;
    [~,c] = find(isnan(activeAngles));
    activeAngles(c) = [];
    plotAngles(c) = [];
    probLick(c) = [];
    
    if(length(activeAngles)>2) %if more than 2 angles were presented
        figure(ff);
        subplot(plotRows,plotCols,currPlot);
        plot(plotAngles,activeAngles,'ok-','MarkerFaceColor','k','LineWidth',2, 'MarkerSize',3);
        hold on
        currPlot=currPlot+1;
        xlabel('Angles');
        ylabel('Performance');
        ylim([0 1])
        xlimit = xlim;
        text(xlimit(1),1.05,basicPropertiesToPlot{numPoints(h),1}.namedata)
        str = strcat('d'' ', num2str(basicPropertiesToPlot{numPoints(h),1}.dprime));
        text(xlimit(1),0.2,num2str(h))
        text(xlimit(1),0.1,str)
        
        if(basicPropertiesToPlot{numPoints(h),1}.negReinforcer==1)
            squareLeg1 =  plot(xlimit(2),0.1,'rs','MarkerFaceColor','r', 'MarkerSize',8);
            if (addedNeg==0)
                legendAdd(addtoLegend,:) = squareLeg1;
                legendTab{addtoLegend} = 'neg reinforcement';
                addedNeg=1;
                addtoLegend=addtoLegend+1;
            end
        end
        
        if(basicPropertiesToPlot{numPoints(h),1}.waterSchedule==1)
            squareLeg2 = plot(xlimit(2),0.2,'gs','MarkerFaceColor','b', 'MarkerSize',8);
            if (addedOpto==0)
                legendAdd(addtoLegend,:) = squareLeg2;
                legendTab{addtoLegend} = 'water schedule';
                addedOpto=1;
                addtoLegend=addtoLegend+1;
            end
        end
        
        if(basicPropertiesToPlot{numPoints(h),1}.optogenetics==1)
            squareLeg3 = plot(xlimit(2),0.3,'bs','MarkerFaceColor','g', 'MarkerSize',8);
            if (addedWater==0)
                legendAdd(addtoLegend,:) = squareLeg3;
                legendTab{addtoLegend} = 'optogenetics';
                addedWater=1;
                addtoLegend=addtoLegend+1;
            end
        end
        
        figure(ff);
        subplot(plotRows,plotCols,currPlot);
        plot(plotAngles,probLick,'o-k','MarkerFaceColor','k','LineWidth',2, 'MarkerSize',3);
        currPlot=currPlot+1;
        xlabel('Angles');
        ylabel('Licking Probability');
        ylim([-0.2 1.2])
        
        subplot(plotRows,plotCols,currPlot)
        bar(trialTypecombo, 'stacked');
        hold on
        XlableAxis = {'HIT/MISS'; 'FA/CR'};
        set(gca,'XTickLabel',XlableAxis);
        currPlot=currPlot+1;
        
        subplot(plotRows,plotCols,currPlot);
        
        if ((basicPropertiesToPlot{numPoints(h),1}.optogenetics==1) && (stim==1))
            dPrimeAngles = basicPropertiesToPlot{numPoints(h),1}.pairsDprimeSTIM;
            
        elseif ((basicPropertiesToPlot{numPoints(h),1}.optogenetics==1) && (stim==2))
            dPrimeAngles = basicPropertiesToPlot{numPoints(h),1}.pairsDprimeNoSTIM;
        else
            dPrimeAngles = basicPropertiesToPlot{numPoints(h),1}.pairsDprime;
        end
        dPrimeAngles(any(isnan(dPrimeAngles),2),:)=[];
        plot(basicPropertiesToPlot{numPoints(h),1}.pairsDiff(1:length(dPrimeAngles)),dPrimeAngles,'o-k','MarkerFaceColor','k','LineWidth',2, 'MarkerSize',3)
        xlabel('Angle Pairs');
        ylabel('dPrime');
        set(gca, 'Xdir', 'reverse');
        ylim([0 length(dPrimeAngles)])
        currPlot=currPlot+1;
        
        if rem(currPlot,plotTally)==1 %if the value is divisable by 4 - open a new plot
            baseFileName = strcat(tempTitle,num2str(numFigs)); %save old figure
            if(length(legendAdd))>0
                hLL = legend([legendAdd(1:length(legendAdd))], legendTab{:});
                newPosition = [0 0 0.2 0.1];
                set(hLL, 'Position', newPosition, 'Box', 'off')
                addtoLegend = 1;
                addedWater = 0;
                addedNeg = 0;
                addedOpto = 0;
                clear legendTab legendAdd;
            end
            saveas(gca,fullfile('C:\Users\adesniklab\Documents\BehaviorRawData\currFigs\psychometricCurves',baseFileName),'jpeg');
            saved = 1;
            if (plotON==1)
                ff=figure;clf
                set(ff,'Position',positionGraph2);
            else
                figure('Visible','off');clf;
            end
            numFigs = numFigs+1;
            currPlot = 1;
            legendAdd = [];
        end
        
        if ((average==1) && ~isempty(selection));
               DATAavg.dprime(h,:) = basicPropertiesToPlot{numPoints(h),1}.dprime;
                 DATAavg.sessionperformance(h,:) = basicPropertiesToPlot{numPoints(h),1}.sessionperformance;
            DATAavg.activeAngles(h,:) = activeAngles;
            DATAavg.plotAngles(h,:) = plotAngles;
            DATAavg.probLick(h,:) = probLick;
            DATAavg.dPrimeAngles(h,:) = dPrimeAngles;
            DATAavg.dPrimeAnglePairs(h,:) = basicPropertiesToPlot{numPoints(h),1}.pairsDiff(1:length(dPrimeAngles));
            DATAavg.HITRate(h,:) = basicPropertiesToPlot{numPoints(h),1}.HIT/(basicPropertiesToPlot{numPoints(h),1}.HIT+basicPropertiesToPlot{numPoints(h),1}.MISS);
             DATAavg.FArate(h,:) = basicPropertiesToPlot{numPoints(h),1}.FA/(basicPropertiesToPlot{numPoints(h),1}.FA+basicPropertiesToPlot{numPoints(h),1}.CR);
        end
    end
    
    if saved==0
        if currPlot>1
            baseFileName = strcat(tempTitle,num2str(numFigs)); %save old figure
            saveas(gca,fullfile('C:\Users\adesniklab\Documents\BehaviorRawData\currFigs\psychometricCurves',baseFileName),'jpeg');
        end
    end
end

if ((average==1) && ~isempty(selection));
    avgLineWidth = 3;
    figure(fff);
    subplot(1,3,1);
    if(individualTraces==1)
        plot(DATAavg.plotAngles',DATAavg.activeAngles','-','Color',[0.4 0.4 0.4])
        hold on
    end
    errorbar(mean(DATAavg.plotAngles)',mean(DATAavg.activeAngles)',std(DATAavg.activeAngles),'k','Linewidth',avgLineWidth)
    hold on
    plot([min(xlim) max(xlim)],percentCorrectChance,'k--','LineWidth',2)
    xlabel('Angles');
    ylabel('Performance');
    ylim([0 1])
    xlimit = xlim;
    
    subplot(1,3,2);
    if(individualTraces==1)
        plot(DATAavg.plotAngles',DATAavg.probLick','-','Color',[0.4 0.4 0.4])
        hold on
    end
    errorbar(mean(DATAavg.plotAngles)',mean(DATAavg.probLick)',std(DATAavg.probLick),'k','Linewidth',avgLineWidth)
    hold on
    plot([min(xlim) max(xlim)],percentCorrectChance,'k--','LineWidth',2)
    xlabel('Angles');
    ylabel('Licking Probability');
    ylim([-0.2 1.2])
    
    subplot(1,3,3);
    if(individualTraces==1)
        plot(DATAavg.dPrimeAnglePairs',DATAavg.dPrimeAngles','-','Color',[0.4 0.4 0.4])
        hold on
    end
    errorbar(mean(DATAavg.dPrimeAnglePairs)',mean(DATAavg.dPrimeAngles)',std(DATAavg.dPrimeAngles),'k','Linewidth',avgLineWidth)
    hold on
    plot([1 max(xlim)],dprimeThreshold,'k--','LineWidth',2)
    ylim([-1 4])
    xlabel('Angle Pairs');
    ylabel('dPrime');
    set(gca, 'Xdir', 'reverse');
    saveas(gca,fullfile('C:\Users\adesniklab\Documents\BehaviorRawData\currFigs\avgPsychometricCurves',baseFileName),'jpeg');
    
end

end
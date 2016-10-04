function [] = JB_plotSessionPerformance(basicPropertiesToPlot,possibleAngles,plotON)
% JB_plotSessionPerformance plots performance and d' for individual
% sessions

%   [DATAtrim] = JB_plotSessionPerformance(basicPropertiesToPlot) - data
%   inputed - returned from JB_basicBehaviorProperties
%   possibleAngles - returned form JB_basicBehaviorProperties of all
%   posible angles
%   plotON = 1 (plots output), = 0 (no plot generated)

% Examples:
%   [] = JB_plotSessionPerformance(basicPropertiesToPlot,possibleAngles,0);
%

positionGraph2 = [602   83   644   869];

plotRows = 4;
plotCols = 4;
plotTally = (plotRows*plotCols);
numFigs = 1;
numPoints = 1:1:length(basicPropertiesToPlot);

if (plotON==1)
    ff=figure;clf
    set(ff,'Position',positionGraph2);
else
    figure('Visible','off');clf;
end

tempTitle = basicPropertiesToPlot{1,1}.mouseID;
tempTitle(findstr(tempTitle,'_'))=[];
set(gcf,'name',tempTitle,'numbertitle','off')
currPlot = 1;

addtoLegend = 1;
addedWater = 0;
addedNeg = 0;
addedOpto = 0;
legendAdd = [];

for h = 1:length(numPoints)
    saved = 0;
    activeAngles = cell2mat(basicPropertiesToPlot{h,1}.performance);
    probLick = cell2mat(basicPropertiesToPlot{h,1}.probLicking);
    plotAngles = possibleAngles-270;
    [~,c] = find(isnan(activeAngles));
    activeAngles(c) = [];
    plotAngles(c) = [];
    probLick(c) = [];
    trialTypecombo = [basicPropertiesToPlot{h,1}.HIT basicPropertiesToPlot{h,1}.MISS;basicPropertiesToPlot{h,1}.FA basicPropertiesToPlot{h,1}.CR];
    
    if(length(activeAngles)>2) %if more than 2 angles were presented
        subplot(plotRows,plotCols,currPlot);
        plot(plotAngles,activeAngles,'o-');
        hold on
        currPlot=currPlot+1;
        xlabel('Angles');
        ylabel('Performance');
        ylim([0 1])
        xlimit = xlim;
        text(xlimit(1),1.05,basicPropertiesToPlot{h,1}.namedata)
        str = strcat('d'' ', num2str(basicPropertiesToPlot{h,1}.dprime));
        text(xlimit(1),0.1,str)
        
        if(basicPropertiesToPlot{h,1}.negReinforcer==1)
            squareLeg1 =  plot(xlimit(2),0.1,'rs','MarkerFaceColor','r', 'MarkerSize',8);
            if (addedNeg==0)
                legendAdd(addtoLegend,:) = squareLeg1;
                legendTab{addtoLegend} = 'neg reinforcement';
                addedNeg=1;
                addtoLegend=addtoLegend+1;
            end
        end
        
        if(basicPropertiesToPlot{h,1}.waterSchedule==1)
            squareLeg2 = plot(xlimit(2),0.2,'gs','MarkerFaceColor','b', 'MarkerSize',8);
            
            if (addedOpto==0)
                legendAdd(addtoLegend,:) = squareLeg2;
                legendTab{addtoLegend} = 'water schedule';
                addedOpto=1;
                addtoLegend=addtoLegend+1;
            end
        end
        
        if(basicPropertiesToPlot{h,1}.optogenetics==1)
            squareLeg3 = plot(xlimit(2),0.3,'bs','MarkerFaceColor','g', 'MarkerSize',8);
            if (addedWater==0)
                legendAdd(addtoLegend,:) = squareLeg3;
                legendTab{addtoLegend} = 'optogenetics';
                addedWater=1;
                addtoLegend=addtoLegend+1;
            end
        end
        
        subplot(plotRows,plotCols,currPlot);
        plot(plotAngles,probLick,'o-');
        currPlot=currPlot+1;
        xlabel('Angles');
        ylabel('Licking Probability');
          xlimit = xlim;
        text(xlimit(1),0.2,num2str(h))
        ylim([0 1])
        
        subplot(plotRows,plotCols,currPlot)
        bar(trialTypecombo, 'stacked');
        hold on
        XlableAxis = {'HIT/MISS'; 'FA/CR'};
        set(gca,'XTickLabel',XlableAxis);
        currPlot=currPlot+1;
        
        subplot(plotRows,plotCols,currPlot);
        dPrimeAngles = basicPropertiesToPlot{h,1}.pairsDprime;
        dPrimeAngles(isnan(dPrimeAngles))=[];

        
        plot(basicPropertiesToPlot{h,1}.pairsDiff(1:length(dPrimeAngles)),dPrimeAngles,'*-')
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
    end
    
    if saved==0
        if currPlot>1
            baseFileName = strcat(tempTitle,num2str(numFigs)); %save old figure
            saveas(gca,fullfile('C:\Users\adesniklab\Documents\BehaviorRawData\currFigs\psychometricCurves',baseFileName),'jpeg');
        end
    end
end
end


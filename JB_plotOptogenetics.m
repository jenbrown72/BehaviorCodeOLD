function [] = JB_plotOptogenetics(basicPropertiesToPlot,possibleAngles,plotON)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    positionGraph3 = [1268   83   644   869];
    positionGraph4 = [32   515   354   438];
    
    plotRows = 4;
    plotCols = 3;
    plotTally = (plotRows*plotCols);
    numFigs = 1;
        numPoints = 1:1:length(basicPropertiesToPlot);
    
    if (plotON==1)
        fff=figure;clf
        set(fff,'Position',positionGraph3);
    else
        figure('Visible','off');clf;
    end
    
        
    tempTitle = basicPropertiesToPlot{1,1}.mouseID;
    tempTitle(findstr(tempTitle,'_'))=[];
    set(gcf,'name',tempTitle,'numbertitle','off')
    currPlot = 1;
    tn=1;
    tempdPrime = nan(length(numPoints),2);
    tempdPerformance = nan(length(numPoints),2);
    
    for h = 1:length(numPoints)

        if (basicPropertiesToPlot{h,1}.optogenetics)==1
            activeAnglesSTIM = cell2mat(basicPropertiesToPlot{h,1}.performanceSTIM);
            activeAnglesnoSTIM = cell2mat(basicPropertiesToPlot{h,1}.performanceNoSTIM);
            probLickSTIM = cell2mat(basicPropertiesToPlot{h,1}.probLickingSTIM);
            probLickNoSTIM = cell2mat(basicPropertiesToPlot{h,1}.probLickingNoSTIM);
            plotAngles = possibleAngles;
            [~,c] = find(isnan(activeAnglesSTIM));
            activeAnglesSTIM(c) = [];
            activeAnglesnoSTIM(c) = [];
            plotAngles(c) = [];
            probLickSTIM(c) = [];
            probLickNoSTIM(c) = [];
            
            subplot(plotRows,plotCols,currPlot);
            line1 = plot(plotAngles,activeAnglesSTIM,'o-r');
            hold on;
            line2 = plot(plotAngles,activeAnglesnoSTIM,'o-k');
            currPlot=currPlot+1;
            xlabel('Angles');
            ylabel('Performance');
            xlimit = xlim;
            ylim([0 1])
            text(xlimit(1),1.05,basicPropertiesToPlot{h,1}.namedata)
            str = ['d'' C/S',' ', num2str(basicPropertiesToPlot{h,1}.dprimeNoSTIM),' / ', num2str(basicPropertiesToPlot{h,1}.dprimeSTIM)];
            text(xlimit(1),0.1,str)
            
            tempdPrime(tn,:) = [basicPropertiesToPlot{h,1}.dprimeNoSTIM, basicPropertiesToPlot{h,1}.dprimeSTIM];
            tempdPerformance(tn,:) = [basicPropertiesToPlot{h,1}.sessionperformanceNoSTIM, basicPropertiesToPlot{h,1}.sessionperformanceSTIM];
            tn=tn+1;
            
            subplot(plotRows,plotCols,currPlot);
            plot(plotAngles,probLickSTIM,'o-r');
            hold on;
            plot(plotAngles,probLickNoSTIM,'o-k');
            currPlot=currPlot+1;
            xlabel('Angles');
            ylabel('Licking Probability');
            ylim([0 1])
            
            subplot(plotRows,plotCols,currPlot);
            plot(basicPropertiesToPlot{h,1}.pairsDiff(1:length(basicPropertiesToPlot{h,1}.pairsDprimeSTIM)),basicPropertiesToPlot{h,1}.pairsDprimeSTIM,'o-r');
            hold on;
            plot(basicPropertiesToPlot{h,1}.pairsDiff(1:length(basicPropertiesToPlot{h,1}.pairsDprimeNoSTIM)),basicPropertiesToPlot{h,1}.pairsDprimeNoSTIM,'o-k');
            currPlot=currPlot+1;
            set(gca, 'Xdir', 'reverse')
            xlabel('Diff Angles');
            ylabel('dPrime');
            
        end
        
        if (currPlot>1)
            if rem(currPlot,plotTally)==1 %if the value is divisable by 4 - open a new plot
                baseFileName = strcat(tempTitle,num2str(numFigs)); %save old figure
                hL = legend([line1, line2],{'Stimulated', 'Control'});
                newPosition = [0 0 0.2 0.1];
                set(hL, 'Position', newPosition, 'Box', 'off')
                
                saveas(gca,fullfile('C:\Users\adesniklab\Documents\BehaviorRawData\currFigs\optogenetics',baseFileName),'jpeg');
                if (plotON==1)
                    fff=figure;clf
                    set(fff,'Position',positionGraph3);
                else
                    figure('Visible','off');clf;
                end
                numFigs = numFigs+1;
                currPlot = 1;
            end
        end
    end

        if (plotON==1)
        if currPlot==1
            close(fff)
        end
    end
    
    currPlot=1;
    
   tempdPrime =  tempdPrime(~any(isnan(tempdPrime),2),:);
   tempdPerformance =  tempdPerformance(~any(isnan(tempdPerformance),2),:);
    
    if any(tempdPrime)
        
        plotRows = 2;
        plotCols = 1;
        numFigs = 1;
        
        if (plotON==1)
            ffff=figure;clf
            set(ffff,'Position',positionGraph4);
        else
            figure('Visible','off');clf;
        end
        
        meandPrime = mean(tempdPrime);
        SEMdPrime = std(tempdPrime)/sqrt(length(tempdPrime));
        [~,p] = ttest(tempdPrime(:,1),tempdPrime(:,2));
        
        subplot(plotRows,plotCols,currPlot);
        plot(tempdPrime','k-')
        hold on
        errorbar(meandPrime,SEMdPrime,'r-','LineWidth',5)
        set(gca,'XTick',[1:2]);
        set(gca,'XTickLabel',['Cntr';'Stim']);
        xlimit = [0.5 2.5];
        set(gca, 'Xlim',[xlimit(1) xlimit(2)])
        ylimit = ylim;
        ylim([0 ylimit(2)]);
        set(gca, 'Ylim',[0 ylimit(2)])
        str = ['P = ', num2str(p)];
        ylabel('d prime');
        text(xlimit(1),0.2,str)
        currPlot=currPlot+1;
        
        
        meandPrime = nanmean(tempdPerformance);
        SEMdPrime = nanstd(tempdPerformance)/sqrt(length(tempdPerformance));
        [~,p] = ttest(tempdPerformance(:,1),tempdPerformance(:,2));
        
        subplot(plotRows,plotCols,currPlot);
        plot(tempdPerformance','k-')
        hold on
        errorbar(meandPrime,SEMdPrime,'r-','LineWidth',5)
        set(gca,'XTick',[1:2]);
        set(gca,'XTickLabel',['Cntr';'Stim']);
        xlimit = [0.5 2.5];
        set(gca, 'Xlim',[xlimit(1) xlimit(2)])
        ylimit = ylim;
        ylim([0 ylimit(2)]);
        set(gca, 'Ylim',[0 ylimit(2)])
        str = ['P = ', num2str(p)];
        ylabel('Performance');
        text(xlimit(1),0.2,str)
        
    end

end

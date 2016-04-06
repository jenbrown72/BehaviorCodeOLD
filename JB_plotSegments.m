function [] = JB_plotSegments(basicPropertiesToPlot,plotON)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
  %end
%     
    positionGraph7 = [754    78   327   869];
        percentCorrectChance = [0.5 0.5];
    dprimeThreshold = [1 1];
    plotRows = 2;
    plotCols = 1;
    numFigs = 1;
    currPlot=1;
    Segtally=1;
        numPoints = 1:1:length(basicPropertiesToPlot);

    if (plotON==1)
        fV = figure;clf
        set(fV,'Position',positionGraph7);
    else
        figure('Visible','off');clf;
    end
    
            tempTitle = basicPropertiesToPlot{1,1}.mouseID;
    tempTitle(findstr(tempTitle,'_'))=[];
    set(gcf,'name',tempTitle,'numbertitle','off');
    
    indSeg = nan(length(numPoints),1);
    for h = 1:length(numPoints)
        if ~isempty(basicPropertiesToPlot{h,1}.segmentPerformance)
            indSeg(h,1)=1;
            Segtally = Segtally+1;
        else
            indSeg(h,1)=0;
        end
    end
    
    

    for h = 1:length(numPoints)
        if (indSeg(h,1)==1)
            subplot(plotRows,plotCols,currPlot);
            plot(basicPropertiesToPlot{h,1}.segmentPerformance,'ok-','LineWidth',5)
            hold on
            
        plot([1 max(xlim)],percentCorrectChance,'k--','LineWidth',2)
         %       line([1 max(xlim)],percentCorrectThreshold,'k--','LineWidth',2)
            % SessionTypes = {'Cntr' ; 'Stim out'; 'Cntr'; 'Stim out'};
            %SessionTypes = {'Block 1' ; 'Stim out'; 'Cntr'; 'Stim out'};
            %NumTicks = length(SessionTypes);
            NumTicks = length(basicPropertiesToPlot{h,1}.segmentPerformance);
            L = get(gca,'XLim');
            set(gca,'XTick',linspace(L(1),L(2),NumTicks));
            currPlot = currPlot+1;
            
            ylim([0 1]);
            title('Performance over blocks','fontWeight','bold');
            xlabel('Block number','fontWeight','bold');
            ylabel('Percent Correct ((hit+CR)/totalTrials', 'fontWeight','bold');
            %set(gca,'XTickLabel',SessionTypes)
            subplot(plotRows,plotCols,currPlot);
            plot(basicPropertiesToPlot{h,1}.segmentdPrime,'ok-','LineWidth',5)
            hold on
                %    line([1 max(xlim)],dprimeChance,'k--','LineWidth',2)
              plot([1 max(xlim)],dprimeThreshold,'k--','LineWidth',2)
            % SessionTypes = {'Cntr' ; 'Stim out'; 'Cntr'; 'Stim out'};
            %SessionTypes = {'Block 1' ; 'Stim out'; 'Cntr'; 'Stim out'};
            %NumTicks = length(SessionTypes);
            NumTicks = length(basicPropertiesToPlot{h,1}.segmentdPrime);
            L = get(gca,'XLim');
            set(gca,'XTick',linspace(L(1),L(2),NumTicks));
            currPlot = currPlot+1;
            ylimit = ylim;

            if ylimit(1)<0;
                yLimitSet = ylimit(1);
            else
                yLimitSet = 0;
            end
            
            ylim([yLimitSet ylimit(2)]);
          %  set(gca, 'Ylim',[yLimitSet max(tempTrimdPrime)])
            title('Performance over blocks','fontWeight','bold');
            xlabel('Block number','fontWeight','bold');
            ylabel('d Prime', 'fontWeight','bold');
            text(1,ylimit(1)+0.1,basicPropertiesToPlot{h,1}.namedata)
            %set(gca,'XTickLabel',SessionTypes)
             baseFileName = strcat(tempTitle,num2str(numFigs)); %save old figure
             saveas(gca,fullfile('C:\Users\adesniklab\Documents\BehaviorRawData\currFigs\segmented',baseFileName),'jpeg');
numFigs = numFigs+1;
 currPlot=1;
  if (plotON==1)
        fV = figure;clf
        set(fV,'Position',positionGraph7);
    else
        figure('Visible','off');clf;
    end
        else
            indSeg(h,1)=0;
        end
    end
   % saveas(gca,fullfile('C:\Users\adesniklab\Documents\BehaviorRawData\currFigs\segmented',tempTitle),'jpeg');
    
   
    if (plotON==1)
        if currPlot==1
            close(fV)
        end
    end

end


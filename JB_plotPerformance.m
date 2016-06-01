function [] = JB_plotPerformance(basicPropertiesToPlot,plotON) 
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

positionGraph1 = [14   83   503   869];
dprimeChance = [1 1];
dprimeThreshold = [1 1];
percentCorrectChance = [0.5 0.5];
percentCorrectThreshold = [0.7 0.7];

    if (plotON==1)
        f = figure;clf
        set(f,'Position',positionGraph1);
    else
        figure('Visible','off');clf;
    end
    
    tempTitle = basicPropertiesToPlot{1,1}.mouseID;
    tempTitle(findstr(tempTitle,'_'))=[];
    set(gcf,'name',tempTitle,'numbertitle','off');
    
    noSubPlots = 4;
    subplot(noSubPlots,1,1)
    numPoints = 1:1:length(basicPropertiesToPlot);
    for j = 1:length(numPoints);
        subplot(noSubPlots,1,1)
        plot(numPoints(j),basicPropertiesToPlot{j}.totalStepsPerMin,'or','MarkerSize', 10,'MarkerFaceColor','r')
        hold on
    end
    
    ylabel('totalStepsPerMin');
    xlabel('Session Number');
    %
    % subplot(noSubPlots,1,2)
    % for j = 1:length(numPoints);
    %     plot(numPoints(j),basicPropertiesToPlot{j}.sessionDuration,'or','MarkerSize', 10,'MarkerFaceColor','r')
    %     hold on
    % end
    %
    % ylabel('sessionDuration');
    % xlabel('Session Number');
    
    %plot primary session type per day
    subplot(noSubPlots,1,2)
    SessionTypes = {'S1auto' ; 'S1'; 'S2'; 'S6'; 'S8'; 'S10'; 'S12'};
    
    for j = 1:length(numPoints);
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
    
    for j = 1:length(numPoints);
        negReinforcer(j,1) = basicPropertiesToPlot{j,1}.negReinforcer;
        waterSchedule(j,1) = basicPropertiesToPlot{j,1}.waterSchedule;
        optogenetics(j,1) = basicPropertiesToPlot{j,1}.optogenetics;
    end
    
    ColorCodeColors = [1.0;0.9;0.8;0.6;0.4;0.2;0];
    
    for j = 1:length(numPoints);
        if waterSchedule(j)==1;
            plot(numPoints(j),colorCode(j,1),'o','MarkerSize', 10, 'linewidth',2,'MarkerEdgeColor',[0,0,waterSchedule(j)], 'MarkerFaceColor',[0.4,ColorCodeColors(colorCode(j)),0.8],'Color', [0.4,ColorCodeColors(colorCode(j)),0.8])
        elseif negReinforcer(j)==1;
            plot(numPoints(j),colorCode(j,1),'o','MarkerSize', 10, 'linewidth',2,'MarkerEdgeColor',[negReinforcer(j),0,0], 'MarkerFaceColor',[0.4,ColorCodeColors(colorCode(j)),0.8],'Color', [0.4,ColorCodeColors(colorCode(j)),0.8])
        else
            plot(numPoints(j),colorCode(j,1),'o','MarkerSize', 10, 'linewidth',2,'MarkerEdgeColor',[0.4,ColorCodeColors(colorCode(j)),0.8], 'MarkerFaceColor',[0.4,ColorCodeColors(colorCode(j)),0.8],'Color', [0.4,ColorCodeColors(colorCode(j)),0.8])
        end
        hold on
        if optogenetics(j)==1;
            plot(numPoints(j),colorCode(j,1),'*y','MarkerSize', 5)
        end
    end
    set(gca,'YTick',[1:length(SessionTypes)]);
    set(gca,'YTickLabel',SessionTypes)
    ylabel('Session Type');
    xlabel('Session Number');
    
    %plot performance
    subplot(noSubPlots,1,3)
    addtoLegend = 1;
    addedWater = 0;
    addedNeg = 0;
    addedOpto = 0;
    
    legendAdd = [];
    legendTab = [];
    for jj = 1:length(numPoints);
%         
%     X = ['Date ', basicPropertiesToPlot{jj,1}.date(1:11), ' Session ', basicPropertiesToPlot{jj,1}.sessionType, ' Performance ', num2str( basicPropertiesToPlot{jj,1}.sessionperformance), ' dPrime ', num2str(basicPropertiesToPlot{jj,1}.dprime)];
%     disp(X)
        
        if waterSchedule(jj)==1;
            circle1 = plot(numPoints(jj),basicPropertiesToPlot{jj,1}.sessionperformance,'o','MarkerSize', 10, 'linewidth',2,'MarkerEdgeColor',[0,0,waterSchedule(jj)], 'MarkerFaceColor',[0.4,ColorCodeColors(colorCode(jj)),0.8],'Color', [0.4,ColorCodeColors(colorCode(jj)),0.8]);
            hold on
            if (addedWater==0)
                legendAdd(addtoLegend,:) = circle1;
                legendTab{addtoLegend} = 'water schedule';
                addedWater=1;
                addtoLegend=addtoLegend+1;
            end
        elseif negReinforcer(jj)==1;
            circle2 = plot(numPoints(jj),basicPropertiesToPlot{jj,1}.sessionperformance,'o','MarkerSize', 10, 'linewidth',2,'MarkerEdgeColor',[negReinforcer(jj),0,0], 'MarkerFaceColor',[0.4,ColorCodeColors(colorCode(jj)),0.8],'Color', [0.4,ColorCodeColors(colorCode(jj)),0.8]);
            hold on
            if (addedNeg==0)
                legendAdd(addtoLegend,:) = circle2;
                legendTab{addtoLegend} = 'neg reinforcement';
                addedNeg=1;
                addtoLegend=addtoLegend+1;
            end
        else
            plot(numPoints(jj),basicPropertiesToPlot{jj,1}.sessionperformance,'o','MarkerSize', 10, 'linewidth',2,'MarkerEdgeColor',[0.4,ColorCodeColors(colorCode(jj)),0.8], 'MarkerFaceColor',[0.4,ColorCodeColors(colorCode(jj)),0.8],'Color', [0.4,ColorCodeColors(colorCode(jj)),0.8]);
        end
        hold on
        
        if optogenetics(jj)==1;
            circle3 = plot(numPoints(jj),basicPropertiesToPlot{jj,1}.sessionperformance,'*y','MarkerSize', 5);
            hold on
            if (addedOpto==0)
                legendAdd(addtoLegend,:) = circle3;
                legendTab{addtoLegend} = 'optogenetics';
                addedOpto = 1;
                addtoLegend=addtoLegend+1;
            end
        end
        
        plot([0 max(numPoints)],percentCorrectChance,'k--')
       % plot([0 max(numPoints)],percentCorrectThreshold,'k--')
        ylim([0 1])
        xlim([0 length(numPoints)])
        ylabel('Performance');
        xlabel('Session Number');
    end
    
    
    %plot Dprime
    subplot(noSubPlots,1,4)
    addtoLegend = 1;
    addedWater = 0;
    addedNeg = 0;
    addedOpto = 0;
    
    legendAdd = [];
    legendTab = [];
    for jj = 1:length(numPoints);
        
        if waterSchedule(jj)==1;
            circle1 = plot(numPoints(jj),basicPropertiesToPlot{jj,1}.dprime,'o','MarkerSize', 10, 'linewidth',2,'MarkerEdgeColor',[0,0,waterSchedule(jj)], 'MarkerFaceColor',[0.4,ColorCodeColors(colorCode(jj)),0.8],'Color', [0.4,ColorCodeColors(colorCode(jj)),0.8]);
            hold on
            if (addedWater==0)
                legendAdd(addtoLegend,:) = circle1;
                legendTab{addtoLegend} = 'water schedule';
                addedWater=1;
                addtoLegend=addtoLegend+1;
            end
        elseif negReinforcer(jj)==1;
            circle2 = plot(numPoints(jj),basicPropertiesToPlot{jj,1}.dprime,'o','MarkerSize', 10, 'linewidth',2,'MarkerEdgeColor',[negReinforcer(jj),0,0], 'MarkerFaceColor',[0.4,ColorCodeColors(colorCode(jj)),0.8],'Color', [0.4,ColorCodeColors(colorCode(jj)),0.8]);
            hold on
            if (addedNeg==0)
                legendAdd(addtoLegend,:) = circle2;
                legendTab{addtoLegend} = 'neg reinforcement';
                addedNeg=1;
                addtoLegend=addtoLegend+1;
            end
        else
            plot(numPoints(jj),basicPropertiesToPlot{jj,1}.dprime,'o','MarkerSize', 10, 'linewidth',2,'MarkerEdgeColor',[0.4,ColorCodeColors(colorCode(jj)),0.8], 'MarkerFaceColor',[0.4,ColorCodeColors(colorCode(jj)),0.8],'Color', [0.4,ColorCodeColors(colorCode(jj)),0.8]);
        end
        hold on
        
        if optogenetics(jj)==1;
            circle3 = plot(numPoints(jj),basicPropertiesToPlot{jj,1}.dprime,'*y','MarkerSize', 5);
            hold on
            if (addedOpto==0)
                legendAdd(addtoLegend,:) = circle3;
                legendTab{addtoLegend} = 'optogenetics';
                addedOpto = 1;
                addtoLegend=addtoLegend+1;
            end
        end

    %    plot([0 max(numPoints)],dprimeChance,'k--')
         plot([0 max(numPoints)],dprimeThreshold,'k--')
        xlim([0 length(numPoints)])
        ylabel('dprime');
        xlabel('Session Number');
    end


    if(length(legendAdd))>0
        hLL = legend([legendAdd(1:length(legendAdd))], legendTab{:});
        newPosition = [0 0 0.2 0.1];
        set(hLL, 'Position', newPosition, 'Box', 'off')
    end
    
    baseFileName = strcat(tempTitle); %save old figure
    saveas(gca,fullfile('C:\Users\adesniklab\Documents\BehaviorRawData\currFigs\basicSessionProperties',baseFileName),'jpeg');
    

end


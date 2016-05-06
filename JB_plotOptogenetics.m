function [DATAstim] = JB_plotOptogenetics(basicPropertiesToPlot,possibleAngles,plotON,subset,nostim)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
positionGraph2 = [34   250   968   692];
positionGraph3 = [152   127   994   869];

if nargin>3
    toPlot = subset;
else
    toPlot = 1:length(basicPropertiesToPlot);
end

if nargin>4
    nostim = nostim;
else
    nostim = 0;
end

percentCorrectChance = [0.5 0.5];
dprimeThreshold = [1 1];
plotRows = 5;
plotCols = 4;
plotTally = (plotRows*plotCols);
numFigs = 1;
numPoints = 1:1:length(toPlot);
possibleAngles = [225;241;254;263;266;268;270;272;274;277;284;299;315];
stimulatedSession = 0;
saved=1;

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
DATAstim.dPrime = nan(length(numPoints),2);

for hh = 1:length(toPlot);
    
    h = toPlot(hh);
    %for h = 1:length(numPoints)
    clear meanLatency stdLatency
    
    if (((basicPropertiesToPlot{h,1}.optogenetics)==1) || (nostim==1))
        saved=0;
        stimulatedSession=1;
        basicPropertiesToPlot{h,1}.pairsDprimeSTIM(isnan(basicPropertiesToPlot{h,1}.pairsDprimeSTIM),:)=[];
        basicPropertiesToPlot{h,1}.pairsDprimeNoSTIM(isnan(basicPropertiesToPlot{h,1}.pairsDprimeNoSTIM),:)=[];
%               basicPropertiesToPlot{h,1}.pairsDprimeNoSTIM(~any(basicPropertiesToPlot{h,1}.pairsDprimeNoSTIM,2),:)=[];
        activeAnglesSTIM = cell2mat(basicPropertiesToPlot{h,1}.performanceSTIM);
        
        if (basicPropertiesToPlot{h,1}.optogenetics)==1
            activeAnglesnoSTIM = cell2mat(basicPropertiesToPlot{h,1}.performanceNoSTIM);
        else
            activeAnglesnoSTIM = cell2mat(basicPropertiesToPlot{h,1}.performance);
        end
        
        
        
        probLickSTIM = cell2mat(basicPropertiesToPlot{h,1}.probLickingSTIM);
        
        if (basicPropertiesToPlot{h,1}.optogenetics)==1
            probLickNoSTIM = cell2mat(basicPropertiesToPlot{h,1}.probLickingNoSTIM);
        else
            probLickNoSTIM = cell2mat(basicPropertiesToPlot{h,1}.probLicking);
        end
        
        plotAngles = possibleAngles-270;
        [~,c] = find(isnan(activeAnglesnoSTIM));
        activeAnglesSTIM(c) = [];
        activeAnglesnoSTIM(c) = [];
        plotAngles(c) = [];
        probLickSTIM(c) = [];
        probLickNoSTIM(c) = [];
        
        if basicPropertiesToPlot{h,1}.orgGOstim==1
            temphitSTIM = nansum(cell2mat(basicPropertiesToPlot{h, 1}.performanceSTIM(1:6)));
            tempFASTIM =  nansum(1-cell2mat(basicPropertiesToPlot{h, 1}.performanceSTIM(7:13)));
            temphitNoSTIM = nansum(cell2mat(basicPropertiesToPlot{h, 1}.performanceNoSTIM(1:6)));
            tempFANoSTIM =  nansum(1-cell2mat(basicPropertiesToPlot{h, 1}.performanceNoSTIM(7:13)));
        else
            temphitSTIM = nansum(cell2mat(basicPropertiesToPlot{h, 1}.performanceSTIM(7:13)));
            tempFASTIM =  nansum(1-cell2mat(basicPropertiesToPlot{h, 1}.performanceSTIM(1:6)));
            temphitNoSTIM = nansum(cell2mat(basicPropertiesToPlot{h, 1}.performanceNoSTIM(7:13)));
            tempFANoSTIM =  nansum(1-cell2mat(basicPropertiesToPlot{h, 1}.performanceNoSTIM(1:6)));
        end
        
        DATAstim.lickBias(tn,2) = (temphitSTIM+tempFASTIM)/sum(~isnan(activeAnglesSTIM),2);
        DATAstim.lickBias(tn,1) = (temphitNoSTIM+tempFANoSTIM)/sum(~isnan(activeAnglesSTIM),2);
        
        %plot session performance
        subplot(plotRows,plotCols,currPlot);
        line1 = plot(plotAngles,activeAnglesSTIM,'o-r','MarkerSize',3,'LineWidth',3);
        hold on;
        line2 = plot(plotAngles,activeAnglesnoSTIM,'o-k','MarkerSize',3,'LineWidth',3);
        currPlot=currPlot+1;
        xlabel('Angles');
        ylabel('Performance');
        plot([min(xlim) max(xlim)],percentCorrectChance,'k--','LineWidth',2)
        xlimit = xlim;
        ylim([0 1])
        text(xlimit(1),1.05,basicPropertiesToPlot{h,1}.namedata)
        if (basicPropertiesToPlot{h,1}.optogenetics)==1
            str = ['d'' C/S',' ', num2str(basicPropertiesToPlot{h,1}.dprimeNoSTIM),' / ', num2str(basicPropertiesToPlot{h,1}.dprimeSTIM)];
        else
            str = ['d'' C',' ', num2str(basicPropertiesToPlot{h,1}.dprime)];
        end
        text(xlimit(1),0.1,str)
        text(xlimit(1),0.2,num2str(h))
        if (basicPropertiesToPlot{h,1}.optogenetics)==1
            DATAstim.dPrime(tn,:) = [basicPropertiesToPlot{h,1}.dprimeNoSTIM, basicPropertiesToPlot{h,1}.dprimeSTIM];
            DATAstim.performance(tn,:) = [basicPropertiesToPlot{h,1}.sessionperformanceNoSTIM, basicPropertiesToPlot{h,1}.sessionperformanceSTIM];
        end
        %plot session lick probability
        subplot(plotRows,plotCols,currPlot);
        plot(plotAngles,probLickSTIM,'o-r','MarkerSize',3,'LineWidth',3);
        hold on;
        plot(plotAngles,probLickNoSTIM,'o-k','MarkerSize',3,'LineWidth',3);
        currPlot=currPlot+1;
        xlabel('Angles');
        ylabel('Licking Probability');
        plot([min(xlim) max(xlim)],percentCorrectChance,'k--','LineWidth',2)
        ylim([0 1])
        
        %plot session dprime
        subplot(plotRows,plotCols,currPlot);
        
        plot(basicPropertiesToPlot{h,1}.pairsDiff(1:length(basicPropertiesToPlot{h,1}.pairsDprimeSTIM)),basicPropertiesToPlot{h,1}.pairsDprimeSTIM,'o-r','MarkerSize',3,'LineWidth',3);
        hold on;
          if (basicPropertiesToPlot{h,1}.optogenetics)==1
        plot(basicPropertiesToPlot{h,1}.pairsDiff(1:length(basicPropertiesToPlot{h,1}.pairsDprimeNoSTIM)),basicPropertiesToPlot{h,1}.pairsDprimeNoSTIM,'o-k','MarkerSize',3,'LineWidth',3);
        
          else
                      plot(basicPropertiesToPlot{h,1}.pairsDiff(1:length(basicPropertiesToPlot{h,1}.pairsDprime)),basicPropertiesToPlot{h,1}.pairsDprime,'o-k','MarkerSize',3,'LineWidth',3);

          end
          
          currPlot=currPlot+1;
        plot([min(xlim) max(xlim)],dprimeThreshold,'k--','LineWidth',2);
        set(gca, 'Xdir', 'reverse')
        xlabel('Diff Angles');
        ylabel('dPrime');
        Ylimit = ylim;
        Ylimit(Ylimit(1)>0)=0;
        Ylimit(Ylimit(2)<3,2) = 3;
        set(gca,'Ylim',[Ylimit]);
        
        %plot session angle dprimes
        subplot(plotRows,plotCols,currPlot);
        for jj = 1:size((basicPropertiesToPlot{h,1}.pairsDprimeSTIM),1);
            tempName = ['anglesdPrime',num2str(basicPropertiesToPlot{h,1}.pairsDiff(jj))];
            angleName{jj} = tempName;
            DATAstim.angles.(tempName){tn,1} = basicPropertiesToPlot{h,1}.pairsDprimeNoSTIM(jj);
            DATAstim.angles.(tempName){tn,2} = basicPropertiesToPlot{h,1}.pairsDprimeSTIM(jj);
            
        end
        
        for gg=1:length(basicPropertiesToPlot{h,1}.angleLick)
            meanLatency(gg,1:2) = [nanmean(basicPropertiesToPlot{h,1}.angleLickLatency{gg,1}),nanmean(basicPropertiesToPlot{h,1}.angleLickLatency{gg,2})];
            stdLatency(gg,1:2) = [nanstd(basicPropertiesToPlot{h,1}.angleLickLatency{gg,1}),nanstd(basicPropertiesToPlot{h,1}.angleLickLatency{gg,2})];
        end
        
        %plot lick latency
        DATAstim.latencysDiff(:,tn) = diff(meanLatency');
        tempangles = basicPropertiesToPlot{h,1}.angleLick-270;
        plot(tempangles,meanLatency(:,1),'o-k','MarkerSize',3,'LineWidth',3)
        hold on
        plot(tempangles,meanLatency(:,2),'o-r','MarkerSize',3,'LineWidth',3)
        errorbar(tempangles,meanLatency(:,1),stdLatency(:,1),'k','Linestyle','none');
        errorbar(tempangles,meanLatency(:,2),stdLatency(:,2),'r','Linestyle','none');
        
        NumTicks = length(tempangles);
        L = [tempangles(1) tempangles(end)];
        % L = get(gca,'XLim');
        labels = num2str(basicPropertiesToPlot{h,1}.angleLick-270);
        set(gca,'XTick',tempangles,'XTickLabel',labels,'XAxisLocation','top','ylim',[0 1000])
        view(-90,90)
        xlabel('LickLatency');
        ylabel('angles');
        currPlot=currPlot+1;
        
        %         %plot session angle trial counts
        %          subplot(plotRows,plotCols,currPlot);
        %         angleCountSTIM =  basicPropertiesToPlot{h,1}.noTrialsSTIM;
        %         angleCountNoSTIM =  basicPropertiesToPlot{h,1}.noTrialsNoSTIM;
        %         idx = isnan(angleCountSTIM);
        %         anglePlot=possibleAngles;
        %         angleCountSTIM(idx)=[];
        %         angleCountNoSTIM(idx)=[];
        %         anglePlot(idx)=[];
        %
        %         DATAstim.angleCountSTIM(:,tn) = angleCountSTIM;
        %         DATAstim.angleCountNoSTIM(:,tn) = angleCountNoSTIM;
        
        %         plot(anglePlot,angleCountSTIM,'o-r','MarkerSize',3,'LineWidth',3);
        %         hold on;
        %         plot(anglePlot,angleCountNoSTIM,'o-k','MarkerSize',3,'LineWidth',3);
        %         currPlot=currPlot+1;
        %         xlabel('Angles');
        %         ylabel('trial Count');
        %         ylimit = ylim;
        %         set(gca,'XAxisLocation','top','ylim',[0 ylimit(2)+5])
        %         view(-90,90)
        tn=tn+1;
    end
    
    if (currPlot>1)
        if rem(currPlot,plotTally)==1 %if the value is divisable by 4 - open a new plot
            baseFileName = strcat(tempTitle,num2str(numFigs)); %save old figure
            hL = legend([line1, line2],{'Stimulated', 'Control'});
            newPosition = [0 0 0.2 0.1];
            set(hL, 'Position', newPosition, 'Box', 'off')
            
            saveas(gca,fullfile('C:\Users\adesniklab\Documents\BehaviorRawData\currFigs\optogenetics',baseFileName),'jpeg');
            saved=1;
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
    if saved==0
        if currPlot>1
            baseFileName = strcat(tempTitle,num2str(numFigs)); %save old figure
            saveas(gca,fullfile('C:\Users\adesniklab\Documents\BehaviorRawData\currFigs\optogenetics',baseFileName),'jpeg');
        end
    end
end

if(stimulatedSession==1)
    
    if (plotON==1)
        if currPlot==1
            close(fff)
        end
    end
    
    DATAstim.dPrime =  DATAstim.dPrime(~any(isnan(DATAstim.dPrime),2),:);
    DATAstim.performance = DATAstim.performance(~any(isnan(DATAstim.performance),2),:);
    
    %plot individual Angles Avg
    noAngles = sum(~structfun(@isempty,DATAstim.angles));
    plotRows = 3;
    plotCols = noAngles;
    currPlot=1;
    
    if (plotON==1)
        ff=figure;clf
        set(ff,'Position',positionGraph2);
    else
        figure('Visible','off');clf;
    end
    set(gcf,'name',tempTitle,'numbertitle','off')
    
    for k = 1:sum(~structfun(@isempty,DATAstim.angles))
        clear p h sigPairsPval sigPairs;
        subplot(plotRows,plotCols,currPlot);
        tempDATA = cell2mat(DATAstim.angles.(angleName{k}));
        meandPrime = mean(tempDATA);
        SEMdPrime = std(tempDATA)/sqrt(length(tempDATA));
        [h,p] = ttest(tempDATA(:,1),tempDATA(:,2));
        if h==1
            sigPairs = [1,2];
            sigPairsPval = p;
        end
        
        plot(tempDATA','k-')
        hold on
        errorbar(meandPrime,SEMdPrime,'r-','LineWidth',5)
        %             if h==1
        %                 sigstar(sigPairs,sigPairsPval,0,0.2)
        %             end
        set(gca,'XTick',[1:2]);
        set(gca,'XTickLabel',['Cntr';'Stim']);
        xlimit = [0.5 2.5];
        set(gca, 'Xlim',[xlimit(1) xlimit(2)])
        ylimit = ylim;
        ylim([0 ylimit(2)]);
        set(gca, 'Ylim',[0 ylimit(2)])
        str = ['P = ', num2str(p)];
        ylabel((angleName{k}));
        text(xlimit(1),0.2,str)
        currPlot=currPlot+1;
        
    end
    
    %plot grouped dPrime, lickBias, performance
    if any(DATAstim.dPrime)
        tempNames = fieldnames(DATAstim);
        tempNames(strcmp('angles',tempNames))=[];
        tempNames(strcmp('latencysDiff',tempNames))=[];
        tempNames(strcmp('angleCountSTIM',tempNames))=[];
        tempNames(strcmp('angleCountNoSTIM',tempNames))=[];
        
        for k=1:length(tempNames);
            tempDATA = DATAstim.(tempNames{k});
            meandPrime = mean(tempDATA);
            SEMdPrime = std(tempDATA)/sqrt(length(tempDATA));
            [~,p] = ttest(tempDATA(:,1),tempDATA(:,2));
            [h,p] = ttest(tempDATA(:,1),tempDATA(:,2));
            %             if h==1
            %                 sigPairs = [1,2];
            %                 sigPairsPval = p;
            %             end
            %
            subplot(plotRows,plotCols,[plotCols+k,plotCols*2+k]);
            
            plot(tempDATA','k-')
            hold on
            errorbar(meandPrime,SEMdPrime,'r-','LineWidth',5)
            
            %             b1 =  bar(meandPrime)
            %             set(b1,'FaceColor','none','EdgeColor','none','LineWidth',2,'BarWidth',0.5)
            %             hold on
            %             errorbar(meandPrime,SEMdPrime,'.','LineWidth',2,'Color','k')
            %             plot(sigPairs,tempDATA','-','LineWidth',2,'Color',[0.4,0.4,0.4])
            %             if h==1
            %                 sigstar(sigPairs,sigPairsPval,0,0.2)
            %             end
            set(gca,'XTick',[1:2]);
            set(gca,'XTickLabel',['Cntr';'Stim']);
            xlimit = [0.5 2.5];
            set(gca, 'Xlim',[xlimit(1) xlimit(2)])
            ylimit = ylim;
            if strcmp('lickBias',tempNames{k}) || strcmp('performance',tempNames{k});
                ylimit = [0 1];
            else
                ylimit = ([0 ylimit(2)]);
            end
            set(gca, 'Ylim',[0 ylimit(2)])
            str = ['P = ', num2str(p)];
            ylabel((tempNames{k}));
            text(xlimit(1),0.2,str)
            
        end
    end
    
    subplot(plotRows,plotCols,[plotCols+k+1,plotCols*2+k+1]);
    meanLatencies = nanmean(DATAstim.latencysDiff');
    stdLatencies = nanstd(DATAstim.latencysDiff');
    plot(tempangles,meanLatencies,'o-k','MarkerSize',3,'LineWidth',3)
    hold on
%     errorbar(tempangles,meanLatencies,stdLatencies,'o-k','MarkerSize',3,'LineWidth',1)
    NumTicks = length(tempangles);
    L = [tempangles(1) tempangles(end)];
    line(L,[0 0],'LineStyle','--','Color','r')
    set(gca,'xTick',tempangles,'XTickLabel',labels,'XAxisLocation','top')
    view(-90,90)
    xlabel('LickLatency');
    ylabel('diff Lick Latency (Stim-Cntr)');
    baseFileName = strcat(baseFileName,'average');
    
    saveas(gca,fullfile('C:\Users\adesniklab\Documents\BehaviorRawData\currFigs\optogenetics',baseFileName),'jpeg');
end
end


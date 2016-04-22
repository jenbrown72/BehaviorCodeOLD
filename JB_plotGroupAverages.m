function [analysedDATA] = JB_plotGroupAverages(AllDATA,normdata,exclude, trimSession, plotON)
%UNTITLED5 Summary of this function goes here
%   AllDATA matrix is generated from : [AllDATA] = JB_groupAverages(AllDATA,listToAnalyse,condition)
%   norm = 1; normalise data norm = 0; raw data
%   exclude = 1 - exclude mice depending on threshold set below
% trimSession=1 (first timetrimmed), =2 (secondtimeTrimmed)
%   plotON =1; plot, plotON = 0, no plot

%load('AllDATA.mat');

positionGraph2 = [1321 548 587 448];
sessionTypeNo = 5; %e.g Full, C1C2C3C4, C2C3, C2, None
trialPerSessionNo = 2; % needs to be set as want the same number of data points per session type.

%exclude data from mice who could perform the task without their whiskers

dataToAnalyse = {'performance';'dPrime'};
thresholdToEx = [0.7; 1]; %performance
% dataToAnalyse = {'anglesdPrime90','anglesdPrime58','anglesdPrime30','anglesdPrime14'};
figureNo=1;

for hh = 1:length(dataToAnalyse)
    plotRows = 2;
    plotCols = 3;
    currPlot = 1;
    
    if (plotON==1)
        ffff=figure(figureNo);clf
        set(ffff,'Position',positionGraph2);
    else
        figure('Visible','off');clf;
    end
    
    tempTitle = [dataToAnalyse{hh},' ', 'Norm', num2str(normdata)];
    set(gcf,'name',tempTitle,'numbertitle','off')
    
    for j = 1:length(AllDATA)
        analysedDATA{j}.parameter{hh}.name = dataToAnalyse{hh};
        analysedDATA{j}.parameter{hh}.tempDATA=nan(5,length(AllDATA{j}.data)*2);
        analysedDATA{j}.parameter{hh}.tempWhiskers = cell(5,length(AllDATA{j}.data));
        k=1;
        kk=1;
        
        for i=1:length(AllDATA{j}.data)
                        analysedDATA{j}.parameter{hh}.tempDATA(1:sessionTypeNo,k:k+1) =  AllDATA{j}.data{i}.(dataToAnalyse{hh})(1:sessionTypeNo,1:trialPerSessionNo);
            analysedDATA{j}.parameter{hh}.tempDATAmean(1:sessionTypeNo,kk) =  mean(AllDATA{j}.data{i}.(dataToAnalyse{hh})(1:sessionTypeNo,1:trialPerSessionNo),2);

       %     analysedDATA{j}.parameter{hh}.tempDATA(1:sessionTypeNo,k:k+1) =  AllDATA{j}.data{i}.(dataToAnalyse{hh}){trimSession}(1:sessionTypeNo,1:trialPerSessionNo);
         %   analysedDATA{j}.parameter{hh}.tempDATAmean(1:sessionTypeNo,kk) =  mean(AllDATA{j}.data{i}.(dataToAnalyse{hh}){trimSession}(1:sessionTypeNo,1:trialPerSessionNo),2);
            tempWhiskers(:,kk) =   AllDATA{j}.data{i}.whiskerID(1:5);
            k= k+2;
            kk = kk+1;
        end
       
        tempDATA = analysedDATA{j}.parameter{hh}.tempDATA;
        tempDATAmean = analysedDATA{j}.parameter{hh}.tempDATAmean;
        
        if (exclude==1)
            idx = tempDATA(find(strcmp('None',tempWhiskers(:,1))),:)>thresholdToEx(hh);
            tempDATA(:,idx) = [];
            tempDATAmean(:,idx) = [];
        end
        
         
        %normalised data
        if normdata==1;
            tempDATA = bsxfun(@rdivide,tempDATA, max(tempDATA));
            tempDATAmean = bsxfun(@rdivide,tempDATAmean, max(tempDATAmean));
        end
        
        k=1;
        clear h p sigPairs sigPairsPval
        
        for kk = 1:size(tempDATA,1)
            [h(kk),p(kk)] = ttest(tempDATA(1,:),tempDATA(kk,:));
            if h(kk)==1
                sigPairs{k} = [1,kk];
                sigPairsPval(k) = [p(kk)];
                k = k+1;
            end
        end
        
        %Plot day one in each condition
        subplot(plotRows,plotCols,currPlot);
        plot(tempDATA(:,1:2:(size(tempDATA,2))),'-ko', 'MarkerSize',3)
        hold on
        set(gca,'XTick',[1:size(tempDATA,1)]);
        set(gca,'XTickLabel',tempWhiskers(1:5));
        ylabelStr = [(dataToAnalyse{hh})];
        ylabel(ylabelStr);
        ylimit = ylim;
        xlimit = xlim;
        text(xlimit(1),(ylimit(2)+0.2),'First Day','HorizontalAlignment','left')
        currPlot=currPlot+1;
        
        %plot average in each condition
        subplot(plotRows,plotCols,currPlot);
        plot(tempDATAmean,'-ko', 'MarkerSize',3)
        hold on
        set(gca,'XTick',[1:size(tempDATA,1)]);
        set(gca,'XTickLabel',tempWhiskers(1:5));
        ylabelStr = [(dataToAnalyse{hh})];
        ylabel(ylabelStr);
        ylimit = ylim;
        xlimit = xlim;
        text(xlimit(1),(ylimit(2)+0.2),'Average over Session Types','HorizontalAlignment','left')
        currPlot=currPlot+1;
        
        %plot mean and stats
        subplot(plotRows,plotCols,currPlot);
        bar(mean(tempDATA,2),0.5,'faceColor',[j-1 j-1 j-1]);
        hold on
        errorbar(mean(tempDATA,2),std(tempDATA,0,2),'.','Color',[0 0 0]);
        sigstar(sigPairs,sigPairsPval);
        set(gca,'XTick',[1:size(tempDATA,1)]);
        set(gca,'XTickLabel',tempWhiskers(1:5));
        ylabelStr = ['average',' ',(dataToAnalyse{hh})];
        ylabel(ylabelStr);
        currPlot=currPlot+1;
    end
    figureNo = figureNo+1;
end

%JB_plotGroupStats(analysedDATA,1)

end

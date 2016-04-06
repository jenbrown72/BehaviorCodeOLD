function [analysedDATA] = JB_plotGroupAverages(AllDATA,plotON)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

positionGraph2 = [1321 548 587 448];

dataToAnalyse = {'Performance';'dPrime'};
figureNo=1;

for hh = 1:length(dataToAnalyse)
    plotRows = 2;
    plotCols = 2;
    currPlot = 1;
    
    if (plotON==1)
        ffff=figure(figureNo);clf
        set(ffff,'Position',positionGraph2);
    else
        figure('Visible','off');clf;
    end
    
    for j = 1:length(AllDATA)
        analysedDATA{j}.parameter{hh}.name = dataToAnalyse{hh};
        analysedDATA{j}.parameter{hh}.tempDATA=nan(5,length(AllDATA{j}.data)*2);
        analysedDATA{j}.parameter{hh}.tempWhiskers = cell(5,length(AllDATA{j}.data));
        k=1;
        kk=1;
        
        for i=1:length(AllDATA{j}.data)
            analysedDATA{j}.parameter{hh}.tempDATA(:,k:k+1) =   AllDATA{j}.data{i}.(dataToAnalyse{hh})(1:5,:);
            tempWhiskers(:,kk) =   AllDATA{j}.data{i}.whiskerTrimPairType(1:5);
            k= k+2;
            kk = kk+1;
        end
        
        tempDATA = analysedDATA{j}.parameter{hh}.tempDATA;
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
        
        % tempTitle = basicPropertiesToPlot{1,1}.mouseID;
        % tempTitle(findstr(tempTitle,'_'))=[];
        % tempTitle = strcat(tempTitle,'_PerformanceTrimming');
        % set(gcf,'name',tempTitle,'numbertitle','off')
        
        subplot(plotRows,plotCols,currPlot);
        plot(tempDATA,'-k')
        hold on
        plot(mean(tempDATA,2),'-r', 'LineWidth',3);
        hold on
        errorbar(mean(tempDATA,2),std(tempDATA,0,2),'r','LineWidth',3);
        set(gca,'XTick',[1:size(tempDATA,1)]);
        %   ylim([0,1])
        set(gca,'XTickLabel',tempWhiskers(1:5));
        ylabel('average Performance')
        currPlot=currPlot+1;
                
        subplot(plotRows,plotCols,currPlot);
        plot(tempDATA,'-k')
        hold on
        bar(mean(tempDATA,2),0.5);
        hold on
        errorbar(mean(tempDATA,2),std(tempDATA,0,2),'.');
        sigstar(sigPairs,sigPairsPval);
        set(gca,'XTick',[1:size(tempDATA,1)]);
        set(gca,'XTickLabel',tempWhiskers(1:5));
        ylabel('average Performance')
        currPlot=currPlot+1;
    end
    figureNo = figureNo+1;
end

%JB_plotGroupStats(analysedDATA,1)

end

function [analysedDATA] = JB_plotGroupAveragesSTIM(AllDATA,normdata,exclude, plotON)
%UNTITLED5 Summary of this function goes here
%   AllDATA matrix is generated from : [AllDATA] = JB_groupAverages(AllDATA,listToAnalyse,condition)
%   norm = 1; normalise data norm = 0; raw data
%   exclude = 1 - exclude mice depending on threshold set below
% trimSession=1 (first timetrimmed), =2 (secondtimeTrimmed)
%   plotON =1; plot, plotON = 0, no plot

%load('AllDATA.mat');

positionGraph2 = [1321 548 587 448];
sessionTypeNo = 2; %e.g Cntr Stim
labels = ['Cntr'; 'Stim']

trialPerSessionNo = 2; % needs to be set as want the same number of data points per session type.

%exclude data from mice who could perform the task without their whiskers

dataToAnalyse = {'performance';'dPrime';'lickBias'};
thresholdToEx = [0.7; 1]; %performance
% dataToAnalyse = {'anglesdPrime90','anglesdPrime58','anglesdPrime30','anglesdPrime14'};
figureNo=1;

%find average of each parameter for control and stim

    plotRows = 2;
    plotCols = 3;
    currPlot = 1;
    
    if (plotON==1)
        ffff=figure(figureNo);clf
        set(ffff,'Position',positionGraph2);
    else
        figure('Visible','off');clf;
    end


for hh = 1:length(dataToAnalyse)

    
    for j = 1:length(AllDATA{1}.data)
        tempDATA(j,:) = nanmean(AllDATA{1}.data{j}.(dataToAnalyse{hh}));
        tempDATAmean(j,:) = nanstd(AllDATA{1}.data{j}.(dataToAnalyse{hh}));
    end
    
    k=1;
    clear h p sigPairs sigPairsPval
    
    [h,p] = ttest(tempDATA(1,:),tempDATA(2,:));
    if h==1
        sigPairs = [1,2];
        sigPairsPval = [p];
    end


%Plot day one in each condition
subplot(plotRows,plotCols,currPlot);
plot([1,2],tempDATA','-ko', 'MarkerSize',3)
hold on
plot([1,2],nanmean(tempDATA)')
set(gca,'XTick',[1:2]);
xlimit = [0.5 2.5];
set(gca,'XTickLabel',labels,'xlim',xlimit);
ylabelStr = [(dataToAnalyse{hh})];
ylabel(ylabelStr);
% ylimit = ylim;
% 
% text(xlimit(1),(ylimit(2)+0.2),'First Day','HorizontalAlignment','left')
currPlot=currPlot+1;
end

%JB_plotGroupStats(analysedDATA,1)

end

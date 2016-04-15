function [] = JB_plotLickLatency(basicPropertiesToPlot,possibleAngles,plotON)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



for h = 1:length(numPoints)
    
    activeAngles = unique(basicPropertiesToPlot{h,1}.frequencyLicks(:,2));
    activeAngles(~any(activeAngles,2),:)=[];
    
    noSubPlots=length(activeAngles);
    currPlot = 1;
    
    
    for kk=1:length(activeAngles)
        
        idxAngle = find(basicPropertiesToPlot{h,1}.frequencyLicks(:,2) == activeAngles(kk));
        tempLickFreq = basicPropertiesToPlot{h,1}.frequencyLicks(idxAngle,1);
        tempLickLatency = basicPropertiesToPlot{h,1}.firstLickLatency(idxAngle,1);
        
        subplot(2,noSubPlots,currPlot)
        hist(tempLickLatency)
        hold on
        
        subplot(2,noSubPlots,currPlot+length(activeAngles))
        hist(tempLickFreq)
        hold on
        currPlot = currPlot+1;
    end
    
end


function [] = JB_basicBehaviorProperties(currDirectory)

load('DATA.mat')
sessionType = 1;
basisPropertiesID{sessionType,1} = 'sessionType';
date = 2;
basisPropertiesID{date,1} = 'date';
datenum = 3;
basisPropertiesID{datenum,1} = 'datenum';
sessionDuration = 4;
basisPropertiesID{sessionDuration,1} = 'sessionDuration';
totalStepsPerMin = 5;
basisPropertiesID{totalStepsPerMin,1} = 'totalStepsPerMin';
performance = 6;
basisPropertiesID{performance,1} = 'performance';
HIT = 7;
basisPropertiesID{HIT,1} = 'performanceTypeHIT';
FA = 8;
basisPropertiesID{FA,1} = 'performanceTypeFA';
MISS = 9;
basisPropertiesID{MISS,1} = 'performanceTypeMISS';
CR = 10;
basisPropertiesID{CR,1} = 'performanceTypeCR';
possibleAngles = [225;241;254;263;266;268;270;272;274;277;284;299;315];

startAngleID = 10;
basisPropertiesID{startAngleID,1} = 'anglePerformanceHIT';
% FAangle = 12;
% basisPropertiesID{FAangle,1} = 'anglePerformanceFA';
% MISSangle = 13;
% basisPropertiesID{MISSangle,1} = 'anglePerformanceMISS';
% CRangle = 14;
% basisPropertiesID{CRangle,1} = 'anglePerformanceCR';

basicProperties = {};
encoder0Pos = 1;
LEDstate = 5;
rawSessionTime = 7;
trialType = 13;
angle = 11;
rewardDuration = 9;

for i=1:length(DATA.allFiles);
    
    metaData = DATA.allFiles{i}.metaData;
    
    %find out if this session used optogenetics
    optoIdx = find(strcmp('Optogenetics = ', metaData));
    tempOpto = metaData{optoIdx,2};
     basicProperties{i,1}.optogenetics = tempOpto;
    
    basicProperties{i,1}.date = DATA.allFiles{i}.date;
    basicProperties{i,1}.datenum = DATA.allFiles{i}.dateFromFile;
    
    tempDATA = DATA.allFiles{i}.rawData;
    basicProperties{i,1}.sessionDuration = ((tempDATA(end,rawSessionTime)-tempDATA(1,rawSessionTime))/1000)/60; % Calculate duration of session in minutes
    threshold = 50000;
    aboveThreshold = find(tempDATA(:,encoder0Pos)>threshold);
    for j = 1:length(aboveThreshold)
        tempDATA(aboveThreshold(j),:)=nan;
    end
    
    totalSteps = nansum(diff(tempDATA(:,1))>0);
    basicProperties{i,1}.totalStepsPerMin = totalSteps/basicProperties{i,1}.sessionDuration; % Calculate running velocity
    
    %Look at raw data to get stats on performance
    findtrialType = diff(DATA.allFiles{i}.rawData(:,trialType));
    [idx ~] = find(findtrialType>0);
    
    temptrialTypes = DATA.allFiles{i}.rawData(idx+2,trialType);
    tempAngleTypes = DATA.allFiles{i}.rawData(idx+2,angle);
    tempRewardDuration = DATA.allFiles{i}.rawData(idx+2,rewardDuration);
    tempOptogenetics = DATA.allFiles{i}.rawData(idx-10,LEDstate);
    % angleUsed = unique(tempAngleTypes);
    % [idx] = find(angleUsed>0); %no 0
    
    for kk=1:length(possibleAngles)
        
        stimCount = 1;
        nostimCount = 1;
        
        idxAngles = find(tempAngleTypes==possibleAngles(kk));
        performanceAngleTemp = temptrialTypes(idxAngles);
        
        if ~isempty(idxAngles)
            if tempOpto==1
                
                for ff = 1:length(idxAngles)
                    
                    if (tempOptogenetics(idxAngles(ff))==1)
                        
                        performanceAngleTempSTIM(stimCount,1) = temptrialTypes(idxAngles(ff));
                        stimCount = stimCount+1;
                        
                    else
                        
                        performanceAngleTempNoSTIM(nostimCount,1) = temptrialTypes(idxAngles(ff));
                        nostimCount = nostimCount+1;
                    end
                    
                    
                end
                
            end
            
        else 
            
            performanceAngleTempSTIM = [];
             performanceAngleTempNoSTIM = [];
            
        end
        
        HITangle = sum(performanceAngleTemp==1);
        FAangle = sum(performanceAngleTemp==2);
        MISSangle = sum(performanceAngleTemp==3);
        CRangle = sum(performanceAngleTemp==4);
        %Performance at each angle
        basicProperties{i,1}.performance{kk} =  (HITangle+CRangle)/(length(performanceAngleTemp));
        %Probability of licking
        basicProperties{i,1}.probLicking{kk} =  (HITangle+FAangle)/(MISSangle+CRangle+HITangle+FAangle);
        
        if (tempOpto==1)
            HITangleSTIM = sum(performanceAngleTempSTIM==1);
            FAangleSTIM = sum(performanceAngleTempSTIM==2);
            MISSangleSTIM = sum(performanceAngleTempSTIM==3);
            CRangleSTIM = sum(performanceAngleTempSTIM==4);
            basicProperties{i,1}.performanceSTIM{kk} =  (HITangleSTIM+CRangleSTIM)/(length(performanceAngleTempSTIM));
            basicProperties{i,1}.probLickingSTIM{kk} =  (HITangleSTIM+FAangleSTIM)/(MISSangleSTIM+CRangleSTIM+HITangleSTIM+FAangleSTIM);
            
            HITanglenoSTIM = sum(performanceAngleTempNoSTIM==1);
            FAanglenoSTIM = sum(performanceAngleTempNoSTIM==2);
            MISSanglenoSTIM = sum(performanceAngleTempNoSTIM==3);
            CRanglenoSTIM = sum(performanceAngleTempNoSTIM==4);
            basicProperties{i,1}.performanceNoSTIM{kk} =  (HITanglenoSTIM+CRanglenoSTIM)/(length(performanceAngleTempNoSTIM));
            basicProperties{i,1}.probLickingNoSTIM{kk} =  (HITanglenoSTIM+FAanglenoSTIM)/(MISSanglenoSTIM+CRanglenoSTIM+HITanglenoSTIM+FAanglenoSTIM);
        end
        
    end
    
    basisPropertiesID{startAngleID,1};
    
    %Calculate how many of each trial type there were
    basicProperties{i,1}.HIT = sum(temptrialTypes==1); %HitTrialCount
    basicProperties{i,1}.FA = sum(temptrialTypes==2); % FATrialCount
    basicProperties{i,1}.MISS = sum(temptrialTypes==3); %MissTrialCount
    basicProperties{i,1}.CR = sum(temptrialTypes==4); %CRTrialCount
    
    basicProperties{i,1}.sessionperformance = (basicProperties{i,1}.HIT+basicProperties{i,1}.CR)/(length(temptrialTypes));
    
    %import as metaDATA cell array
    n=1;
    
    
    
    for j=1:length(metaData)
        
        match(j,1) = strcmp('Orientation Selected = ', metaData{j,1});
        
        if(strcmp('Auto reward = ', metaData{j,1}));
            autoSession = metaData{j,2};
        end
        
        
    end
    
    [ind] = cell2mat(metaData(match,2));
    tempAngles = ind(ind>0);
    %         if (match(j,1)==1) && (metaData{j,2}>0)
    %             tempAngles(n,1) = metaData{j,2};
    %             n=n+1;
    %         end
    %
    %         %AutoReward
    %         if(strcmp('Auto reward = ', metaData{j,1}));
    %             autoSession = metaData{j,2};
    %         end
    %
    %     end
    
    %Define session type
    tempsessionType = strcat('S',num2str(length(tempAngles)));
    
    if autoSession==1
        tempsessionType = strcat(tempsessionType,'auto');
    end
    
    basicProperties{i,1}.sessionType = tempsessionType;
    
end

for k = 1:length(basicProperties),
    
    days{k,1} = basicProperties{k,1}.datenum(1:8);
    
end

plotNumber = 1;

[un idx_last idx] = unique(days(:,1));
idx = sort(idx);
unique_idx = accumarray(idx(:),(1:length(idx))',[],@(x) {sort(x)});

for k = 1:length(unique_idx)
    tempRunVal=[];
    if length(unique_idx{k})==1
        basicPropertiesToPlot{plotNumber,1} = basicProperties{unique_idx{k}};
        plotNumber = plotNumber+1;
    elseif length(unique_idx{k})>1
        t = unique_idx{k};
        
        for kk = 1:length(t)
            tempRunVal(kk,1) = basicProperties{t(kk)}.sessionDuration;
        end
        
        [id di] = max(tempRunVal); %find the data set where the mouse ran the most
        basicPropertiesToPlot{plotNumber,1} = basicProperties{t(di)};
        plotNumber = plotNumber+1;
    end
end
%%

figure;clf
%figure('Visible','off');clf;
tempTitle = DATA.mouseID;
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

subplot(noSubPlots,1,2)
for j = 1:length(numPoints);
    plot(numPoints(j),basicPropertiesToPlot{j}.sessionDuration,'or','MarkerSize', 10,'MarkerFaceColor','r')
    hold on
end

ylabel('sessionDuration');
xlabel('Session Number');

%plot primary session type per day
subplot(noSubPlots,1,3)
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

ColorCodeColors = [1.0;0.9;0.8;0.6;0.4;0.2;0];

for j = 1:length(numPoints);
    plot(numPoints(j),colorCode(j,1),'o','MarkerSize', 10, 'MarkerFaceColor',[0.4,ColorCodeColors(colorCode(j)),0.8],'Color', [0.4,ColorCodeColors(colorCode(j)),0.8])
    hold on
end

NumTicks = max(colorCode);
L = get(gca,'YLim');
set(gca,'YTick',linspace(L(1),L(2),NumTicks))
set(gca,'YTickLabel',SessionTypes(1:max(colorCode)))
%set(gca,'YTickLabel',{'S1Auto','S1', 'S2'})
ylabel('Session Type');
xlabel('Session Number');

%plot performance
subplot(noSubPlots,1,4)

for j = 1:length(numPoints);
    plot(numPoints(j),basicPropertiesToPlot{j,1}.sessionperformance,'or','MarkerSize', 10,'MarkerFaceColor',[0.4,ColorCodeColors(colorCode(j)),0.8],'Color', [0.4,ColorCodeColors(colorCode(j)),0.8])
    hold on
end

ylim([0 1])
xlim([0 length(numPoints)])
ylabel('Performance');
xlabel('Session Number');

baseFileName = strcat(tempTitle); %save old figure
saveas(gca,fullfile('C:\Users\adesniklab\Documents\BehaviorRawData\currFigs\basicSessionProperties',baseFileName),'jpeg');

%plot sessionType performance to get an idea of how inclinded to lick the
%mouse is

plotTotal = 4;
numFigs = 1;
%figure('Visible','off');clf;
figure;clf;
set(gcf,'name',tempTitle,'numbertitle','off')
currPlot = 1;

for h = 1:length(numPoints)
    saved = 0;
    activeAngles = cell2mat(basicPropertiesToPlot{h,1}.performance);
    probLick = cell2mat(basicPropertiesToPlot{h,1}.probLicking);
    plotAngles = possibleAngles;
    [~,c] = find(isnan(activeAngles));
    activeAngles(c) = [];
    plotAngles(c) = [];
    probLick(c) = [];
    
    trialTypecombo = [basicPropertiesToPlot{h,1}.HIT basicPropertiesToPlot{h,1}.MISS;basicPropertiesToPlot{h,1}.FA basicPropertiesToPlot{h,1}.CR];
    
    if(length(activeAngles)>2) %if more than 2 angles were presented
        subplot(plotTotal,3,currPlot);
        plot(plotAngles,activeAngles,'o-');
        currPlot=currPlot+1;
        xlabel('Angles');
        ylabel('Performance');
        ylim([0 1])
        
        subplot(plotTotal,3,currPlot);
        plot(plotAngles,probLick,'o-');
        currPlot=currPlot+1;
        xlabel('Angles');
        ylabel('Licking Probability');
        ylim([0 1])
        
        subplot(plotTotal,3,currPlot)
        bar(trialTypecombo, 'stacked');
        hold on
        XlableAxis = {'HIT/MISS'; 'FA/CR'};
        set(gca,'XTickLabel',XlableAxis);
        currPlot=currPlot+1;
        
        if rem(currPlot,plotTotal)==1 %if the value is divisable by 4 - open a new plot
            baseFileName = strcat(tempTitle,num2str(numFigs)); %save old figure
            saveas(gca,fullfile('C:\Users\adesniklab\Documents\BehaviorRawData\currFigs\psychometricCurves',baseFileName),'jpeg');
            saved = 1;
            figure;clf;
            numFigs = numFigs+1;
            currPlot = 1;
        end
        
    end
    
    if saved==0
        
        baseFileName = strcat(tempTitle,num2str(numFigs)); %save old figure
        saveas(gca,fullfile('C:\Users\adesniklab\Documents\BehaviorRawData\currFigs\psychometricCurves',baseFileName),'jpeg');
        
    end
    
    
    
end


plotTotal = 4;
numFigs = 1;
%figure('Visible','off');clf;
figure;clf;
set(gcf,'name',tempTitle,'numbertitle','off')
currPlot = 1;

for h = 1:length(numPoints)
    
    if (basicPropertiesToPlot{h,1}.optogenetics)
        
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

  subplot(plotTotal,2,currPlot);
        plot(plotAngles,activeAnglesSTIM,'o-r');
        hold on;
         plot(plotAngles,activeAnglesnoSTIM,'o-k');
        currPlot=currPlot+1;
        xlabel('Angles');
        ylabel('Performance');
        ylim([0 1])
        
        subplot(plotTotal,2,currPlot);
        plot(plotAngles,probLickSTIM,'o-r');
        hold on;
          plot(plotAngles,probLickNoSTIM,'o-k');
        currPlot=currPlot+1;
        xlabel('Angles');
        ylabel('Licking Probability');
        ylim([0 1])
        
    end
    
            baseFileName = strcat(tempTitle); %save old figure
        saveas(gca,fullfile('C:\Users\adesniklab\Documents\BehaviorRawData\currFigs\optogenetics',baseFileName),'jpeg');
    
end

        
end
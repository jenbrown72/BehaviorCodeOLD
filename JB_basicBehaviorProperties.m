function [] = JB_basicBehaviorProperties

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
rawSessionTime = 7;
trialType = 13;
angle = 11;
rewardDuration = 9;

for i=1:length(DATA.allFiles);
    
    basicProperties{date,i} = DATA.allFiles{i}.date;
    basicProperties{datenum,i} = DATA.allFiles{i}.dateFromFile;
    
    tempDATA = DATA.allFiles{i}.rawData;
    basicProperties{sessionDuration,i} = ((tempDATA(end,rawSessionTime)-tempDATA(1,rawSessionTime))/1000)/60; % Calculate duration of session in minutes
    threshold = 50000;
    aboveThreshold = find(tempDATA(:,encoder0Pos)>threshold);
    for j = 1:length(aboveThreshold)
        tempDATA(aboveThreshold(j),:)=nan;
    end
    
    totalSteps = nansum(diff(tempDATA(:,1))>0);
    basicProperties{totalStepsPerMin,i} = totalSteps/basicProperties{sessionDuration,i}; % Calculate running velocity
    
    %Look at raw data to get stats on performance
    findtrialType = diff(DATA.allFiles{i}.rawData(:,trialType));
    [idx idx2] = find(findtrialType>0);
    
    temptrialTypes = DATA.allFiles{i}.rawData(idx+2,trialType);
    tempAngleTypes = DATA.allFiles{i}.rawData(idx+2,angle);
    tempRewardDuration = DATA.allFiles{i}.rawData(idx+2,rewardDuration);
    % angleUsed = unique(tempAngleTypes);
    % [idx] = find(angleUsed>0); %no 0
    
    for kk=1:length(possibleAngles)
        idxAngles = find(tempAngleTypes==possibleAngles(kk));
        performanceAngleTemp = temptrialTypes(idxAngles);
        
        HITangle = sum(performanceAngleTemp==1);
        FAangle = sum(performanceAngleTemp==2);
        MISSangle = sum(performanceAngleTemp==3);
        CRangle = sum(performanceAngleTemp==4);
        basicProperties{startAngleID+kk,i} =  (HITangle+CRangle)/(length(performanceAngleTemp));
        basicProperties{startAngleID+length(possibleAngles)+kk,i} =  (HITangle+FAangle)/(MISSangle+CRangle+HITangle+FAangle);

    end
    
    basisPropertiesID{startAngleID,1};
    
    %Calculate how many of each trial type there were
    basicProperties{HIT,i} = sum(temptrialTypes==1); %HitTrialCount
    basicProperties{FA,i} = sum(temptrialTypes==2); % FATrialCount
    basicProperties{MISS,i} = sum(temptrialTypes==3); %MissTrialCount
    basicProperties{CR,i} = sum(temptrialTypes==4); %CRTrialCount
    
    basicProperties{performance,i} = (basicProperties{HIT,i}+basicProperties{CR,i})/(sum(temptrialTypes));

    %import as metaDATA cell array
    n=1;
    
    metaData = DATA.allFiles{i}.metaData;
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
    
    basicProperties{sessionType,i} = tempsessionType;
    
end

for k = 1:size((basicProperties),2)
    
    days{k} = basicProperties{3,k}(1:8);
    
end

plotNumber = 1;

[un idx_last idx] = unique(days(1,:));
idx = sort(idx);
unique_idx = accumarray(idx(:),(1:length(idx))',[],@(x) {sort(x)});

for k = 1:length(unique_idx)
    tempRunVal=[];
    if length(unique_idx{k})==1   
        basicPropertiesToPlot(:,plotNumber) = basicProperties(:,(unique_idx{k}));    
        plotNumber = plotNumber+1;    
    elseif length(unique_idx{k})>1
        t = unique_idx{k};
        
        for kk = 1:length(t)
            tempRunVal(kk,1) = basicProperties{4,t(kk)};   
        end
        
        [id di] = max(tempRunVal); %find the data set where the mouse ran the most
        basicPropertiesToPlot(:,plotNumber) = basicProperties(:,t(di));
        plotNumber = plotNumber+1; 
    end
end
%%

figure(1);clf
noSubPlots = 4;
subplot(noSubPlots,1,1)


numPoints = 1:1:size((basicPropertiesToPlot), 2);
for j = 1:length(numPoints);
    subplot(noSubPlots,1,1)
    plot(numPoints(j),basicPropertiesToPlot{5,j},'or','MarkerSize', 10,'MarkerFaceColor','r')
    hold on
end

ylabel('totalStepsPerMin');
xlabel('Session Number');

subplot(noSubPlots,1,2)
for j = 1:length(numPoints);
    plot(numPoints(j),basicPropertiesToPlot{4,j},'or','MarkerSize', 10,'MarkerFaceColor','r')
    hold on
end

ylabel('sessionDuration');
xlabel('Session Number');

%plot primary session type per day
subplot(noSubPlots,1,3)
SessionTypes = {'S1auto' ; 'S1'; 'S2'; 'S6'; 'S8'; 'S10'; 'S12'};

for j = 1:length(numPoints);
    if strcmp('S1auto', basicPropertiesToPlot{1,j})
        colorCode(j,1) = 1;
    elseif strcmp('S1',basicPropertiesToPlot{1,j})
        colorCode(j,1) = 2;
    elseif strcmp('S2',basicPropertiesToPlot{1,j})
        colorCode(j,1) = 3;
    elseif strcmp('S6',basicPropertiesToPlot{1,j})
        colorCode(j,1) = 4;
    elseif strcmp('S8',basicPropertiesToPlot{1,j})
        colorCode(j,1) = 5;
           elseif strcmp('S10',basicPropertiesToPlot{1,j})
        colorCode(j,1) = 6;
           elseif strcmp('S12',basicPropertiesToPlot{1,j})
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
    plot(numPoints(j),basicPropertiesToPlot{6,j},'or','MarkerSize', 10,'MarkerFaceColor',[0.4,ColorCodeColors(colorCode(j)),0.8],'Color', [0.4,ColorCodeColors(colorCode(j)),0.8])
    hold on
end

ylim([0 1])
xlim([0 length(numPoints)])
ylabel('Performance');
xlabel('Session Number');

%plot sessionType performance to get an idea of how inclinded to lick the
%mouse is


plotTotal = 0;

for h = 1:size((basicPropertiesToPlot), 2)

    activeAngles = [basicPropertiesToPlot{11:23,h}];
    
if(sum(~isnan(activeAngles),2)>1) %if more than 2 angles were presented
    plotTotal=plotTotal+1; 
end

end


figure(2); clf;
currPlot = 1;

for h = 1:size((basicPropertiesToPlot), 2)

    activeAngles = [basicPropertiesToPlot{11:23,h}];
    probLick = [basicPropertiesToPlot{24:36,h}];
    plotAngles = possibleAngles;
    [r,c] = find(isnan(activeAngles));
    activeAngles(c) = [];
    plotAngles(c) = [];
    probLick(c) = [];
    
    trialTypecombo = [basicPropertiesToPlot{HIT,h} basicPropertiesToPlot{MISS,h};basicPropertiesToPlot{FA,h} basicPropertiesToPlot{CR,h}];

    
if(sum(~isnan(activeAngles),2)>1) %if more than 2 angles were presented
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
    
    
end

end
end
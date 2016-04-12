function [basicPropertiesToPlot possibleAngles] = JB_basicBehaviorPropertiesNeat(plotON,checkIndividualPlots,cleanUp)

%plotON =1, plot fig/ =0 no plot
%cleanUp = 1: stop analysis of session early if critiria is reached:

if nargin <3;
    cleanUp=0;
end

if nargin <2;
    checkIndividualPlots=0;
end


% analysis stop early criteria
stopExThreshold = 32; % no of consequtive misses and correct rejections to say the mouse wasn't engaged and stop analysis
minNoTrials = 0;

load('DATA.mat')
possibleAngles = [225;241;254;263;266;268;270;272;274;277;284;299;315];
basicProperties = cell(100,1);
%define which txt file column the different digital inputs are read in from
%processing

positionGraph1 = [1555 293 343 702];

encoder0Pos = 1;
LEDstate = 5;
count = 6;
rawSessionTime = 7;
trialType = 13;
angle = 11;
rewardDuration = 9;

for i=1:length(DATA.allFiles);
    
    %Extract basic informaiton from the data file
    basicProperties{i,1}.date = DATA.allFiles{i}.date;
    basicProperties{i,1}.datenum = DATA.allFiles{i}.dateFromFile;
    basicProperties{i,1}.sessionDuration = ((DATA.allFiles{i}.rawData(end,rawSessionTime)-DATA.allFiles{i}.rawData(1,rawSessionTime))/1000)/60; % Calculate duration of session in minutes
    basicProperties{i,1}.fileLength = length(DATA.allFiles{i}.rawData);
    basicProperties{i,1}.trimType = NaN;
    
    tempName = DATA.allFiles{i}.name(findstr('_201',DATA.allFiles{i}.name):findstr('_Box',DATA.allFiles{i}.name));
    tempLocation = findstr(tempName,'_');
    for hh = 1:length(tempLocation)-1;
        if(tempLocation(hh+1)-tempLocation(hh))==2;
            tempName = strcat( tempName(1:tempLocation(hh)),'0',tempName(tempLocation(hh)+1:end));
            tempLocation = tempLocation+1;
        end
    end
    
    toDelete = strfind(tempName,'_');
    tempName(toDelete)='.';
    basicProperties{i,1}.namedata = tempName;
    
    %get running velocity
    threshold = 50000;
    aboveThreshold = find(DATA.allFiles{i}.rawData(:,encoder0Pos)>threshold);
    for j = 1:length(aboveThreshold)
        DATA.allFiles{i}.rawData(aboveThreshold(j),:)=nan;
    end
    
    basicProperties{i,1}.totalStepsPerMin = nansum(diff(DATA.allFiles{i}.rawData(:,1))>0)/basicProperties{i,1}.sessionDuration; % Calculate running velocity
    
    %Extract information from the metaData file
    metaData = DATA.allFiles{i}.metaData;
    %find out if this session included a negative stim
    puffIdx = find(strcmp('Negative Reinforcer = ', metaData));
    puffUsed = metaData{puffIdx,2};
    if ischar(puffUsed)==1
        basicProperties{i,1}.negReinforcer = str2num(puffUsed);
    else
        basicProperties{i,1}.negReinforcer = puffUsed;
    end
    
    %find out if this session included a water schedule
    waterSchIdx = find(strcmp('Water Schedule = ', metaData));
    waterSchUsed = metaData{waterSchIdx,2};
    if ischar(waterSchUsed)==1
        basicProperties{i,1}.waterSchedule = str2num(waterSchUsed);
    else
        basicProperties{i,1}.waterSchedule = waterSchUsed;
    end
    
    %find out if this session used optogenetics
    optoIdx = find(strcmp('Optogenetics = ', metaData));
    tempOpto = metaData{optoIdx,2};
    if ischar(tempOpto)==1
        basicProperties{i,1}.optogenetics = str2num(tempOpto);
    else
        basicProperties{i,1}.optogenetics = tempOpto;
    end
    
    %find out if this session was a trimming session
    trimIdx = find(strcmp('Trimming = ', metaData));
    if ~isempty(trimIdx)
        tempTrim = metaData{trimIdx,2};
        basicProperties{i,1}.trimming = tempTrim;
        trimPattern = {basicProperties{i,1}.trimming};
        
        countTrim = strfind(trimPattern,'C');
        if length(countTrim{1})==4; %this is a row data)
            basicProperties{i,1}.trimType = 'Row';
        end
    else
        basicProperties{i,1}.trimming = [];
    end
    
    %find out if this session had segments
    segIdx = find(strcmp('Segments = ', metaData));
    if ~isempty(segIdx)
        tempSeg = metaData{segIdx,2};
        tempSegNum = str2num(tempSeg);
        if (tempSegNum(1)==1) && (tempSegNum(2)==1)
            tempData = DATA.allFiles{i}.rawData;
            index = find(all(tempData==1,2));
            index = [index; length(tempData)]; %add a one to find the start of the first segment
            segments(1,1) = 1; %add first trial
            for kk=1:length(index)-1;
                countCummulative = tempData(index(kk)-1,6);
                tempData(index(kk)+1:index(kk+1),6) = tempData(index(kk)+1:index(kk+1),6)+countCummulative;
                basicProperties{i,1}.segments(kk,2) = countCummulative;
                basicProperties{i,1}.segments(kk+1,1) = countCummulative+1;
            end
            basicProperties{i,1}.segments(length(basicProperties{i,1}.segments),2) = tempData(length(tempData),6); %add last trial
        else
            basicProperties{i,1}.segments = str2num(tempSeg);
        end
    else
        basicProperties{i,1}.segments = [];
    end
    
    %Look at raw data to get stats on performance
    [idx, ~] = find(diff(DATA.allFiles{i}.rawData(:,trialType))>0);
    temptrialTypes = DATA.allFiles{i}.rawData(idx+2,trialType);
    
    %plot the session performance
    if ~isempty(temptrialTypes)
        if (exist('checkIndividualPlots','var') && checkIndividualPlots==1)
            f = figure;clf
            set(f,'Position',positionGraph1);
        else
            figure('Visible','off');clf;
        end
        
        for kk=1:length(temptrialTypes)
            if temptrialTypes(kk)==1; % Hit
                line([0 1], [kk kk], 'Color', 'g','LineWidth',2);
                hold on;
            elseif temptrialTypes(kk)==3; % No Lick on GO stimulus = Miss
                line([0 1], [kk kk], 'Color', 'r','LineWidth',2);
                hold on;
            elseif temptrialTypes(kk)==4; % No Lick on Distractor stimulus = Correct Rejection
                line([2 3], [kk kk], 'Color', 'g','LineWidth',2);
                hold on;
            elseif temptrialTypes(kk)==2; % Lick on Distractor stimulus = False Alarm
                line([2 3], [kk kk], 'Color', 'r','LineWidth',2);
                hold on;
            end
        end
        
        xlim([-0.5 3.5])
     %   ylim([1 length(temptrialTypes)])
        set(gca,'YDir','reverse');
        
        %find how many consequtive trials with no lick there were
        ff = (temptrialTypes==3) | (temptrialTypes==4); %miss and correct rejection
        noresponseCount = diff(find(~(ff==1)));
        tempCum = 1;
        idxStop = [];
        if ~isempty(ff)
            cumCount = nan(length(ff),1);
            for v = 1:length(ff)-1
                if (ff(v) && ff(v+1)) ==1;
                    tempCum = tempCum+1;
                else
                    tempCum = 1;
                end
                cumCount(v,1) = tempCum;
            end
            idxStop = (find(cumCount>stopExThreshold));
            if ~isempty(idxStop) && (cleanUp == 1)
                stopIdx = idxStop(1)-stopExThreshold;
                if stopIdx>minNoTrials
                    X = ['Date ', basicProperties{i,1}.date(1:11), ' Stopping Session Analysis Early By', num2str(length(idx)- stopIdx), 'Trials', 'On trial No.', num2str(stopIdx),'Out Of',num2str(length(idx))] ;
                    line([xlim], [stopIdx stopIdx], 'Color', 'k','LineWidth',2);
                    idx(stopIdx:end)=[];
                    disp(X)
                end
            end
        end
        if (checkIndividualPlots==1)
        waitforbuttonpress;
        end
    end
    
    temptrialTypes = DATA.allFiles{i}.rawData(idx+2,trialType);
    tempCount = DATA.allFiles{i}.rawData(idx,count);
    tempAngleTypes = DATA.allFiles{i}.rawData(idx+2,angle);
    tempRewardDuration = DATA.allFiles{i}.rawData(idx+2,rewardDuration);
    tempOptogenetics = DATA.allFiles{i}.rawData(idx-10,LEDstate);
    
    %if this was a segmented session, calculate performance in each session
    if ~isempty( basicProperties{i,1}.segments);
        sessionDivisions = basicProperties{i,1}.segments;
        for j = 1:length(sessionDivisions)
            if sessionDivisions(j,1)==0;
                idxStart = 1;
            else
                idxStart = sessionDivisions(j,1);
            end
            idxEnd = sessionDivisions(j,2);
            tempSegmentData = temptrialTypes(idxStart:idxEnd,:);
            [Hit, FA, Miss,CR, performance, licking, dprime] = JB_countTrialTypes(tempSegmentData);
            basicProperties{i,1}.segmentPerformance(j,1) =  performance;
            basicProperties{i,1}.segmentdPrime(j,1) = dprime;
        end
    else
        basicProperties{i,1}.segmentPerformance = [];
        basicProperties{i,1}.segmentdPrime = [];
    end
    
    %Calculare performance for each angle presented and overall
    %performance
    
    for kk=1:length(possibleAngles)
        stimCount = 1;
        nostimCount = 1;
        idxAngles = find(tempAngleTypes==possibleAngles(kk));
        performanceAngleTemp = temptrialTypes(idxAngles);
        [Hit, FA, Miss,CR, performance, licking, dprime] = JB_countTrialTypes(performanceAngleTemp);
        basicProperties{i,1}.performance{kk} =  performance;
        basicProperties{i,1}.probLicking{kk} =  licking;
        basicProperties{i,1}.trialTypesAngle{:,kk} = performanceAngleTemp;
        basicProperties{i,1}.performanceSTIM{kk} = nan;
        basicProperties{i,1}.probLickingSTIM{kk} = nan;
        basicProperties{i,1}.trialTypesAngleSTIM{kk} = nan;
        basicProperties{i,1}.performanceNoSTIM{kk} = nan;
        basicProperties{i,1}.probLickingNoSTIM{kk} = nan;
        basicProperties{i,1}.trialTypesAngleNoSTIM{kk} = nan;
        
        if ~isempty(idxAngles) %if this angle was used
            performanceAngleTempSTIM = nan;
            performanceAngleTempNoSTIM = nan;
            if basicProperties{i,1}.optogenetics==1 %if optogenetics was used calculate performance with and without stim
                for ff = 1:length(idxAngles)
                    if (tempOptogenetics(idxAngles(ff))==1)
                        performanceAngleTempSTIM(stimCount,1) = temptrialTypes(idxAngles(ff));
                        stimCount = stimCount+1;
                    else
                        performanceAngleTempNoSTIM(nostimCount,1) = temptrialTypes(idxAngles(ff));
                        nostimCount = nostimCount+1;
                    end
                end
                
                [Hit, FA, Miss,CR, performance, licking, dprime] = JB_countTrialTypes(performanceAngleTempSTIM);
                basicProperties{i,1}.trialTypesAngleSTIM{:,kk} = performanceAngleTempSTIM;
                basicProperties{i,1}.performanceSTIM{kk} = performance;
                basicProperties{i,1}.probLickingSTIM{kk} =  licking;
                
                [Hit, FA, Miss,CR, performance, licking, dprime] = JB_countTrialTypes(performanceAngleTempNoSTIM);
                basicProperties{i,1}.trialTypesAngleNoSTIM{:,kk} = performanceAngleTempNoSTIM;
                basicProperties{i,1}.performanceNoSTIM{kk} =  performance;
                basicProperties{i,1}.probLickingNoSTIM{kk} =  licking;
            end
        end
    end
    
    %Find angle pairs
    
    for d = 1:(floor(length(possibleAngles)/2));
        pair1 = basicProperties{i,1}.trialTypesAngle{1,(1+d-1)};
        pair2 = basicProperties{i,1}.trialTypesAngle{1,end-(d-1)};
        if (~isempty(pair1) || ~isempty(pair2))
            %find which was the go and nogo
            if any(pair1==1) || any(pair1==3) %pair 1 go stim
                [Hit, ~, Miss,~, ~, ~, ~] = JB_countTrialTypes(pair1);
                [~, FA, ~,CR, ~, ~, ~] = JB_countTrialTypes(pair2);
            else
                [Hit, ~, Miss,~, ~, ~, ~] = JB_countTrialTypes(pair2);
                [~, FA, ~,CR, ~, ~, ~] = JB_countTrialTypes(pair1);
            end
            basicProperties{i,1}.pairsDprime(d,1) =  JB_dPrime(Hit,Miss,CR, FA);
        else
            basicProperties{i,1}.pairsDprime(d,1) =  NaN;
        end
        pair1 = basicProperties{i,1}.trialTypesAngleSTIM{1,(1+d-1)};
        pair2 = basicProperties{i,1}.trialTypesAngleSTIM{1,end-(d-1)};
        if (~isempty(pair1) || ~isempty(pair2))
            %find which was the go and nogo
            if any(pair1==1) || any(pair1==3) %pair 1 go stim
                [Hit, ~, Miss,~, ~, ~, ~] = JB_countTrialTypes(pair1);
                [~, FA, ~,CR, ~, ~, ~] = JB_countTrialTypes(pair2);
            else
                [Hit, ~, Miss,~, ~, ~, ~] = JB_countTrialTypes(pair2);
                [~, FA, ~,CR, ~, ~, ~] = JB_countTrialTypes(pair1);
            end
            basicProperties{i,1}.pairsDprimeSTIM(d,1) =  JB_dPrime(Hit,Miss,CR, FA);
        else
            basicProperties{i,1}.pairsDprimeSTIM(d,1) =  NaN;
        end
        pair1 = basicProperties{i,1}.trialTypesAngleNoSTIM{1,(1+d-1)};
        pair2 = basicProperties{i,1}.trialTypesAngleNoSTIM{1,end-(d-1)};
        if (~isempty(pair1) || ~isempty(pair2))
            %find which was the go and nogo
            if any(pair1==1) || any(pair1==3) %pair 1 go stim
                [Hit, ~, Miss,~, ~, ~, ~] = JB_countTrialTypes(pair1);
                [~, FA, ~,CR, ~, ~, ~] = JB_countTrialTypes(pair2);
            else
                [Hit, ~, Miss,~, ~, ~, ~] = JB_countTrialTypes(pair2);
                [~, FA, ~,CR, ~, ~, ~] = JB_countTrialTypes(pair1);
            end
            basicProperties{i,1}.pairsDprimeNoSTIM(d,1) =  JB_dPrime(Hit,Miss,CR, FA);
        else
            basicProperties{i,1}.pairsDprimeSTIM(d,1) =  NaN;
        end
        basicProperties{i,1}.pairsDiff(d,1) = possibleAngles(end-(d-1))-possibleAngles(d);
    end
    
    basicProperties{i,1}.pairsDprimeNoSTIM(isnan(basicProperties{i,1}.pairsDprimeNoSTIM(:,1)),:)=[];
    basicProperties{i,1}.pairsDprimeSTIM(isnan(basicProperties{i,1}.pairsDprimeSTIM(:,1)),:)=[];
    basicProperties{i,1}.pairsDprime(isnan(basicProperties{i,1}.pairsDprime(:,1)),:)=[];
    
    %Calculate how many of each trial type there were
    [Hit, FA, Miss,CR, performance, licking, dprime] = JB_countTrialTypes(temptrialTypes); 
    basicProperties{i,1}.HIT = Hit; %HitTrialCount
    basicProperties{i,1}.FA = FA; % FATrialCount
    basicProperties{i,1}.MISS = Miss; %MissTrialCount
    basicProperties{i,1}.CR = CR; %CRTrialCount
    basicProperties{i,1}.sessionperformance = performance;
    basicProperties{i,1}.dprime=dprime;
    
    %now look at optogenetic sessions
    if (basicProperties{i,1}.optogenetics==1)
        basicProperties{i,1}.HITStim = sum((temptrialTypes==1)&(tempOptogenetics==1))+1; %HitTrialCount
        basicProperties{i,1}.FAStim = sum((temptrialTypes==2)&(tempOptogenetics==1))+1; % FATrialCount
        basicProperties{i,1}.MISSStim = sum((temptrialTypes==3)&(tempOptogenetics==1))+1; %MissTrialCount
        basicProperties{i,1}.CRStim = sum((temptrialTypes==4)&(tempOptogenetics==1))+1; %CRTrialCount
        basicProperties{i,1}.HITNoStim = sum((temptrialTypes==1)&(tempOptogenetics==0))+1; %HitTrialCount
        basicProperties{i,1}.FANoStim = sum((temptrialTypes==2)&(tempOptogenetics==0))+1; % FATrialCount
        basicProperties{i,1}.MISSNoStim = sum((temptrialTypes==3)&(tempOptogenetics==0))+1; %MissTrialCount
        basicProperties{i,1}.CRNoStim = sum((temptrialTypes==4)&(tempOptogenetics==0))+1; %CRTrialCount
        basicProperties{i,1}.dprimeSTIM=norminv(( basicProperties{i,1}.HITStim/( basicProperties{i,1}.HITStim+basicProperties{i,1}.MISSStim)),0,1)-norminv((basicProperties{i,1}.FAStim/(basicProperties{i,1}.FAStim+basicProperties{i,1}.CRStim)),0,1);
        basicProperties{i,1}.dprimeNoSTIM=norminv(( basicProperties{i,1}.HITNoStim/( basicProperties{i,1}.HITNoStim+basicProperties{i,1}.MISSNoStim)),0,1)-norminv((basicProperties{i,1}.FANoStim/(basicProperties{i,1}.FANoStim+basicProperties{i,1}.CRNoStim)),0,1);
        basicProperties{i,1}.sessionperformanceSTIM = (basicProperties{i,1}.HITStim+basicProperties{i,1}.CRStim)/(sum(tempOptogenetics==1)+4);
        basicProperties{i,1}.sessionperformanceNoSTIM = (basicProperties{i,1}.HITNoStim+basicProperties{i,1}.CRNoStim)/(sum(tempOptogenetics==0)+4);
    end
    
    %import as metaDATA cell array
    for j=1:length(metaData)
        if(strcmp('Auto reward = ', metaData{j,1}));
            autoSession = metaData{j,2};
        end
    end
    
    match = strcmp('Orientation Selected = ', metaData(:,1));
    clear ind
    
    for kk=1:length(match)
        if match(kk,1)==1;
            if isstr(metaData{kk,2})
                ind(kk,1) = str2num(metaData{kk,2});
            else
                ind(kk,1) =(metaData{kk,2});
            end
        end
    end
    tempAngles = ind(ind>0);
    
    %Define session type
    tempsessionType = strcat('S',num2str(length(tempAngles)));
    basicProperties{i,1}.stimNumber = length(tempAngles);
    if autoSession=='1'
        tempsessionType = strcat(tempsessionType,'auto');
        basicProperties{i,1}.stimNumber = 0;
    end
    basicProperties{i,1}.mouseID = DATA.mouseID;
    basicProperties{i,1}.sessionType = tempsessionType;
    
end

%delete empty cells
basicProperties(any(cellfun(@isempty,basicProperties),2),:)=[];

%Now go along and extra only a single session type from each day, e.g. if
%both S2 and S8 was run.. we look to see which file is the longest, and
%then use that
basicPropertiesToPlot = cell(100,1);

for k = 1:length(basicProperties),
    days{k,1} = basicProperties{k,1}.datenum(1:8);
end

plotNumber = 1;
[~,~,idx] = unique(days(:,1));
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
            tempRunVal(kk,1) = basicProperties{t(kk)}.fileLength;
        end
        
        [~,di] = max(tempRunVal); %find the data set where the mouse ran the most
        basicPropertiesToPlot{plotNumber,1} = basicProperties{t(di)};
        plotNumber = plotNumber+1;
    end
end

%delete empty cells
basicPropertiesToPlot(any(cellfun(@isempty,basicPropertiesToPlot),2),:)=[];
%%
JB_plotPerformance(basicPropertiesToPlot,1)
% JB_plotSessionPerformance(basicPropertiesToPlot,possibleAngles,1)
%

% [DATAavg] = JB_plotSelectionPerformance(basicPropertiesToPlot,possibleAngles,plotON,0,1,1);
% %JB_plotSessionPerformance(basicPropertiesToPlot,possibleAngles,plotON)
% JB_plotOptogenetics(basicPropertiesToPlot,possibleAngles,plotON)
% JB_plotTrimming(basicPropertiesToPlot, plotON)
% JB_plotSegments(basicPropertiesToPlot,plotON)
%

end
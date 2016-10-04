function [data,latencytempCatch] = JB_trialPerformance(DATA,sessionNo,stim,plotAngles,includeVelocityMap,cleanUp, noTrialsToPlot)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% analysis stop early criteria

% %%example
% sessionNo =46;
% includeVelocityMap = 0;
% cleanUp = 0;
% stim =1;
% plotAngles=0;
% noTrialsToPlot=[]; %plot only a subset of trials - these will be the middle trials unless otherwise specified.
% load('DATA.mat');

%%
minNoTrialsTotal = 200; % no of trials required to have been completed for this criteria to be used
stopExThreshold = 32; % no of consequtive misses and correct rejections to say the mouse wasn't engaged and stop analysis (4xS8)

% load DATA.mat
positionGraph1 = [1501 71 384 902];
positionGraph2 = [680 288 314 690];
positionGraph3 = [2 42 958 954];
positionGraph4 = [1400 36 485 954];

plotColor = {'Red','Crimson','Salmon','Pink','LightSteelBlue','LightSkyBlue','DodgerBlue','MediumBlue'};
PreStim=2000; %time stamps are in 10ms bins, so 100 = 1 second 200 = 2seconds
PostStim=2000;
figure(99);clf;

tempDATA = DATA.allFiles{1,sessionNo}.rawData;
metaData = DATA.allFiles{1,sessionNo}.metaData;

%find out if this session has positive or neg GO stim
orgGOIdx = find(strcmp('orgGOangleTab = ', metaData));

if isempty(orgGOIdx)
    orgGOstim = 1;
else
    orgGOstim = metaData{orgGOIdx,2};
    if ischar(orgGOstim)==1
        orgGOstim = str2num(orgGOstim);
    else
        orgGOstim = orgGOstim;
    end
end

% specify what type of session we are going to be analysing
if (exist('stim','var') && stim==0)
    varToPlot = {'fullSession'};
elseif (exist('stim','var') && stim==1)
    varToPlot = {'stimTrials' 'cntrTrials'};
elseif (exist('stim','var') && stim==2)
    varToPlot = {'GoTrials' 'NoGoTrials'};
else
    varToPlot = {'GoTrials' 'NoGoTrials'};
end

if (exist('plotAngles','var'))
else
    plotAngles = 0;
end

if (exist('includeVelocityMap','var'))
else
    includeVelocityMap = 0;
end
if (exist('cleanUp','var'))
else
    cleanUp = 0;
end

if (exist('noTrialsToPlot','var'))
    noTrialsToPlot;
else
    noTrialsToPlot = [];
end

% find out if this was a 2AFC
AFCIdx = find(strcmp('2ACF = ', metaData));
if ~isempty(AFCIdx)
    AFCsession = 1;
else
    AFCsession = 0;
end

% find out if this was an autoReward Session
autoIdx = find(strcmp('Auto reward = ', metaData));
if ~isempty(autoIdx)
    data.autoSession = str2num(metaData{autoIdx,2});
else
    data.autoSession = 0;
end

%define which txt file column the different digital inputs are read in from
%processing
encoder0Pos = 1;
targetArea = 2;
reward = 3;
licks = 4;
LEDstate = 5;
count = 6;
millis = 7;
maskingLight = 9;
angle = 11;
catchTrial = 12;
trialType = 13;
licksR = 16;
rewardR = 17;

midPoint=270;
threshold = 50000;
aboveThreshold = find(tempDATA(:,1)>threshold);

% if mouse ran backwards before starting the trial, the step number will be
% 65540 - delete these points

for i = 1:length(aboveThreshold)
    tempDATA(aboveThreshold(i),:)=nan;
end
tempDATA = tempDATA(~any(isnan(tempDATA),2),:);

%Delete any blips in the recordings - missaligned txt input - find these by
%searching for non binary values in the columns where these are expected
[ind, ~] = find(tempDATA(:,[2,3,4,5,8])>1); %all these inputs should be binary
if ~isempty(ind) %if there was a blip, delete this row
    UniqInd = unique(ind);
    for n=1:length(UniqInd);
        tempDATA(UniqInd(n),:)=nan;
    end
    tempDATA = tempDATA(~any(isnan(tempDATA),2),:);
end

% Delete any blips where the clock went backwards
negTime = find(diff(tempDATA(:,7))<=0) ;
for i = 1:length(negTime)
    tempDATA(negTime(i),:)=nan;
end
tempDATA = tempDATA(~any(isnan(tempDATA),2),:);

% Delete any blips where the clock didn't increase
negTime = find(diff(tempDATA(:,7))==0);
for i = 1:length(negTime)
    tempDATA(negTime(i),:)=nan;
end
tempDATA = tempDATA(~any(isnan(tempDATA),2),:);

%extract all the angles used in this session
C = unique(tempDATA(:,angle));
data.angleOrientations = C(C>100);

%% Divide up data into single revolutions/trials

revolutionsDATA = find((diff(tempDATA(:,1)))<-100); % find where encoder is reset back to 0 to mark the end of a trial/start of enw trial
for j=1:length(revolutionsDATA)-1; %Dont include first revolution as wont be properly callibrated.
    data.trial{j}.raw(:,:) = tempDATA(revolutionsDATA(j)+1:revolutionsDATA(j+1),:);
end

% Check Thats Detection of stimulus worked
toDelete = [];
p=1;
for j = 1:length(data.trial);
    [loc2,~] = find((data.trial{j}.raw(:,targetArea)>0)); %did the target enter the target region - should do on all trials, if not - delete
    if isempty(loc2) || loc2(1)==1;
        toDelete(p) = j;
        p=p+1;
    end
end

if ~isempty(toDelete)
    data.trial(toDelete)=[];
end

% Check Thats target entered and exsited stim area
toDelete = [];
p=1;
for j = 1:length(data.trial);
    
    [startTrial, ~] = find(diff(data.trial{j}.raw(:,targetArea))>0);
    [endTrial, ~] = find(diff(data.trial{j}.raw(:,targetArea))<0);
    
    if isempty(startTrial) || isempty(endTrial);
        toDelete(p) = j;
        p=p+1;
    end
end

if ~isempty(toDelete)
    data.trial(toDelete)=[];
end

%% Define Target Time, Position, Correct (Hit, correct Rejection) or Incorrect (Miss, or False Alarm) data.trial

sampleRateSpacing = 1; %1ms sampling
hitCount = 0;
crCount = 0;
hitCountCatch = 0;
crCountCatch = 0;
stopIdx=0;

for j = 1: length(data.trial);
    % make sample times uniform and zeroed
    SampleTimes = data.trial{j}.raw(:,millis)-data.trial{j}.raw(1,millis);
    clear tempFiltered;
    
    %make each data point uniformly sampled
    for i = 1:size((data.trial{j}.raw),2);
        tempTrial = data.trial{j}.raw(:,i);
        filteredSampleTimes = 0:sampleRateSpacing:max(SampleTimes); %1msec spacing
        tempFiltered(:,i)=interp1(SampleTimes,tempTrial,filteredSampleTimes);
    end
    
    % make sure all values are integers
    data.trial{j}.filteredDATA = [tempFiltered(:,1) (round(tempFiltered(:,2:end)))];
    
    %LEFT LICK
    % make sure licks are represented as just one tick mark
    B = diff(data.trial{j}.filteredDATA(:,licks))==1;
    B(numel(data.trial{j}.filteredDATA(:,licks))) = 0;
    data.trial{j}.filteredDATA(:,licks)=B;
    
    %RIGHT LICK
    % make sure right licks are represented as just one tick mark
    B = diff(data.trial{j}.filteredDATA(:,licksR))==1;
    B(numel(data.trial{j}.filteredDATA(:,licksR))) = 0;
    data.trial{j}.filteredDATA(:,licksR)=B;
    
    data.trial{j}.timeStamps(:,1) = [0;(cumsum(diff(data.trial{j}.filteredDATA(:,millis)))/1000)];
    
    %%NEW
    [startTrial(j,1), ~] = find(diff(data.trial{j}.filteredDATA(:,targetArea))>0);
    [endTrial(j,1), ~] = find(diff(data.trial{j}.filteredDATA(:,targetArea))<0);
    tempAngleTypes = data.trial{j}.filteredDATA(1,angle);
    data.trial{j}.Performance.StimType =tempAngleTypes;
    
    data.trial{j}.Stimulus.AnswerPeriodStartIdx =startTrial(j,1);
    data.trial{j}.Stimulus.AnswerPeriodEndIdx = endTrial(j,1);
    
    clear temptrialTypes;
    clear temptrialTypesCatch;
    
    
    if ~isempty(endTrial(j,1))
        lick = find((data.trial{j}.filteredDATA(startTrial(j,1):endTrial(j,1),licks))==1);
        if isempty(lick);
            licksTrial=0;
        else
            licksTrial=1;
        end
        if ((tempAngleTypes<270) &&  (orgGOstim==1)) %origional Go trial
            if (licksTrial==1) %HIT
                temptrialTypes = 1;
            elseif (licksTrial==0) %MISS
                temptrialTypes = 3;
            end
        elseif ((tempAngleTypes>270) &&  (orgGOstim==0)) %switched NoGo trial
            if (licksTrial==1) %False Alarm
                temptrialTypes = 2;
            elseif (licksTrial==0) %Correct Rejection
                temptrialTypes = 4;
            end
        elseif ((tempAngleTypes>270) &&  (orgGOstim==1))%origional NoGo trial
            if (licksTrial==1) %False Alarm
                temptrialTypes = 2;
            elseif (licksTrial==0) %Correct Rejection
                temptrialTypes = 4;
            end
        elseif ((tempAngleTypes>270) &&  (orgGOstim==0))%switched Go trial
            if (licksTrial==1) %HIT
                temptrialTypes = 1;
            elseif (licksTrial==0) %MISS
                temptrialTypes= 3;
            end
        end
    end
    
    data.trial{j}.Performance.StimType = data.trial{j}.filteredDATA(1,angle);
    data.trial{j}.Performance.TrialType = temptrialTypes;
    data.trial{j}.Performance.TrialCatch = data.trial{j}.filteredDATA(startTrial(j,1),catchTrial);
    %  data.trial{j}.Performance.TrialType = data.trial{j}.raw(find(data.trial{j}.raw(:,trialType)>0),trialType);
    
    %  [loc2,~] = find((data.trial{j}.filteredDATA(:,targetArea)>0)); % Identify when stim entered TargetArea
    rewardIdx = find((data.trial{j}.filteredDATA(:,reward)==1)); % Identify RewardTime if Hit data.trial
    rewardRIdx = find((data.trial{j}.filteredDATA(:,rewardR)==1)); % Identify RewardTime if Hit data.trial
    
    % If Reward was given, update information on record properties
    if ~isempty(rewardIdx); %
        rewardIdx = rewardIdx(1);
    elseif ~isempty(rewardRIdx);
        rewardIdx = rewardRIdx(1);
    else
        rewardIdx = nan;
    end
    
    %Gather information about Answer Period
    answerPeriodDuration = (endTrial(j,1)-startTrial(j,1))*sampleRateSpacing;
    
    % If Optogenetic Stim turned on, identify ON and OFF time
    tempStimIdx = find((data.trial{j}.filteredDATA(:,LEDstate)==1));
    tempStimTime(j,1) = NaN;
    if ~isempty(tempStimIdx); %
        data.trial{j}.Performance.StimTrial = 1;
        tempStimTime(j,1) = data.trial{j}.filteredDATA(tempStimIdx(1),1)-data.trial{j}.filteredDATA(startTrial(j,1),1);
        stimIdx(j,1) = tempStimIdx(1);
    else
        data.trial{j}.Performance.StimTrial=0;
        stimIdx(j,1) = nan;
    end
    
    % If masking light was used identify ON and OFF time
    tempMaskingIdx = find((data.trial{j}.filteredDATA(:,maskingLight)==1));
    tempMaskingTime(j,1) = NaN;
    if ~isempty(tempMaskingIdx); %
        data.trial{j}.Performance.MaskTrial = 1;
        tempMaskingTime(j,1) = data.trial{j}.filteredDATA(tempMaskingIdx(1),1)-data.trial{j}.filteredDATA(startTrial(j,1),1);
        maskIdx(j,1) = tempMaskingIdx(1);
    else
        data.trial{j}.Performance.MaskTrial=0;
        maskIdx(j,1) = nan;
    end
    
    % Determine if animal licked during Answer Period
    if ~isnan(rewardIdx) %If there was a reward given
        data.trial{j}.Performance.LickFrequency = sum(data.trial{j}.filteredDATA(startTrial(j,1):rewardIdx,licks)==1)/(answerPeriodDuration/1000);
        tempFirstLickIdx = (startTrial(j,1)+(find((data.trial{j}.filteredDATA(startTrial(j,1):rewardIdx,licks)==1),1)));
        tempLickLatency = data.trial{j}.timeStamps(tempFirstLickIdx)-data.trial{j}.timeStamps(startTrial(j,1));
        if ~isempty(tempLickLatency)
            data.trial{j}.Performance.FirstLickLatency = tempLickLatency;
        else
            data.trial{j}.Performance.FirstLickLatency = nan;
        end
    else
        if length(data.trial{j}.filteredDATA)>endTrial(j,1);%if reward port didn't open find out if any licks occured during whloe answer period
            maxTimeIdx = endTrial(j,1);
        else
            maxTimeIdx = length(data.trial{j}.filteredDATA);
        end
        data.trial{j}.Performance.LickFrequency = sum(data.trial{j}.filteredDATA(startTrial(j,1):maxTimeIdx,licks)==1)/(answerPeriodDuration/1000);
        tempFirstLickIdx = (startTrial(j,1)+(find((data.trial{j}.filteredDATA(startTrial(j,1):maxTimeIdx,licks)==1),1)));
        tempLickLatency = data.trial{j}.timeStamps(tempFirstLickIdx)-data.trial{j}.timeStamps(startTrial(j,1));
        if ~isempty(tempLickLatency)
            data.trial{j}.Performance.FirstLickLatency = tempLickLatency;
        else
            data.trial{j}.Performance.FirstLickLatency = nan;
        end
    end
    
    % Update Trial Type Results
    if (mode(data.trial{j}.Performance.TrialType) == 1); % If this was the GoTrials stimulus - Hit trial
        if (data.trial{j}.Performance.TrialCatch==0);
            hitCount = hitCount+1;
        end
        data.trial{j}.Performance.OutcomeLick = 1; % Anticipatory lick - Hit trial
        data.trial{j}.Performance.OutcomeCode = 1;
    elseif (mode(data.trial{j}.Performance.TrialType) == 2)% If this was a distractor stimulus - false alarm
        data.trial{j}.Performance.OutcomeLick = 1;% False Alarm - incorrect lick
        data.trial{j}.Performance.OutcomeCode = 2;
    elseif (data.trial{j}.Performance.StimType == midPoint) % If this was exactly in the middle
        data.trial{j}.Performance.OutcomeLick = 1;% licked for 50% verticle
        data.trial{j}.Performance.OutcomeCode = 5; % 5 if licked, 6 if didn't
    elseif (mode(data.trial{j}.Performance.TrialType) == 3); % If this was the no GoTrials stimulus
        data.trial{j}.Performance.OutcomeLick = 0;  % No anticipatory lick - Miss trial
        data.trial{j}.Performance.OutcomeCode = 3;
    elseif (mode(data.trial{j}.Performance.TrialType) == 4) % If this was a distractor stimulus and no lick
        if (data.trial{j}.Performance.TrialCatch==0);
            crCount = crCount+1;
        end
        data.trial{j}.Performance.OutcomeLick = 0; % Correctly withheld lick
        data.trial{j}.Performance.OutcomeCode = 4;
    elseif (data.trial{j}.Performance.StimType == midPoint) % If this was exactly in the middle
        data.trial{j}.Performance.OutcomeLick = 0;% licked for 50% verticle
        data.trial{j}.Performance.OutcomeCode = 6; % 5 if licked, 6 if didn't
    elseif ((data.autoSession == 1) && (data.trial{j}.Performance.StimType==225)) %if this was a autorewards session
        data.trial{j}.Performance.OutcomeLick = 1;% licked for 50% verticle
        data.trial{j}.Performance.OutcomeCode = 1; % 5 if licked, 6 if didn't
    elseif ((data.autoSession == 1) && (data.trial{j}.Performance.StimType==315)) %if this was a autorewards session
        data.trial{j}.Performance.OutcomeLick = 1;% licked for 50% verticle
        data.trial{j}.Performance.OutcomeCode = 4; % 5 if licked, 6 if didn't
    end
    cumulativePerformance(j,1) = ((hitCount+crCount)/j);
end

tempStimTime(tempStimTime==0)=NaN;
meanStimONSET = nanmean(tempStimTime);

tempMaskingTime(tempMaskingTime==0)=NaN;
meanMaskONSET = nanmean(tempMaskingTime);

%get the average Stim time onset
for j = 1: length(data.trial);
    %         [startTrial, ~] = find(diff(data.trial{j}.filteredDATA(:,targetArea))>0);
    %         [endTrial, ~] = find(diff(data.trial{j}.filteredDATA(:,targetArea))<0);
    %          stimIdx = find((data.trial{j}.filteredDATA(:,LEDstate)==1));
    tempLickLatencySTIM=[];
    if length(data.trial{j}.filteredDATA)>endTrial(j,1);%if reward port didn't open find out if any licks occured during whloe answer period
        maxTimeIdx = endTrial(j,1);
    else
        maxTimeIdx = length(data.trial{j}.filteredDATA);
    end
    
    if any(data.trial{j}.raw(:,LEDstate),1); %USE REAL STIM ON TIME
        tempFirstLickIdxSTIM = (stimIdx(j,1)+(find((data.trial{j}.filteredDATA(stimIdx(j,1):maxTimeIdx,licks)==1),1)));
        if ~isempty(tempFirstLickIdxSTIM)
            tempLickLatencySTIM = data.trial{j}.timeStamps(tempFirstLickIdxSTIM(1))-data.trial{j}.timeStamps(stimIdx(j,1));
            if ~isempty(tempLickLatencySTIM)
                data.trial{j}.Performance.FirstLickLatencySTIM = tempLickLatencySTIM;
            end
            
        else
            data.trial{j}.Performance.FirstLickLatencySTIM = nan;
        end
    else
        answerTime = data.trial{j}.timeStamps(startTrial(j,1));
        estStimTime = answerTime+(meanStimONSET/1000);
        estStimTimeIdx = find(data.trial{j}.timeStamps>=estStimTime,1);
        tempFirstLickIdxSTIM = (estStimTimeIdx+(find((data.trial{j}.filteredDATA(estStimTimeIdx:maxTimeIdx,licks)==1),1)));
        %         tempLickLatencySTIM = data.trial{j}.timeStamps(tempFirstLickIdxSTIM)-data.trial{j}.timeStamps(estStimTimeIdx);
        if ~isempty(tempLickLatencySTIM)
            data.trial{j}.Performance.FirstLickLatencySTIM = data.trial{j}.timeStamps(tempFirstLickIdxSTIM)-data.trial{j}.timeStamps(estStimTimeIdx);
        else
            data.trial{j}.Performance.FirstLickLatencySTIM = nan;
        end
        
    end
end

%%Get latency to first lick from masking light on
for j = 1: length(data.trial);
    if any(data.trial{j}.raw(:,maskingLight),1); %USE REAL masking light ON TIME
        maskingLightIdx = find(diff((data.trial{j}.filteredDATA(:,maskingLight))>0),1);
        targetOFFIdx = find(diff(data.trial{j}.filteredDATA(:,targetArea))<0);
        tempFirstLickIdxMask = maskingLightIdx+(find((data.trial{j}.filteredDATA(maskingLightIdx:targetOFFIdx,licks)==1),1));
        if ~isempty(tempFirstLickIdxMask)
            tempLickLatencyMask = data.trial{j}.timeStamps(tempFirstLickIdxMask)-data.trial{j}.timeStamps(maskingLightIdx);
            data.trial{j}.Performance.FirstLickLatencyMask = tempLickLatencyMask;
        else
            data.trial{j}.Performance.FirstLickLatencyMask = nan;
        end
        
    else
        data.trial{j}.Performance.FirstLickLatencyMask = nan;
    end
end


for j=1:length(data.trial);
    data.CumulativePerformance(j,1) = data.trial{j}.Performance.StimType;
    data.CumulativePerformance(j,2) = data.trial{j}.Performance.OutcomeLick;
    data.CumulativePerformance(j,3) = data.trial{j}.Performance.LickFrequency;
    data.CumulativePerformance(j,4) = data.trial{j}.Performance.FirstLickLatency;
    data.CumulativePerformance(j,5) = data.trial{j}.Performance.OutcomeCode;
    data.CumulativePerformance(j,6) = data.trial{j}.Performance.StimTrial;
    data.CumulativePerformance(j,7) = data.trial{j}.Performance.FirstLickLatencySTIM;
    data.CumulativePerformance(j,8) = data.trial{j}.Performance.TrialCatch;
    data.CumulativePerformance(j,9) = data.trial{j}.Performance.FirstLickLatencyMask;
end


%find how many consequtive trials with no lick there were
ff = (data.CumulativePerformance(minNoTrialsTotal:end,5)==3) | (data.CumulativePerformance(minNoTrialsTotal:end,5)==4); %miss and correct rejection
tempCum = 1;
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
    idxStop = (find(cumCount>stopExThreshold))+(minNoTrialsTotal);
    if ~isempty(idxStop) && (cleanUp == 1)
        stopIdx = idxStop(1)-stopExThreshold;
        if stopIdx>0
            X = [' Stopping Session Analysis Early By', num2str(length(data.CumulativePerformance)- stopIdx), 'Trials', 'On trial No.', num2str(stopIdx),'Out Of',num2str(length(data.CumulativePerformance))] ;
            disp(X)
        end
    else
        stopIdx=0;
    end
end


%% Test Pad data.trial and align to reward onset TEST

fullSampleTime = PreStim+PostStim;
TimePointsMin = -PreStim:1:PostStim-1;
% varToAlign = {'AnswerPeriodEndIdx'};
varToAlign = {'AnswerPeriodStartIdx'};
%varToAlign = {'AnswerPeriodStartIdx' 'AnswerPeriodEndIdx'};

data.AlignedDATA.TrialType.HIT = [];
data.AlignedDATA.TrialType.MISS = [];
data.AlignedDATA.TrialType.CorrectRejection = [];
data.AlignedDATA.TrialType.FalseAlarm = [];
data.AlignedDATA.TrialType.HITStim = [];
data.AlignedDATA.TrialType.MISSStim = [];
data.AlignedDATA.TrialType.CorrectRejectionStim = [];
data.AlignedDATA.TrialType.FalseAlarmStim = [];
data.AlignedDATA.TrialType.HITNoStim = [];
data.AlignedDATA.TrialType.MISSNoStim = [];
data.AlignedDATA.TrialType.CorrectRejectionNoStim = [];
data.AlignedDATA.TrialType.FalseAlarmNoStim = [];
data.stimTrials = [];

for f = 1:length(varToAlign)
    for k = 1:length(data.angleOrientations) % create empty field for each stimulus orientation
        orientationfieldname{k,1} = ['stimAngle' int2str(data.angleOrientations(k))];
        data.(orientationfieldname{k,1}).(varToAlign{f}) = [];
        HCRorientationfieldname{k,1} = ['HCRstimAngle' int2str(data.angleOrientations(k))];
        data.(HCRorientationfieldname{k,1}).(varToAlign{f}) = [];
        MFAorientationfieldname{k,1} = ['MFAstimAngle' int2str(data.angleOrientations(k))];
        data.(MFAorientationfieldname{k,1}).(varToAlign{f}) = [];
    end
    
    clear data.fullSession.(varToAlign{f});
    %Counters for Performance outcome
    h1=1;
    m1=1;
    cr1=1;
    fa1=1;
    %Counters for Glabal Stimulus type
    d1=1;
    s1=1;
    mp1=1;
    mnl1=1;
    ml1=1;
    auto1=1;
    auto11=1;
    auto111=1;
    
    %Counters for Performance outcome with and without Stimulation
    hStim1 = 1;
    hNoStim1 = 1;
    mStim1 = 1;
    mNoStim1 = 1;
    crStim1 = 1;
    crNoStim1 = 1;
    faStim1 = 1;
    faNoStim1 = 1;
    
    %Counters for Trials with and without Stimulation
    stimON = 1;
    stimOFF = 1;
    stimONCatch = 1;
    stimOFFCatch = 1;
    gg = 1;
    jj=1;
    block=1;
    
    for j=1:length(data.trial)-stopIdx;
        raster_matrix = [];
        velocityTemp = diff(data.trial{j}.filteredDATA(:,1));
        
        addPreStim = [];
        addPostStim = [];
        
        if (data.trial{j}.Stimulus.(varToAlign{f})-PreStim<1)
            addPreStim=NaN*ones(1,PreStim - (data.trial{j}.Stimulus.(varToAlign{f})-1)); %%add nans to start of vector
            currPreStim =  PreStim-(PreStim - data.trial{j}.Stimulus.(varToAlign{f})+1);
        else
            currPreStim = PreStim;
        end
        
        if length(velocityTemp)<data.trial{j}.Stimulus.(varToAlign{f})+PostStim-1
            addPostStim=NaN*ones(1,data.trial{j}.Stimulus.(varToAlign{f})+PostStim-1-length(velocityTemp)); %%add nans to start of vector
            currPostStim =  PostStim-((data.trial{j}.Stimulus.(varToAlign{f})+PostStim-1)-length(velocityTemp));
        else
            currPostStim = PostStim;
        end
        
        %Get timestamps from End of sample time (target time)
        %ALL data
        tempVelocityStim = velocityTemp(data.trial{j}.Stimulus.(varToAlign{f})-currPreStim:data.trial{j}.Stimulus.(varToAlign{f})+currPostStim-1,encoder0Pos)';
        tempLickStim = (data.trial{j}.filteredDATA(data.trial{j}.Stimulus.(varToAlign{f})-currPreStim:data.trial{j}.Stimulus.(varToAlign{f})+currPostStim-1,licks))';
        IdxTarget = (data.trial{j}.filteredDATA(data.trial{j}.Stimulus.(varToAlign{f})-currPreStim:data.trial{j}.Stimulus.(varToAlign{f})+currPostStim-1,targetArea))';
        IdxReward = (data.trial{j}.filteredDATA(data.trial{j}.Stimulus.(varToAlign{f})-currPreStim:data.trial{j}.Stimulus.(varToAlign{f})+currPostStim-1,reward))';
        IdxStim = (data.trial{j}.filteredDATA(data.trial{j}.Stimulus.(varToAlign{f})-currPreStim:data.trial{j}.Stimulus.(varToAlign{f})+currPostStim-1,LEDstate))';
        IdxRewardR = (data.trial{j}.filteredDATA(data.trial{j}.Stimulus.(varToAlign{f})-currPreStim:data.trial{j}.Stimulus.(varToAlign{f})+currPostStim-1,rewardR))';
        tempLickRStim = (data.trial{j}.filteredDATA(data.trial{j}.Stimulus.(varToAlign{f})-currPreStim:data.trial{j}.Stimulus.(varToAlign{f})+currPostStim-1,licksR))';
        IdxMask = (data.trial{j}.filteredDATA(data.trial{j}.Stimulus.(varToAlign{f})-currPreStim:data.trial{j}.Stimulus.(varToAlign{f})+currPostStim-1,maskingLight))';
        
        %  if (data.trial{j}.Performance.TrialCatch==0)
        data.fullSession.(varToAlign{f})(1:fullSampleTime,gg,1) = [addPreStim tempVelocityStim addPostStim];
        data.fullSession.(varToAlign{f})(1:fullSampleTime,gg,2) = [addPreStim tempLickStim addPostStim];
        data.fullSession.(varToAlign{f})(1,gg,3) = data.trial{j}.Performance.OutcomeLick;
        data.fullSession.(varToAlign{f})(1,gg,4) = data.trial{j}.Performance.StimType;
        data.fullSession.(varToAlign{f})(1:fullSampleTime,gg,5) = [addPreStim IdxTarget addPostStim];
        data.fullSession.(varToAlign{f})(1:fullSampleTime,gg,6) = [addPreStim IdxReward addPostStim];
        data.fullSession.(varToAlign{f})(1,gg,7) = data.trial{j}.Performance.OutcomeCode;
        data.fullSession.(varToAlign{f})(1:fullSampleTime,gg,8) = [addPreStim IdxStim addPostStim];
        data.fullSession.(varToAlign{f})(1:fullSampleTime,gg,9) = [addPreStim IdxRewardR addPostStim];
        data.fullSession.(varToAlign{f})(1:fullSampleTime,gg,10) = [addPreStim tempLickRStim addPostStim];
        data.fullSession.(varToAlign{f})(1,gg,11) = data.trial{j}.Performance.TrialCatch;
        data.fullSession.(varToAlign{f})(1:fullSampleTime,gg,12) = [addPreStim IdxMask addPostStim];
        gg=gg+1;
        
        %  for j=1:length(data.trial)-stopIdx; %only non catch trials
        %  if (data.trial{j}.Performance.TrialCatch==0)
        
        
        blockPerformance(1,block) = data.fullSession.(varToAlign{f})(1,j,7);
        block=block+1;
        
        if (data.trial{j}.Performance.StimTrial ==1); %Stimulated Session
            data.stimTrials.(varToAlign{f})(:,stimON,:) = data.fullSession.(varToAlign{f})(:,j,:);
            stimON = stimON+1;
        else
            data.cntrTrials.(varToAlign{f})(:,stimOFF,:) = data.fullSession.(varToAlign{f})(:,j,:);
            stimOFF = stimOFF+1;
        end
        
        if (data.autoSession == 1)
            if (data.trial{j}.Performance.StimType==225) %Left Lick
                data.GoTrials.(varToAlign{f})(:,auto1,:) = data.fullSession.(varToAlign{f})(:,j,:);
                auto1 = auto1+1;
            elseif (data.trial{j}.Performance.StimType==315) %Right Lick
                data.NoGoTrials.(varToAlign{f})(:,auto11,:) = data.fullSession.(varToAlign{f})(:,j,:);
                auto11 = auto11+1;
            end
            
        elseif (data.trial{j}.Performance.TrialType(1) == 1 || (data.trial{j}.Performance.TrialType(1) == 3)); %GoTrials Stim
            data.GoTrials.(varToAlign{f})(:,s1,:) = data.fullSession.(varToAlign{f})(:,j,:);
            s1=s1+1;
            if data.trial{j}.Performance.OutcomeLick==1; %HIT trial
                data.AlignedDATA.TrialType.HIT.(varToAlign{f})(:,h1,:) = data.fullSession.(varToAlign{f})(:,j,:);
                h1 = h1+1;
                if data.trial{j}.Performance.StimTrial==1; %Stimulation
                    data.AlignedDATA.TrialType.HITStim.(varToAlign{f})(:,hStim1,:) = data.fullSession.(varToAlign{f})(:,j,:);
                    hStim1 = hStim1+1;
                else
                    data.AlignedDATA.TrialType.HITNoStim.(varToAlign{f})(:,hNoStim1,:) = data.fullSession.(varToAlign{f})(:,j,:);
                    hNoStim1 = hNoStim1+1;
                end
            elseif data.trial{j}.Performance.OutcomeLick==0 %Miss trial
                data.AlignedDATA.TrialType.MISS.(varToAlign{f})(:,m1,:) = data.fullSession.(varToAlign{f})(:,j,:);
                m1 = m1+1;
                if data.trial{j}.Performance.StimTrial==1; %Stimulation
                    data.AlignedDATA.TrialType.MISSStim.(varToAlign{f})(:,mStim1,:) = data.fullSession.(varToAlign{f})(:,j,:);
                    mStim1 = mStim1+1;
                else
                    data.AlignedDATA.TrialType.MISSNoStim.(varToAlign{f})(:,mNoStim1,:) = data.fullSession.(varToAlign{f})(:,j,:);
                    mNoStim1 = mNoStim1+1;
                end
            end
        elseif (data.trial{j}.Performance.TrialType(1) ==2 || (data.trial{j}.Performance.TrialType(1) == 4)); %NoGoTrials Stim
            data.NoGoTrials.(varToAlign{f})(:,d1,:) = data.fullSession.(varToAlign{f})(:,j,:);
            d1=d1+1;
            if data.trial{j}.Performance.OutcomeLick==0; %Correct Rejection Trial
                data.AlignedDATA.TrialType.CorrectRejection.(varToAlign{f})(:,cr1,:) = data.fullSession.(varToAlign{f})(:,j,:);
                cr1 = cr1+1;
                if data.trial{j}.Performance.StimTrial==1; %Stimulation
                    data.AlignedDATA.TrialType.CorrectRejectionStim.(varToAlign{f})(:,crStim1,:) = data.fullSession.(varToAlign{f})(:,j,:);
                    crStim1 = crStim1+1;
                else
                    data.AlignedDATA.TrialType.CorrectRejectionNoStim.(varToAlign{f})(:,crNoStim1,:) = data.fullSession.(varToAlign{f})(:,j,:);
                    crNoStim1 = crNoStim1+1;
                end
                
            elseif data.trial{j}.Performance.OutcomeLick==1; %false alarm trial
                data.AlignedDATA.TrialType.FalseAlarm.(varToAlign{f})(:,fa1,:) = data.fullSession.(varToAlign{f})(:,j,:);
                fa1 = fa1+1;
                if data.trial{j}.Performance.StimTrial==1; %Stimulation
                    data.AlignedDATA.TrialType.FalseAlarmStim.(varToAlign{f})(:,faStim1,:) = data.fullSession.(varToAlign{f})(:,j,:);
                    faStim1 = faStim1+1;
                else
                    data.AlignedDATA.TrialType.FalseAlarmNoStim.(varToAlign{f})(:,faNoStim1,:) = data.fullSession.(varToAlign{f})(:,j,:);
                    faNoStim1 = faNoStim1+1;
                end
            end
            
        elseif (data.trial{j}.Performance.StimType == midPoint); %midpoint
            data.AlignedDATAallDATAmidPoint.(varToAlign{f})(:,mp1,:) = data.fullSession.(varToAlign{f})(:,j,:);
            mp1=mp1+1;
            if data.trial{j}.Performance.OutcomeLick==0; %didn't lick
                data.AlignedDATA.TrialType.midPointNoLickCount.(varToAlign{f})(:,mnl1,:) = data.fullSession.(varToAlign{f})(:,j,:);
                mnl1 = mnl1+1;
            elseif data.trial{j}.Performance.OutcomeLick==1; %false alarm trial
                data.AlignedDATA.TrialType.midPointLickCount.(varToAlign{f})(:,ml1,:) = data.fullSession.(varToAlign{f})(:,j,:);
                ml1 = ml1+1;
            end
            
            %                     elseif (data.autoSession == 1); %autoReward Session
            %             data.AlignedDATAallDATAmidPoint.(varToAlign{f})(:,auto1,:) = data.fullSession.(varToAlign{f})(:,j,:);
            %             auto1=auto1+1;
            %             if data.trial{j}.Performance.OutcomeLick==0; %didn't lick
            %                 data.AlignedDATA.TrialType.midPointNoLickCount.(varToAlign{f})(:,auto11,:) = data.fullSession.(varToAlign{f})(:,j,:);
            %                 auto11 = auto11+1;
            %             elseif data.trial{j}.Performance.OutcomeLick==1; %false alarm trial
            %                 data.AlignedDATA.TrialType.midPointLickCount.(varToAlign{f})(:,auto111,:) = data.fullSession.(varToAlign{f})(:,j,:);
            %                 auto111 = auto111+1;
            %             end
        end
        
        for k = 1:length(data.angleOrientations)-stopIdx
            if data.trial{j}.Performance.StimType==data.angleOrientations(k); %
                data.(orientationfieldname{k}).(varToAlign{f})(:,size((data.(orientationfieldname{k}).(varToAlign{f})),2)+1,:) = data.fullSession.(varToAlign{f})(:,j,:);
            end
            
            if ((data.trial{j}.Performance.StimType==data.angleOrientations(k)) && ((data.trial{j}.Performance.OutcomeCode==1) || data.trial{j}.Performance.OutcomeCode==4)); %
                data.(HCRorientationfieldname{k}).(varToAlign{f})(:,size((data.(HCRorientationfieldname{k}).(varToAlign{f})),2)+1,:) = data.fullSession.(varToAlign{f})(:,j,:);
            end
            
            if ((data.trial{j}.Performance.StimType==data.angleOrientations(k)) && ((data.trial{j}.Performance.OutcomeCode==2) || data.trial{j}.Performance.OutcomeCode==3)); %
                data.(MFAorientationfieldname{k}).(varToAlign{f})(:,size((data.(MFAorientationfieldname{k}).(varToAlign{f})),2)+1,:) = data.fullSession.(varToAlign{f})(:,j,:);
            end
        end
        
        
        %         elseif (data.trial{j}.Performance.TrialCatch==1) %catch trials
        %             if (data.trial{j}.Performance.StimTrial ==1); %Stimulated Session
        %                 data.stimTrialsCatch.(varToAlign{f})(:,stimONCatch,:) = data.fullSession.(varToAlign{f})(:,j,:);
        %                 stimONCatch = stimONCatch+1;
        %             else
        %                 data.cntrTrialsCatch.(varToAlign{f})(:,stimOFFCatch,:) = data.fullSession.(varToAlign{f})(:,j,:);
        %                 stimOFFCatch = stimOFFCatch+1;
        %             end
        %         end
    end
    
    data = orderfields(data);
    
    %Calculate Block Performance
    blockBinSize = 12;
    blockPerformance(blockPerformance==4)=1;
    blockPerformance(blockPerformance==1)=1;
    blockPerformance(blockPerformance==2)=2;
    blockPerformance(blockPerformance==3)=2;
    binspacing = 1:blockBinSize:size((blockPerformance),2);
    
    for s = 1:length(binspacing)-1;
        blockPerformancePlot(s) = size((find(blockPerformance(1,binspacing(s):binspacing(s+1)-1)==1)),2)/blockBinSize;
    end
    
    if (plotAngles==1)
        for tt = 1:length(orientationfieldname)
            varToPlot{1,length(varToPlot)+1} = orientationfieldname{tt};
        end
        
        for tt = 1:length(HCRorientationfieldname)
            varToPlot{1,length(varToPlot)+1} = HCRorientationfieldname{tt};
        end
        
        for tt = 1:length(MFAorientationfieldname)
            varToPlot{1,length(varToPlot)+1} = MFAorientationfieldname{tt};
        end
    end
    %     end
end

%%

fontSize = 16;
% smoothVal=201;
colorTally=1;
colorTallyHCR=1;
colorTallyMFA=1;
lineWidthPlots = 4;
f=1;
TimeScaled= TimePointsMin;
spaceMatrix = NaN(1,length(TimeScaled));
currFig = nan;
% legendTab = [];
addtoLegend = 1;
binSize = 20;
baselinePeriod = 500/binSize; %500ms for baseline

%%
for vp = 1:length(varToPlot);
    
    currFig(vp) = figure(vp+1);
    set(currFig(vp),'name',varToPlot{vp},'numbertitle','off');
    set(currFig(vp),'Position',positionGraph1)
    clf;
    %raster of licking
    clear raster_matrix raster_matrixR raster_matrixCatch raster_matrixRCatch stimStartLine stimEndLine lickLine lickLineR performance performanceCatch tempVelDATA tempVelDATACatch tempLickDATA tempLickDATACatch;
    iControl=1;
    iCatch=1;
    
    if ~isempty(data.(varToPlot{vp}).(varToAlign{f}))
        spike_matrix = data.(varToPlot{vp}).(varToAlign{f})(:,:,2)';
        spike_matrixR = data.(varToPlot{vp}).(varToAlign{f})(:,:,10)';
    else
        spike_matrix = [];
        spike_matrixR = [];
    end
    
    if (isempty(spike_matrix))
        continue
    else
        
        for i=1:size((spike_matrix),1)
            if (data.(varToPlot{vp}).(varToAlign{f})(1,i,11)==0) %if not a catch trial
                raster_matrix(iControl,:) = spike_matrix(i,:);
%                 if ~isnan(spike_matrixR)
                    raster_matrixR(iControl,:) = spike_matrixR(i,:);
%                 end
                startLine(iControl,1) = min(find(data.(varToPlot{vp}).(varToAlign{f})(:,i,5)==1));
                endLine(iControl,1) = max(find(data.(varToPlot{vp}).(varToAlign{f})(:,i,5)==1));
                stimTimes = find(data.(varToPlot{vp}).(varToAlign{f})(:,i,8)==1); %index to LED on times
                if ~isempty(stimTimes)
                    stimStartLine(iControl,1) = min(stimTimes);
                    stimEndLine(iControl,1) = max(stimTimes);
                else
                    stimStartLine(iControl,1) = nan;
                    stimEndLine(iControl,1) = nan;
                end
                
                if ~isempty(find(data.(varToPlot{vp}).(varToAlign{f})(:,i,6)==1,1))
                    lickLine(iControl,1) = min(find(data.(varToPlot{vp}).(varToAlign{f})(:,i,6)==1)); %rewardIdx
                else
                    lickLine(iControl,1) = nan;
                end
                
                if ~isempty(find(data.(varToPlot{vp}).(varToAlign{f})(:,i,9)==1,1))
                    lickLineR(iControl,1) = min(find(data.(varToPlot{vp}).(varToAlign{f})(:,i,9)==1));
                else
                    lickLineR(iControl,1) = nan;
                end
                
                if ~isempty(find(data.(varToPlot{vp}).(varToAlign{f})(:,i,12)==1,1))
                    maskLine(iControl,1) = min(find(data.(varToPlot{vp}).(varToAlign{f})(:,i,12)==1));
                else
                    maskLine(iControl,1) = nan;
                end
                
                if data.(varToPlot{vp}).(varToAlign{f})(1,i,3)==1 && data.(varToPlot{vp}).(varToAlign{f})(1,i,7)==1; % Correct Lick on GoTrials stimulus = Hit
                    performance(iControl,1) = 1;
                elseif data.(varToPlot{vp}).(varToAlign{f})(1,i,3)==0 && data.(varToPlot{vp}).(varToAlign{f})(1,i,7)==3; % No Lick on GoTrials stimulus = Miss
                    performance(iControl,1) = 3;
                elseif data.(varToPlot{vp}).(varToAlign{f})(1,i,3)==0 && data.(varToPlot{vp}).(varToAlign{f})(1,i,7)==4; % No Lick on Distractor stimulus = Correct Rejection
                    performance(iControl,1) = 4;
                elseif data.(varToPlot{vp}).(varToAlign{f})(1,i,3)==1 && data.(varToPlot{vp}).(varToAlign{f})(1,i,7)==2; % Lick on Distractor stimulus = False Alarm
                    performance(iControl,1) = 2;
                end
                iControl = iControl+1;
            else
                
                raster_matrixCatch(iCatch,:) = spike_matrix(i,:);
                raster_matrixRCatch(iCatch,:) = spike_matrixR(i,:);
                
                startLineCatch(iCatch,1) = min(find(data.(varToPlot{vp}).(varToAlign{f})(:,i,5)==1));
                endLineCatch(iCatch,1) = max(find(data.(varToPlot{vp}).(varToAlign{f})(:,i,5)==1));
                stimTimes = find(data.(varToPlot{vp}).(varToAlign{f})(:,i,8)==1);
                if ~isempty(stimTimes)
                    stimStartLineCatch(iCatch,1) = min(stimTimes);
                    stimEndLineCatch(iCatch,1) = max(stimTimes);
                else
                    stimStartLineCatch(iCatch,1) = nan;
                    stimEndLineCatch(iCatch,1) = nan;
                end
                
                if ~isempty(find(data.(varToPlot{vp}).(varToAlign{f})(:,i,6)==1,1))
                    lickLineCatch(iCatch,1) = min(find(data.(varToPlot{vp}).(varToAlign{f})(:,i,6)==1));
                else
                    lickLineCatch(iCatch,1) = nan;
                end
                
                if ~isempty(find(data.(varToPlot{vp}).(varToAlign{f})(:,i,9)==1,1))
                    lickLineRCatch(iCatch,1) = min(find(data.(varToPlot{vp}).(varToAlign{f})(:,i,9)==1));
                else
                    lickLineRCatch(iCatch,1) = nan;
                end
                
                if ~isempty(find(data.(varToPlot{vp}).(varToAlign{f})(:,i,12)==1,1))
                    maskLineCatch(iCatch,1) = min(find(data.(varToPlot{vp}).(varToAlign{f})(:,i,12)==1));
                else
                    maskLineCatch(iCatch,1) = nan;
                end
                
                if data.(varToPlot{vp}).(varToAlign{f})(1,i,3)==1 && data.(varToPlot{vp}).(varToAlign{f})(1,i,7)==1; % Correct Lick on GoTrials stimulus = Hit
                    performanceCatch(iCatch,1) = 1;
                    
                elseif data.(varToPlot{vp}).(varToAlign{f})(1,i,3)==0 && data.(varToPlot{vp}).(varToAlign{f})(1,i,7)==3; % No Lick on GoTrials stimulus = Miss
                    performanceCatch(iCatch,1) = 3;
                elseif data.(varToPlot{vp}).(varToAlign{f})(1,i,3)==0 && data.(varToPlot{vp}).(varToAlign{f})(1,i,7)==4; % No Lick on Distractor stimulus = Correct Rejection
                    performanceCatch(iCatch,1) = 4;
                elseif data.(varToPlot{vp}).(varToAlign{f})(1,i,3)==1 && data.(varToPlot{vp}).(varToAlign{f})(1,i,7)==2; % Lick on Distractor stimulus = False Alarm
                    performanceCatch(iCatch,1) = 2;
                end
                iCatch = iCatch+1;
            end
        end
        
        %only plot subset of trials if value given for noTrialsToPlot - taken
        %from middle of session unless otherwise specified
        if ~isempty(noTrialsToPlot)
            val = round(size((raster_matrix),1)/2)-(noTrialsToPlot/2);
            raster_matrix = raster_matrix(val:val+noTrialsToPlot,:);
            raster_matrixR = raster_matrixR(val:val+noTrialsToPlot,:);
            startLine = startLine(val:val+noTrialsToPlot,:);
            endLine = endLine(val:val+noTrialsToPlot,:);
            stimStartLine = stimStartLine(val:val+noTrialsToPlot,:);
            stimEndLine = stimEndLine(val:val+noTrialsToPlot,:);
            lickLine = lickLine(val:val+noTrialsToPlot,:);
            lickLineR = lickLineR(val:val+noTrialsToPlot,:);
            maskLine = maskLine(val:val+noTrialsToPlot,:);
            performance = performance(val:val+noTrialsToPlot,:);
        end
        
        if (exist('raster_matrixCatch','var'))
            raster_matrix = [raster_matrix;spaceMatrix;raster_matrixCatch];
      raster_matrixR = [raster_matrixR;spaceMatrix;raster_matrixRCatch]; %concaternate to catch trials appear at bottom
            startLine = [startLine;spaceMatrix(:,1);startLineCatch];
            endLine = [endLine;spaceMatrix(:,1);endLineCatch];
            stimStartLine = [stimStartLine;spaceMatrix(:,1);stimStartLineCatch];
            stimEndLine = [stimEndLine;spaceMatrix(:,1);stimEndLineCatch];
            lickLine = [lickLine;spaceMatrix(:,1);lickLineCatch];
            lickLineR = [lickLineR;spaceMatrix(:,1);lickLineRCatch];
            maskLine = [maskLine;spaceMatrix(:,1);maskLineCatch];
            performance = [performance;spaceMatrix(:,1);performanceCatch];
        end
        figure( currFig(vp))
        subplot(3,3,[1 2 4 5 7 8])
        %Velocity
        if (includeVelocityMap==1)
            imagesc(TimeScaled,1,data.(varToPlot{vp}).(varToAlign{f})(:,:,1)')
            hold on
        end
        
        for i=1:size((raster_matrix),1)
            raster_matrix(i,:) = raster_matrix(i,:).*i;
            raster_matrixR(i,:) = raster_matrixR(i,:).*i;
        end
        
        avgstartLine = TimeScaled(round(nanmean(startLine)));
        endLineTally = TimeScaled(round((nanmean(endLine))));
        NoTrials = size((raster_matrix),1);
        
        avgMaskLine = TimeScaled(round(nanmean(maskLine)));
        NoTrials = size((raster_matrix),1);
        
        pp = patch([avgstartLine endLineTally endLineTally avgstartLine],[0 0 NoTrials NoTrials],[.8 .8 .8],'EdgeColor','none');
        hold on
        set(pp,'FaceAlpha',0.5)
        
        pp = patch([avgMaskLine endLineTally endLineTally avgMaskLine],[0 0 NoTrials NoTrials],[.8 .8 .8],'EdgeColor','none');
        hold on
        set(pp,'FaceAlpha',0.25)
        
        plot(TimeScaled,raster_matrix,'.','MarkerEdgeColor',rgb('Black'),'MarkerSize',2);
        hold on
        plot(TimeScaled,raster_matrixR,'.','MarkerEdgeColor',rgb('Black'));
        
        
        
        for k = 1:NoTrials;
            if (isnan(lickLine(k,1))==0) %if it doesn't equal nan
                line([TimeScaled(lickLine(k,1)) TimeScaled(lickLine(k,1))], [k-0.3 k+0.3],'Color',rgb('Blue'), 'LineWidth',3);
            end
            if (isnan(lickLineR(k,1))==0) %if it doesn't equal nan
                line([TimeScaled(lickLineR(k,1)) TimeScaled(lickLineR(k,1))], [k-0.3 k+0.3],'Color',rgb('Blue'), 'LineWidth',3);
            end
            
            %                         if (isnan(maskLine(k,1))==0) %if it doesn't equal nan
            %                 line([TimeScaled(maskLine(k,1)) TimeScaled(maskLine(k,1))], [k-0.3 k+0.3],'Color',rgb('Black'), 'LineWidth',3);
            %             end
            %REAL TIMES
%             if (isnan(stimStartLine(k,1))==0) %if it doesn't equal nan
%                 line([TimeScaled(stimStartLine(k)) TimeScaled(stimStartLine(k))], [k-0.3 k+0.3],'Color',rgb('Cyan'), 'LineWidth',2);
%                 line([TimeScaled(stimEndLine(k)) TimeScaled(stimEndLine(k))], [k-0.3 k+0.3],'Color',rgb('Cyan'), 'LineWidth',2);
%             end
%             
            
            %NEATTIMES
            if (isnan(stimStartLine(k,1))==0) %if it doesn't equal nan
                line([avgMaskLine avgMaskLine], [k-0.3 k+0.3],'Color',rgb('Cyan'), 'LineWidth',2);
                line([endLineTally endLineTally], [k-0.3 k+0.3],'Color',rgb('Cyan'), 'LineWidth',2);
            end
        end
        %line where catch trials start
        CatchTrailNo = size((raster_matrix),1)-size((raster_matrixCatch),1);
        line([TimeScaled(1) TimeScaled(end)],[CatchTrailNo CatchTrailNo],'Color',rgb('Black'),'LineWidth',2)
        
        set(gca,'YDir','reverse');
        set(gca,'FontSize',16);
        xlim([min(TimeScaled), max(TimeScaled)]);
        y2 = size((raster_matrix),1);
        ylim([1 y2]);
        ylabel('Trial number', 'FontSize', fontSize,'fontWeight','bold');
        xlabel('Time (ms)', 'FontSize', fontSize,'fontWeight','bold');
        
        noCatch=1;
        Catch=1;
        for i=1:size((spike_matrix),1)
            if (data.(varToPlot{vp}).(varToAlign{f})(1,i,11)==0) %if not a catch trial
                tempVelDATA(:,noCatch) = data.(varToPlot{vp}).(varToAlign{f})(:,i,1);
                noCatch=noCatch+1;
            else
                tempVelDATACatch(:,Catch) = data.(varToPlot{vp}).(varToAlign{f})(:,i,1);
                Catch=Catch+1;
            end
        end
        
        AvgVelocity = nanmean(tempVelDATA,2);
        semVelocity = nanstd(tempVelDATA,0,2)/sqrt(size(tempVelDATA,2));
        
        if (exist('tempVelDATACatch','var'));
            AvgVelocityCatch = nanmean(tempVelDATACatch,2);
            semVelocityCatch = nanstd(tempVelDATACatch,0,2)/sqrt(size(tempVelDATACatch,2));
        end
        if ~isempty(AvgVelocity);
            [running_speed,time_bins] = JB_computeSpeed(AvgVelocity,sampleRateSpacing,TimeScaled,1,binSize);
            [running_speedSEM] = JB_computeSpeed(semVelocity,sampleRateSpacing,TimeScaled,1);
            running_speed = smooth(running_speed);
            
            if (exist('AvgVelocityCatch','var'));
                [running_speedCatch,time_bins] = JB_computeSpeed(AvgVelocityCatch,sampleRateSpacing,TimeScaled,1,binSize);
                [running_speedSEMCatch] = JB_computeSpeed(semVelocityCatch,sampleRateSpacing,TimeScaled,1);
                running_speedCatch = smooth(running_speedCatch);
            end
            
            chargeSpeed = diff(running_speed);
            changeTime = diff(time_bins)';
            acc = chargeSpeed./changeTime;
            [gg] = interp1(time_bins(1:end-1),acc,time_bins(1:end-1));
            [~,peakMin] = min(gg(10:end-10));
            [~,peakMax] = max(gg(10:end-10));
            
            gg_subset = gg(1,peakMin:end);
            time_subset = time_bins(1,peakMin:end-1);
            tallyC=1;
            clear idxIn;
            idxIn = NaN;
            for kk = 1:length(gg_subset)-1
                if ((gg_subset(kk)<0) && (gg_subset(kk+1)>0))
                    idxIn(tallyC,1) = kk;
                    tallyC = tallyC+1;
                end
            end
            %     idxIn = dsearchn(gg_subset',0);
            if ~isnan(idxIn)
                inflectionPoint = time_subset(idxIn(1));
            else
                inflectionPoint = NaN;
            end
            %smooth data and find slopes
            %low pass filter acceleration data
            smoothed_gg = smooth(gg);
            slope_gg = diff(smoothed_gg)/diff(time_bins(1:end-1)');
            %max deceleration
            [~,min_slopeIdx] = min(slope_gg(10:end-10,1));
            [~,m] = min(acc);
            [~,mA] = max(acc);
            
            %Average licking
            noCatch=1;
            Catch=1;
            for i=1:size((spike_matrix),1)
                if (data.(varToPlot{vp}).(varToAlign{f})(1,i,11)==0) %if not a catch trial
                    tempLickDATA(:,noCatch) = data.(varToPlot{vp}).(varToAlign{f})(:,i,2);
                    noCatch = noCatch+1;
                else
                    tempLickDATACatch(:,Catch) = data.(varToPlot{vp}).(varToAlign{f})(:,i,2);
                    Catch=Catch+1;
                end
            end
            
            AvgLicking = nanmean(tempLickDATA,2);
            semLicking = nanstd(tempLickDATA,0,2)/sqrt(size(tempLickDATA,2));
            
            if (exist('tempLickDATACatch','var'));
                AvgLickingCatch = nanmean(tempLickDATACatch,2);
                semLickingCatch = nanstd(tempLickDATACatch,0,2)/sqrt(size(tempLickDATACatch,2));
            end
            
            [licking_rate,time_bins] = JB_computeSpeed(AvgLicking,sampleRateSpacing,TimeScaled,0);
            [licking_rateSEM] = JB_computeSpeed(semLicking,sampleRateSpacing,TimeScaled,0);
            licking_rate = smooth(licking_rate,20);
            licking_rateSEM = smooth(licking_rateSEM,20);
            
            if (exist('AvgLickingCatch','var'));
                [licking_rateCatch,time_bins] = JB_computeSpeed(AvgLickingCatch,sampleRateSpacing,TimeScaled,0);
                [licking_rateSEMCatch] = JB_computeSpeed(semLickingCatch,sampleRateSpacing,TimeScaled,0);
                licking_rateSEMCatch = smooth(licking_rateSEMCatch,20);
                licking_rateCatch = smooth(licking_rateCatch,20);
            end
            
            subplot(3,3,[3 6 9])
            lineWidthToSet = .3;
            
            for i=1:length(performance)
                if (performance(i)==1); % Correct Lick on GoTrials stimulus = Hit
                    line([1-lineWidthToSet 1+lineWidthToSet], [i i], 'Color', rgb('Lime'),'LineWidth',2);
                    hold on;
                elseif (performance(i)==3); % No Lick on GoTrials stimulus = Miss
                    line([1-lineWidthToSet 1+lineWidthToSet], [i i], 'Color', rgb('LightGrey'),'LineWidth',2);
                    hold on;
                elseif  (performance(i)==4); % No Lick on Distractor stimulus = Correct Rejection
                    line([2-lineWidthToSet 2+lineWidthToSet], [i i], 'Color', rgb('LightGrey'),'LineWidth',2);
                    hold on;
                elseif  (performance(i)==2); % Lick on Distractor stimulus = False Alarm
                    line([2-lineWidthToSet 2+lineWidthToSet], [i i], 'Color',  rgb('Red'),'LineWidth',2);
                    hold on;
                end
            end
            
            
            ylimT = length(performance);
            if ylimT==1;
                ylimT=2;
            end
            
            ylim([1 ylimT]);
            xlim([0 3])
            line([0 3],[CatchTrailNo CatchTrailNo],'Color',rgb('Black'),'LineWidth',2)
            set(gca,'YDir','reverse');
            set(gca,'XTick',[1 2])
            xlabels = {'Left';'Right'};
            set(gca,'XTickLabel',xlabels)
        end
        
        figure(99)
        set(gcf,'Position',positionGraph3)
        subplot(2,2,1) %Plot Average Velocity
        if strcmp(varToPlot{vp},'fullSession')
            %plotCatchVelocity
            if ~isempty(AvgVelocityCatch);
                plot(time_bins,running_speedCatch, 'Color', rgb('PaleGreen'));
                hold on
                shadedErrorBar(time_bins,running_speedCatch,running_speedSEMCatch,{'-k'});
            end
            %plotsessionVelocity
            h1 = plot(time_bins,running_speed, 'Color', rgb('Lime'));
            hold on
            shadedErrorBar(time_bins,running_speed,running_speedSEM,{'-g'});
            %             plot(time_bins(m),running_speed(m),'og');
            %             plot(time_bins(mA),running_speed(mA),'ob')
            legendAdd{addtoLegend} = 'h1';
            
        elseif strcmp(varToPlot{vp},'stimTrials') || strcmp(varToPlot{vp},'NoGoTrials')
            
            hold on
            shadedErrorBar(time_bins,running_speed,running_speedSEM, {'-r'});
            h2 = plot(time_bins,running_speed, 'r','LineWidth',3);
            %             plot(time_bins(m),running_speed(m),'og');
            %             plot(time_bins(mA),running_speed(mA),'ob')
            legend(h2,varToPlot{vp})
        elseif strcmp(varToPlot{vp},'cntrTrials') || strcmp(varToPlot{vp},'GoTrials')
            
            hold on
            shadedErrorBar(time_bins,running_speed,running_speedSEM,{'-g'});
            h3 = plot(time_bins,running_speed, 'g','LineWidth',3);
            %             plot(time_bins(m),running_speed(m),'or');
            %             plot(time_bins(mA),running_speed(mA),'or')
            legend(h3,varToPlot{vp})
        end
        hold on
        set(gca,'FontSize',16);
        maxY1 = ylim;
        % ylim([-5 25]);
        % xlim([min(TimeScaled), max(TimeScaled)])
        ylabel('Velocity cm/s', 'FontSize', fontSize,'fontWeight','bold');
        xlabel('Time (ms)', 'FontSize', fontSize,'fontWeight','bold');
        hold on
        %         set(gca,'FontSize',16);
        %         ylabel('Licks per s', 'FontSize', fontSize,'fontWeight','bold');
        %         xlabel('Time (ms)', 'FontSize', fontSize,'fontWeight','bold');
        markerSize = 10;
        
        
        subplot(2,2,3) %Plot Average Velocity
        if strcmp(varToPlot{vp},'stimTrials') || strcmp(varToPlot{vp},'NoGoTrials')
            
            hold on
            shadedErrorBar(time_bins,running_speedCatch,running_speedSEMCatch, {'-r'});
            h2 = plot(time_bins,running_speedCatch, 'r','LineWidth',3);
            legend(h2,varToPlot{vp})
        elseif strcmp(varToPlot{vp},'cntrTrials') || strcmp(varToPlot{vp},'GoTrials')
            
            hold on
            shadedErrorBar(time_bins,running_speedCatch,running_speedSEMCatch,{'-g'});
            h3 = plot(time_bins,running_speedCatch, 'g','LineWidth',3);
            %             plot(time_bins(m),running_speed(m),'or');
            %             plot(time_bins(mA),running_speed(mA),'or')
            legend(h3,varToPlot{vp})
        end
        hold on
        set(gca,'FontSize',16);
        maxY2 = ylim;
        
        
        if maxY1(1)<maxY2(1)
            yToUse(1) = maxY1(1);
        else
            yToUse(1) = maxY2(1);
        end
        
        if maxY1(2)>maxY2(2)
            yToUse(2) = maxY1(2);
        else
            yToUse(2) = maxY2(2);
        end
        
        subplot(2,2,1)
        ylim([yToUse(1) yToUse(2)])
        subplot(2,2,3)
        ylim([yToUse(1) yToUse(2)])
        
        
        
        % ylim([-5 25]);
        % xlim([min(TimeScaled), max(TimeScaled)])
        ylabel('Velocity cm/s CATCH', 'FontSize', fontSize,'fontWeight','bold');
        xlabel('Time (ms)', 'FontSize', fontSize,'fontWeight','bold');
        hold on
        %         set(gca,'FontSize',16);
        %         ylabel('Licks per s', 'FontSize', fontSize,'fontWeight','bold');
        %         xlabel('Time (ms)', 'FontSize', fontSize,'fontWeight','bold');
        markerSize = 10;
        
        
        
        
        subplot(2,2,2) %Average licking
        if strcmp(varToPlot{vp},'fullSession')
            if ~isempty(AvgLickingCatch);
                shadedErrorBar(time_bins,licking_rateCatch,licking_rateSEMCatch,{'-b'})
                hold on
            end
            shadedErrorBar(time_bins,licking_rate,licking_rateSEM,{'-k'})
            hold on
            
        elseif strcmp(varToPlot{vp},'stimTrials') || strcmp(varToPlot{vp},'NoGoTrials')
            plot(time_bins,licking_rate, 'r','LineWidth',3);
            %             shadedErrorBar(time_bins,licking_rate,licking_rateSEM,{'-b'})
            hold on
        elseif strcmp(varToPlot{vp},'cntrTrials')  || strcmp(varToPlot{vp},'GoTrials')
            plot(time_bins,licking_rate, 'g','LineWidth',3);
            %             shadedErrorBar(time_bins,licking_rate,licking_rateSEM,{'-k'})
            hold on
        end
        hold on
        set(gca,'FontSize',16);
        ylabel('Licks per s', 'FontSize', fontSize,'fontWeight','bold');
        xlabel('Time (ms)', 'FontSize', fontSize,'fontWeight','bold');
        markerSize = 10;
        maxY1 = ylim;
        
        subplot(2,2,4) %Average licking Catch
        if strcmp(varToPlot{vp},'stimTrials') || strcmp(varToPlot{vp},'NoGoTrials')
            plot(time_bins,licking_rateCatch, 'r','LineWidth',3)
            hold on
            %             shadedErrorBar(time_bins,licking_rateCatch,licking_rateSEMCatch,{'-b'})
        elseif strcmp(varToPlot{vp},'cntrTrials')  || strcmp(varToPlot{vp},'GoTrials')
            plot(time_bins,licking_rateCatch, 'g','LineWidth',3)
            hold on
            %              shadedErrorBar(time_bins,licking_rateCatch,licking_rateSEMCatch,{'-k'})
        end
        
        maxY2 = ylim;
        
        
        if maxY1(1)<maxY2(1)
            yToUse(1) = maxY1(1);
        else
            yToUse(1) = maxY2(1);
        end
        
        if maxY1(2)>maxY2(2)
            yToUse(2) = maxY1(2);
        else
            yToUse(2) = maxY2(2);
        end
        
        subplot(2,2,2)
        ylim([yToUse(1) yToUse(2)])
        subplot(2,2,4)
        ylim([yToUse(1) yToUse(2)])
        
        
        hold on
        set(gca,'FontSize',16);
        ylabel('Licks per s CATCH', 'FontSize', fontSize,'fontWeight','bold');
        xlabel('Time (ms)', 'FontSize', fontSize,'fontWeight','bold');
        markerSize = 10;
        
        if (plotAngles==1)
            figure(100);
            set(gcf,'Position',positionGraph4)
            set(gcf,'name','allDATA','numbertitle','off');
            subplot(2,1,1) %Plot Average Velocity
            if strncmp(varToPlot{vp},'stimAngle',9)
                plot(time_bins,running_speed,'k','LineWidth',lineWidthPlots,'Color',  rgb(plotColor{colorTally}))
                hold on
                plot([inflectionPoint inflectionPoint],[min(ylim) max(ylim)],'--','LineWidth',2,'Color', rgb(plotColor{colorTally}))
                plot(time_bins(min_slopeIdx),running_speed(min_slopeIdx),'o','Color', rgb(plotColor{colorTally}),'MarkerSize',markerSize)
                plot(time_bins(peakMin),running_speed(peakMin),'+','Color', rgb(plotColor{colorTally}),'MarkerSize',markerSize)
                plot(time_bins(peakMax),running_speed(peakMax),'*','Color', rgb(plotColor{colorTally}),'MarkerSize',markerSize)
                %                         plot(time_bins(m),running_speed(m),'og');
                %                         plot(time_bins(mA),running_speed(mA),'ob')
                set(gca,'FontSize',16);
                ylabel('Velocity cm/s', 'FontSize', fontSize,'fontWeight','bold');
                xlabel('Time (ms)', 'FontSize', fontSize,'fontWeight','bold');
                h_legend = legend(orientationfieldname)
                set(h_legend,'FontSize',6)
                %ylim([0 20]);
                
                subplot(2,1,2) %Average licking
                plot(time_bins,licking_rate,'k','LineWidth',lineWidthPlots,'Color', rgb(plotColor{colorTally}))
                hold on
                %             plot([inflectionPoint inflectionPoint],[min(ylim) max(ylim)],'b--','LineWidth',2)
                %             plot(time_bins(min_slopeIdx),running_speed(min_slopeIdx),'*m')
                %             plot(time_bins(peakMin),running_speed(peakMin),'*g')
                %             plot(time_bins(peakMax),running_speed(peakMax),'*r')
                %                         plot(time_bins(m),running_speed(m),'og');
                %                         plot(time_bins(mA),running_speed(mA),'ob')
                colorTally = colorTally+1;
                hold on
                set(gca,'FontSize',16);
                ylabel('Licks per s', 'FontSize', fontSize,'fontWeight','bold');
                xlabel('Time (ms)', 'FontSize', fontSize,'fontWeight','bold');
                h_legend = legend(orientationfieldname)
                set(h_legend,'FontSize',6)
            end
        end
        if (plotAngles==1)
            figure(200);
            
            set(gcf,'Position',positionGraph4)
            set(gcf,'name','Hit&CR','numbertitle','off');
            subplot(2,1,1) %Plot Average Velocity
            if strncmp(varToPlot{vp},'HCRstimAngle',9)
                plot(time_bins,running_speed,'k','LineWidth',lineWidthPlots,'Color',  rgb(plotColor{colorTallyHCR}))
                hold on
                plot([inflectionPoint inflectionPoint],[min(ylim) max(ylim)],'--','LineWidth',2,'Color', rgb(plotColor{colorTallyHCR}))
                plot(time_bins(min_slopeIdx),running_speed(min_slopeIdx),'o','Color', rgb(plotColor{colorTallyHCR}),'MarkerSize',markerSize)
                plot(time_bins(peakMin),running_speed(peakMin),'+','Color', rgb(plotColor{colorTallyHCR}),'MarkerSize',markerSize)
                plot(time_bins(peakMax),running_speed(peakMax),'*','Color', rgb(plotColor{colorTallyHCR}),'MarkerSize',markerSize)
                % plot(time_bins(m),running_speed(m),'og');
                % plot(time_bins(mA),running_speed(mA),'ob')
                set(gca,'FontSize',16);
                ylabel('Velocity cm/s', 'FontSize', fontSize,'fontWeight','bold');
                xlabel('Time (ms)', 'FontSize', fontSize,'fontWeight','bold');
                h_legend = legend(orientationfieldname)
                set(h_legend,'FontSize',6)
                %ylim([0 20]);
                
                subplot(2,1,2) %Average licking
                plot(time_bins,licking_rate,'k','LineWidth',lineWidthPlots,'Color', rgb(plotColor{colorTallyHCR}))
                hold on
                %             plot([inflectionPoint inflectionPoint],[min(ylim) max(ylim)],'--','LineWidth',2,'Color', rgb(plotColor{colorTallyHCR}))
                %             plot(time_bins(min_slopeIdx),running_speed(min_slopeIdx),'o','Color', rgb(plotColor{colorTallyHCR}))
                %             plot(time_bins(peakMin),running_speed(peakMin),'+','Color', rgb(plotColor{colorTallyHCR}))
                %             plot(time_bins(peakMax),running_speed(peakMax),'*','Color', rgb(plotColor{colorTallyHCR}))
                % plot(time_bins(m),running_speed(m),'og');
                % plot(time_bins(mA),running_speed(mA),'ob')
                colorTallyHCR = colorTallyHCR+1;
                hold on
                set(gca,'FontSize',16);
                ylabel('Licks per s', 'FontSize', fontSize,'fontWeight','bold');
                xlabel('Time (ms)', 'FontSize', fontSize,'fontWeight','bold');
                h_legend = legend(orientationfieldname)
                set(h_legend,'FontSize',6)
            end
        end
        if (plotAngles==1)
            figure(201);
            
            set(gcf,'Position',positionGraph4)
            set(gcf,'name','Miss&FA','numbertitle','off');
            subplot(2,1,1) %Plot Average Velocity
            if strncmp(varToPlot{vp},'MFAstimAngle',9)
                plot(time_bins,running_speed,'k','LineWidth',lineWidthPlots,'Color',  rgb(plotColor{colorTallyMFA}))
                hold on
                plot([inflectionPoint inflectionPoint],[min(ylim) max(ylim)],'--','LineWidth',2,'Color', rgb(plotColor{colorTallyMFA}))
                plot(time_bins(min_slopeIdx),running_speed(min_slopeIdx),'o','Color', rgb(plotColor{colorTallyMFA}),'MarkerSize',markerSize)
                plot(time_bins(peakMin),running_speed(peakMin),'+','Color', rgb(plotColor{colorTallyMFA}),'MarkerSize',markerSize)
                plot(time_bins(peakMax),running_speed(peakMax),'*','Color', rgb(plotColor{colorTallyMFA}),'MarkerSize',markerSize)
                % plot(time_bins(m),running_speed(m),'og');
                % plot(time_bins(mA),running_speed(mA),'ob')
                set(gca,'FontSize',16);
                ylabel('Velocity cm/s', 'FontSize', fontSize,'fontWeight','bold');
                xlabel('Time (ms)', 'FontSize', fontSize,'fontWeight','bold');
                h_legend = legend(orientationfieldname)
                set(h_legend,'FontSize',6)
                
                %ylim([0 20]);
                
                subplot(2,1,2) %Average licking
                plot(time_bins,licking_rate,'k','LineWidth',lineWidthPlots,'Color', rgb(plotColor{colorTallyMFA}))
                hold on
                %             plot([inflectionPoint inflectionPoint],[min(ylim) max(ylim)],'b--','LineWidth',2)
                %             plot(time_bins(min_slopeIdx),running_speed(min_slopeIdx),'*m')
                %             plot(time_bins(peakMin),running_speed(peakMin),'*g')
                %             plot(time_bins(peakMax),running_speed(peakMax),'*r')
                %             plot(time_bins(m),running_speed(m),'og');
                %             plot(time_bins(mA),running_speed(mA),'ob')
                colorTallyMFA = colorTallyMFA+1;
                hold on
                set(gca,'FontSize',16);
                ylabel('Licks per s', 'FontSize', fontSize,'fontWeight','bold');
                xlabel('Time (ms)', 'FontSize', fontSize,'fontWeight','bold');
                h_legend = legend(orientationfieldname)
                set(h_legend,'FontSize',6)
            end
        end
        
        figure(101)
        subplot(1,2,1)
        tempangles = data.angleOrientations;
        latencytemp = nan(100,length(tempangles));
        latencytempStim = nan(100,length(tempangles));
        latencytempCatch = nan(100,length(tempangles));
        latencytempCatchStim = nan(100,length(tempangles));
        
        for gg = 1:length(tempangles)
            idx =  find(data.CumulativePerformance(:,1)==(tempangles(gg)));
            for ggh = 1:length(idx)
                latencyFromMasking = data.CumulativePerformance(idx(ggh),9);
                if  (data.CumulativePerformance(idx(ggh),8)==0); %non Catch Trials
                    if  (data.CumulativePerformance(idx(ggh),6)==0); %non stim Trials
                        latencytemp(ggh,gg) = latencyFromMasking; %latency from Masking Light
                    elseif (data.CumulativePerformance(idx(ggh),6)==1) %stim trial
                        latencytempStim(ggh,gg) = latencyFromMasking;
                    end
                elseif (data.CumulativePerformance(idx(ggh),8)==1); %Catch Trials
                    if  (data.CumulativePerformance(idx(ggh),6)==0); %non stim Trials
                        latencytempCatch(ggh,gg) = latencyFromMasking; %latency from Masking Light
                    elseif (data.CumulativePerformance(idx(ggh),6)==1) %latency from stim
                        latencytempCatchStim(ggh,gg) = latencyFromMasking;
                    end
                end
            end
        end
        
        if (~isempty(latencytemp))
            [meanData,~,~,semData] = JB_calBasicStats(latencytemp);
            latencyAVG(1,1) =  nanmean(meanData);
            latencyAVG(1,2) = nanmean(semData);
            plot(tempangles,meanData,'o-k','MarkerSize',3,'LineWidth',3)
            hold on
            errorbar(tempangles,meanData,semData,'k','Linestyle','none');
        end
        
        if (~isempty(latencytempStim))
            [meanData,~,~,semData] = JB_calBasicStats(latencytempStim);
            latencyAVG(2,1) =  nanmean(meanData);
            latencyAVG(2,2) = nanmean(semData);
            plot(tempangles,meanData,'o-r','MarkerSize',3,'LineWidth',3)
            hold on
            errorbar(tempangles,meanData,semData,'r','Linestyle','none');
        end
        
        if (~isempty(latencytempCatch))
            [meanData,~,~,semData] = JB_calBasicStats(latencytempCatch);
            latencyAVG(3,1) =  nanmean(meanData);
            latencyAVG(3,2) = nanmean(semData);
            %         plot(tempangles,meanData,'o-b','MarkerSize',3,'LineWidth',3)
            %         hold on
            %         errorbar(tempangles,meanData,semData,'k','Linestyle','none');
        end
        
        if (~isempty(latencytempCatchStim))
            [meanData,~,~,semData] = JB_calBasicStats(latencytempCatchStim);
            latencyAVG(4,1) =  nanmean(meanData);
            latencyAVG(4,2) = nanmean(semData);
            %         plot(tempangles,meanData,'o-g','MarkerSize',3,'LineWidth',3)
            %         hold on
            %         errorbar(tempangles,meanData,semData,'k','Linestyle','none');
        end
        
        
        NumTicks = length(tempangles);
        L = [tempangles(1) tempangles(end)];
        % L = get(gca,'XLim');
        labels = num2str(tempangles-270);
        set(gca,'XTick',tempangles,'XTickLabel',labels,'XAxisLocation','top','ylim',[0 2],'Ydir','reverse')
        view(-90,90)
        xlabel('angles');
        ylabel('LickLatency from masking Light onset');
        
        subplot(1,2,2)
        plot(latencyAVG(:,1),'ob','MarkerFaceColor','b')
        hold on
        errorbar(latencyAVG(:,1),latencyAVG(:,2),'ko-','MarkerFaceColor','k','LineWidth',2)
%         labelsBar = {'go';'goCatch';'NoGo';'NoGoCatch'};
         labelsBar = {'control';'stim';'controlCatch';'stimCatch'};
        set(gca,'xticklabels',labelsBar)
        set(gca,'XTick',[1 2 3 4],'XTickLabel',labelsBar,'XAxisLocation','top','ylim',[0 3],'Ydir','reverse')
        view(-90,90)
        xlabel('angles');
        ylabel('LickLatency from masking Light onset');
        
        
        %plot latency to first lick from stim on (or estimated stim on)
        %         if   (stim==1);
        %             subplot(1,2,2)
        %             tempangles = data.angleOrientations;
        %             latencytemp = nan(100,length(tempangles));
        %             for gg = 1:length(tempangles)
        %                 tempCountSTIM = 1;
        %                 tempCountNoSTIM = 1;
        %                 %for Stim trials
        %                 for kk = 1:length(data.CumulativePerformance)
        %                     if (data.CumulativePerformance(kk,1)==(tempangles(gg)) && (data.CumulativePerformance(kk,6)==1));
        %
        %                         if (data.CumulativePerformance(kk,8)==0) %exclude Catch Trials
        %                             latencytempSTIM(tempCountSTIM,gg) = data.CumulativePerformance(kk,7);
        %                             tempCountSTIM = tempCountSTIM+1;
        %                         end
        %                     elseif (data.CumulativePerformance(kk,1)==(tempangles(gg)) && (data.CumulativePerformance(kk,6)==0));
        %                         if (data.CumulativePerformance(kk,8)==0) %exclude Catch Trials
        %                             latencytempNoSTIM(tempCountNoSTIM,gg) = data.CumulativePerformance(kk,7);
        %                             tempCountNoSTIM = tempCountNoSTIM+1;
        %                         end
        %                     end
        %                 end
        %             end
        %
        %             [meanDataSTIM,~,~,semDataSTIM] = JB_calBasicStats(latencytempSTIM);
        %             [meanDataNoSTIM,~,~,semDataNoSTIM] = JB_calBasicStats(latencytempNoSTIM);
        %
        %             plot(tempangles,meanDataNoSTIM,'o-k','MarkerSize',3,'LineWidth',3)
        %             hold on
        %             errorbar(tempangles,meanDataNoSTIM,semDataNoSTIM,'k','Linestyle','none');
        %
        %             plot(tempangles,meanDataSTIM,'o-r','MarkerSize',3,'LineWidth',3)
        %             errorbar(tempangles,meanDataSTIM,semDataSTIM,'r','Linestyle','none');
        %             NumTicks = length(tempangles);
        %             L = [tempangles(1) tempangles(end)];
        %             % L = get(gca,'XLim');
        %             labels = num2str(tempangles-270);
        %             set(gca,'XTick',tempangles,'XTickLabel',labels,'XAxisLocation','top','ylim',[0 3],'Ydir','reverse')
        %             view(-90,90)
        %             xlabel('angles');
        %             ylabel('LickLatency from stim Onset');
        %
        %
        %             %     figure(103)
        %             %     tempangles = data.angleOrientations;
        %             %     accelerationTemp = nan(100,length(tempangles));
        %             %     for gg = 1:length(tempangles)
        %             %         %only for NoGo trials
        %             %         idx =  find(stimType(:,2)==(tempangles(gg)));
        %             %         for ggh = 1:length(idx)
        %             %             accelerationTemp(ggh,gg) = accelerationTime(idx(ggh),2);
        %             %         end
        %             %     end
        %             %
        %             %     meanLatency = nanmean(accelerationTemp);
        %             %     stdLatency = nanstd(accelerationTemp);
        %             %     plot(tempangles,meanLatency,'o-k','MarkerSize',3,'LineWidth',3)
        %             %     hold on
        %             %     errorbar(tempangles,meanLatency,stdLatency,'k','Linestyle','none');
        %             %     NumTicks = length(tempangles);
        %             %     L = [tempangles(1) tempangles(end)];
        %             %     % L = get(gca,'XLim');
        %             %     labels = num2str(tempangles-270);
        %             %     %         set(gca,'XTick',tempangles,'XTickLabel',labels,'XAxisLocation','top','ylim',[0 1],'Ydir','reverse')
        %             %
        %             %     set(gca,'XTick',tempangles,'XTickLabel',labels,'XAxisLocation','top','Ydir','reverse')
        %             %     view(-90,90)
        %             %     xlabel('angles');
        %             %     ylabel('accelerationTemp');
        %
        %         end
        
        % if (length(varToPlot)>1)
        %     figure(102)
        %     bar(nanmean(tempAccAll))
        %     set(gca,'XTick',[1:2]);
        %     set(gca,'XTickLabel',(varToPlot));
        %     ylabel('Latency to acceleration')
        %
        % end
        
        figure(21)
        set(gcf,'Position',positionGraph2)
        subplot(2,1,1);
        plot(cumulativePerformance,'*b');
        ylabel('cumulative performance', 'FontSize', fontSize,'fontWeight','bold');
        xlabel('Trial number', 'FontSize', fontSize,'fontWeight','bold');
        ylim([0,1])
        
    end
    
    %     subplot(2,1,2);
    %     plot(blockPerformancePlot,'-b','LineWidth',2);
    %     ylabel('block performance', 'FontSize', fontSize,'fontWeight','bold');
    %     xlabel('block (12) number', 'FontSize', fontSize,'fontWeight','bold');
    %     ylim([0,1]);
end




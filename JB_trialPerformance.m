function [data] = JB_trialPerformance(DATA,sessionNo,stim,plotAngles,includeVelocityMap)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% load DATA.mat
positionGraph2 = [680 288 314 690];
positionGraph1 = [1501 71 384 902];
positionGraph3 = [2 42 958 954];
% plotColor = {'k','b','r','g','m','y',[.5 .6 .7],[.8 .2 .6]};
plotColor = {'Red','Crimson','Salmon','Pink','LightSteelBlue','LightSkyBlue','DodgerBlue','MediumBlue'};
% plotColor = {'DarkRed','Red','Salmon','LightSalmon','LightBlue','MediumTurquoise','DodgerBlue','Blue'};
PreStim=2000; %time stamps are in 10ms bins, so 100 = 1 second 200 = 2seconds
PostStim=2000;
figure(99);clf;

tempDATA = DATA.allFiles{1,sessionNo}.rawData;

if (stim==0)
    varToPlot = {'fullSession'};
elseif (stim==1)
    varToPlot = {'stimTrials' 'cntrTrials'};
elseif (stim==2)
    varToPlot = {'GoTrials' 'NoGoTrials'};
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
angle = 11;
trialType = 13;

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

C = unique(tempDATA(:,angle));
data.angleOrientations = C(C>0);
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

%% Define Target Time, Position, Correct (Hit, correct Rejection) or Incorrect (Miss, or False Alarm) data.trial

sampleRateSpacing = 1; %2ms sampling
hitCount = 0;
crCount = 0;

for j = 1: length(data.trial);
    
    % make sample times uniform and zeroed
    SampleTimes = data.trial{j}.raw(:,millis)-data.trial{j}.raw(1,millis);
    clear tempFiltered;
    
    %make each data point uniformly sampled
    for i = 1:size((data.trial{j}.raw),2);
        tempTrial = data.trial{j}.raw(:,i);
        filteredSampleTimes = 0:sampleRateSpacing:max(SampleTimes); %2msec spacing
        tempFiltered(:,i)=interp1(SampleTimes,tempTrial,filteredSampleTimes);
    end
    
    % make sure all values are integers
    data.trial{j}.filteredDATA = [tempFiltered(:,1) (round(tempFiltered(:,2:end)))];
    
    % make sure licks are represented as just one tick mark
    B = diff(data.trial{j}.filteredDATA(:,licks))==1; %%make sure licks are represented as just one tick mark
    B(numel(data.trial{j}.filteredDATA(:,licks))) = 0;
    data.trial{j}.filteredDATA(:,licks)=B;
    
    data.trial{j}.timeStamps(:,1) = [0;(cumsum(diff(data.trial{j}.filteredDATA(:,millis)))/1000)];
    data.trial{j}.Performance.StimType = data.trial{j}.filteredDATA(1,angle);
    data.trial{j}.Performance.TrialType = data.trial{j}.raw(find(data.trial{j}.raw(:,trialType)>0),trialType);
    
    [loc2,~] = find((data.trial{j}.filteredDATA(:,targetArea)>0)); % Identify when stim entered TargetArea
    rewardIdx = find((data.trial{j}.filteredDATA(:,reward)==1)); % Identify RewardTime if Hit data.trial
    
    % If Reward was given, update information on record properties
    if ~isempty(rewardIdx); %
        data.trial{j}.Stimulus.RewardIdx = rewardIdx(1);
    else
        data.trial{j}.Stimulus.RewardIdx = nan;
    end
    
    %Gather information about Answer Period
    data.trial{j}.Stimulus.AnswerPeriodDuration = (loc2(end)-loc2(1))*sampleRateSpacing;
    data.trial{j}.Stimulus.AnswerPeriodStartIdx = loc2(1);
    data.trial{j}.Stimulus.AnswerPeriodEndIdx = loc2(end);
    
    % Determine if animal licked during Answer Period
    
    %Look at time of first lick for all stim and all trial types cover both
    %the answer period and 1 second after this
    SampleLickTime = 2000/sampleRateSpacing;  % = 2seconds
    EndTime = data.trial{j}.Stimulus.AnswerPeriodStartIdx+SampleLickTime;
    if EndTime>length(data.trial{j}.filteredDATA)
        EndTime = length(data.trial{j}.filteredDATA);
    end
    
    %if reward port opened find licks before this happened
    if ~isnan(data.trial{j}.Stimulus.RewardIdx) %If there was a reward given
        data.trial{j}.Performance.LickFrequency = sum(data.trial{j}.filteredDATA(data.trial{j}.Stimulus.AnswerPeriodStartIdx:data.trial{j}.Stimulus.RewardIdx,licks)==1)/(data.trial{j}.Stimulus.AnswerPeriodDuration/1000);
        tempFirstLickIdx = (data.trial{j}.Stimulus.AnswerPeriodStartIdx+(find((data.trial{j}.filteredDATA(data.trial{j}.Stimulus.AnswerPeriodStartIdx:data.trial{j}.Stimulus.RewardIdx,licks)==1),1)));
        tempLickLatency = data.trial{j}.timeStamps(tempFirstLickIdx)-data.trial{j}.timeStamps(data.trial{j}.Stimulus.AnswerPeriodStartIdx);
        if ~isempty(tempLickLatency)
            data.trial{j}.Performance.FirstLickLatency = tempLickLatency;
        else
            data.trial{j}.Performance.FirstLickLatency = nan;
        end
    else
        if length(data.trial{j}.filteredDATA)>data.trial{j}.Stimulus.AnswerPeriodEndIdx;%if reward port didn't open find out if any licks occured during whloe answer period
            maxTimeIdx = data.trial{j}.Stimulus.AnswerPeriodEndIdx;
        else
            maxTimeIdx = length(data.trial{j}.filteredDATA);
        end
        data.trial{j}.Performance.LickFrequency = sum(data.trial{j}.filteredDATA(data.trial{j}.Stimulus.AnswerPeriodStartIdx:maxTimeIdx,licks)==1)/(data.trial{j}.Stimulus.AnswerPeriodDuration/1000);
        tempFirstLickIdx = (data.trial{j}.Stimulus.AnswerPeriodStartIdx+(find((data.trial{j}.filteredDATA(data.trial{j}.Stimulus.AnswerPeriodStartIdx:maxTimeIdx,licks)==1),1)));
        tempLickLatency = data.trial{j}.timeStamps(tempFirstLickIdx)-data.trial{j}.timeStamps(data.trial{j}.Stimulus.AnswerPeriodStartIdx);
        if ~isempty(tempLickLatency)
            data.trial{j}.Performance.FirstLickLatency = tempLickLatency;
        else
            data.trial{j}.Performance.FirstLickLatency = nan;
        end
    end
    
    if any(data.trial{j}.raw(:,LEDstate),1);
        data.trial{j}.Performance.StimTrial = 1;
    else
        data.trial{j}.Performance.StimTrial = 0;
    end
    
    % Update Trial Type Results
    if data.trial{j}.Performance.LickFrequency>0; % If there was licking..
        if (data.trial{j}.Performance.StimType < midPoint); % If this was the GoTrials stimulus - Hit trial
            hitCount = hitCount+1;
            data.trial{j}.Performance.OutcomeLick = 1; % Anticipatory lick - Hit trial
            data.trial{j}.Performance.OutcomeCode = 1;
        elseif (data.trial{j}.Performance.StimType > midPoint) % If this was a distractor stimulus
            data.trial{j}.Performance.OutcomeLick = 1;% False Alarm - incorrect lick
            data.trial{j}.Performance.OutcomeCode = 2;
        elseif (data.trial{j}.Performance.StimType == midPoint) % If this was exactly in the middle
            data.trial{j}.Performance.OutcomeLick = 1;% licked for 50% verticle
            data.trial{j}.Performance.OutcomeCode = 5; % 5 if licked, 6 if didn't
        end
    else
        if (data.trial{j}.Performance.StimType < midPoint); % If this was the no GoTrials stimulus
            data.trial{j}.Performance.OutcomeLick = 0;  % No anticipatory lick - Miss trial
            data.trial{j}.Performance.OutcomeCode = 3;
        elseif (data.trial{j}.Performance.StimType > midPoint) % If this was a distractor stimulus
            crCount = crCount+1;
            data.trial{j}.Performance.OutcomeLick = 0; % Correctly withheld lick
            data.trial{j}.Performance.OutcomeCode = 4;
        elseif (data.trial{j}.Performance.StimType == midPoint) % If this was exactly in the middle
            data.trial{j}.Performance.OutcomeLick = 0;% licked for 50% verticle
            data.trial{j}.Performance.OutcomeCode = 6; % 5 if licked, 6 if didn't
        end
    end
    cumulativePerformance(j,1) = ((hitCount+crCount)/j);
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
    
    for j=1:length(data.trial);
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
        tempVelocityStim = velocityTemp(data.trial{j}.Stimulus.(varToAlign{f})-currPreStim:data.trial{j}.Stimulus.(varToAlign{f})+currPostStim-1,1)';
        tempLickStim = (data.trial{j}.filteredDATA(data.trial{j}.Stimulus.(varToAlign{f})-currPreStim:data.trial{j}.Stimulus.(varToAlign{f})+currPostStim-1,4))';
        IdxTarget = (data.trial{j}.filteredDATA(data.trial{j}.Stimulus.(varToAlign{f})-currPreStim:data.trial{j}.Stimulus.(varToAlign{f})+currPostStim-1,2))';
        IdxReward = (data.trial{j}.filteredDATA(data.trial{j}.Stimulus.(varToAlign{f})-currPreStim:data.trial{j}.Stimulus.(varToAlign{f})+currPostStim-1,3))';
        
        data.fullSession.(varToAlign{f})(1:fullSampleTime,j,1) = [addPreStim tempVelocityStim addPostStim];
        data.fullSession.(varToAlign{f})(1:fullSampleTime,j,2) = [addPreStim tempLickStim addPostStim];
        data.fullSession.(varToAlign{f})(1,j,3) = data.trial{j}.Performance.OutcomeLick;
        data.fullSession.(varToAlign{f})(1,j,4) = data.trial{j}.Performance.StimType;
        data.fullSession.(varToAlign{f})(1:fullSampleTime,j,5) = [addPreStim IdxTarget addPostStim];
        data.fullSession.(varToAlign{f})(1:fullSampleTime,j,6) = [addPreStim IdxReward addPostStim];
        data.fullSession.(varToAlign{f})(1,j,7) = data.trial{j}.Performance.OutcomeCode;
        
        
    end
    
    for j=1:length(data.trial);
        
        %Separate Data output into GoTrials Stim and Distractor and Hit Miss Correct Rejection False Alarm
        %Trials
        
        %           if (data.trial{j}.Performance.StimTrial ==1); %Stimulated Session
        %             data.stimTrials.(varToAlign{f})(:,size(data.stimTrials.(varToAlign{f}),2),:) = data.fullSession.(varToAlign{f})(:,j,:);
        %         else
        %             data.cntrTrials.(varToAlign{f})(:,size(data.cntrTrials.(varToAlign{f}),2),:) = data.fullSession.(varToAlign{f})(:,j,:);
        %         end
        %
        %         if (data.trial{j}.Performance.StimType < midPoint); %GoTrials Stim
        %             data.GoTrials.(varToAlign{f})(:,j,:) = data.fullSession.(varToAlign{f})(:,j,:);
        %             if data.trial{j}.Performance.OutcomeLick==1; %HIT trial
        %                 data.AlignedDATA.TrialType.HIT.(varToAlign{f})(:,j,:) = data.fullSession.(varToAlign{f})(:,j,:);
        %                 if data.trial{j}.Performance.StimTrial==1; %Stimulation
        %                     data.AlignedDATA.TrialType.HITStim.(varToAlign{f})(:,j,:) = data.fullSession.(varToAlign{f})(:,j,:);
        %                 else
        %                     data.AlignedDATA.TrialType.HITNoStim.(varToAlign{f})(:,j,:) = data.fullSession.(varToAlign{f})(:,j,:);
        %                 end
        %             elseif data.trial{j}.Performance.OutcomeLick==0 %Miss trial
        %                 data.AlignedDATA.TrialType.MISS.(varToAlign{f})(:,j,:) = data.fullSession.(varToAlign{f})(:,j,:);
        %                 if data.trial{j}.Performance.StimTrial==1; %Stimulation
        %                     data.AlignedDATA.TrialType.MISSStim.(varToAlign{f})(:,j,:) = data.fullSession.(varToAlign{f})(:,j,:);
        %                 else
        %                     data.AlignedDATA.TrialType.MISSNoStim.(varToAlign{f})(:,j,:) = data.fullSession.(varToAlign{f})(:,j,:);
        %                 end
        %             end
        %         elseif (data.trial{j}.Performance.StimType > midPoint); %NoGoTrials Stim
        %             data.NoGoTrials.(varToAlign{f})(:,j,:) = data.fullSession.(varToAlign{f})(:,j,:);
        %             if data.trial{j}.Performance.OutcomeLick==0; %Correct Rejection Trial
        %                 data.AlignedDATA.TrialType.CorrectRejection.(varToAlign{f})(:,j,:) = data.fullSession.(varToAlign{f})(:,j,:);
        %                 if data.trial{j}.Performance.StimTrial==1; %Stimulation
        %                     data.AlignedDATA.TrialType.CorrectRejectionStim.(varToAlign{f})(:,j,:) = data.fullSession.(varToAlign{f})(:,j,:);
        %                 else
        %                     data.AlignedDATA.TrialType.CorrectRejectionNoStim.(varToAlign{f})(:,j,:) = data.fullSession.(varToAlign{f})(:,j,:);
        %                 end
        %             elseif data.trial{j}.Performance.OutcomeLick==1; %false alarm trial
        %                 data.AlignedDATA.TrialType.FalseAlarm.(varToAlign{f})(:,j,:) = data.fullSession.(varToAlign{f})(:,j,:);
        %                 if data.trial{j}.Performance.StimTrial==1; %Stimulation
        %                     data.AlignedDATA.TrialType.FalseAlarmStim.(varToAlign{f})(:,j,:) = data.fullSession.(varToAlign{f})(:,j,:);
        %                 else
        %                     data.AlignedDATA.TrialType.FalseAlarmNoStim.(varToAlign{f})(:,j,:) = data.fullSession.(varToAlign{f})(:,j,:);
        %                 end
        %             end
        %         elseif (data.trial{j}.Performance.StimType == midPoint); %midpoint
        %             data.AlignedDATAallDATAmidPoint.(varToAlign{f})(:,j,:) = data.fullSession.(varToAlign{f})(:,j,:);
        %             if data.trial{j}.Performance.OutcomeLick==0; %didn't lick
        %                 data.AlignedDATA.TrialType.midPointNoLickCount.(varToAlign{f})(:,j,:) = data.fullSession.(varToAlign{f})(:,j,:);
        %             elseif data.trial{j}.Performance.OutcomeLick==1; %false alarm trial
        %                 data.AlignedDATA.TrialType.midPointLickCount.(varToAlign{f})(:,j,:) = data.fullSession.(varToAlign{f})(:,j,:);
        %             end
        %         end
        
        %Here
        if (data.trial{j}.Performance.StimTrial ==1); %Stimulated Session
            data.stimTrials.(varToAlign{f})(:,stimON,:) = data.fullSession.(varToAlign{f})(:,j,:);
            stimON = stimON+1;
        else
            data.cntrTrials.(varToAlign{f})(:,stimOFF,:) = data.fullSession.(varToAlign{f})(:,j,:);
            stimOFF = stimOFF+1;
        end
        
        if (data.trial{j}.Performance.StimType < midPoint); %GoTrials Stim
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
        elseif (data.trial{j}.Performance.StimType > midPoint); %NoGoTrials Stim
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
        end
        
        for k = 1:length(data.angleOrientations)
            if data.trial{j}.Performance.StimType==data.angleOrientations(k); % distractor 1
                data.(orientationfieldname{k}).(varToAlign{f})(:,size((data.(orientationfieldname{k}).(varToAlign{f})),2)+1,:) = data.fullSession.(varToAlign{f})(:,j,:);
            end
        end
    end
    
    %     fieldLabels = fieldnames(data.AlignedDATA.TrialType);
    %     for jj = 1:length(fieldnames( data.AlignedDATA.TrialType))
    %
    %         if ~isempty(data.AlignedDATA.TrialType.(fieldLabels{jj}))
    %         A = data.AlignedDATA.TrialType.(fieldLabels{jj}).(varToAlign{f});
    %         A(:,find(sum(abs(A))==0),:)=[];
    %         data.AlignedDATA.TrialType.(fieldLabels{jj}).(varToAlign{f})=A;
    %
    %         end
    %     end
    %
    data = orderfields(data);
    
    %Calculate Block Performance
    blockBinSize = 12;
    blockPerformance = data.fullSession.(varToAlign{f})(1,:,7);
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
    end
    
end

%%

fontSize = 16;
smoothVal=201;
colorTally=1;
lineWidthPlots = 4;
f=1;
TimeScaled= TimePointsMin;
currFig = nan;

for vp = 1:length(varToPlot);
    currFig(vp) = figure(vp+1);
    set(currFig(vp),'name',varToPlot{vp},'numbertitle','off');
    set(currFig(vp),'Position',positionGraph1)
    clf;
    %raster of licking
    clear raster_matrix;
    spike_matrix = data.(varToPlot{vp}).(varToAlign{f})(:,:,2)';
    
    for i=1:size((spike_matrix),1)
        raster_matrix(i,:) = spike_matrix(i,:).*i;
    end
    
    subplot(5,3,[1 2 4 5 7 8])
    %Velocity
    if (includeVelocityMap==1)
        imagesc(TimeScaled,1,data.(varToPlot{vp}).(varToAlign{f})(:,:,1)')
        hold on
    end
    plot(TimeScaled,raster_matrix,'k.','LineWidth',4);
    hold on
    startLine=nan(size(data.(varToPlot{vp}).(varToAlign{f}),2),1);
    endLine=nan(size(data.(varToPlot{vp}).(varToAlign{f}),2),1);
    lickLine=nan(size(data.(varToPlot{vp}).(varToAlign{f}),2),1);
    
    for k = 1:size(data.(varToPlot{vp}).(varToAlign{f}),2);
        startLine(k,1) = min(find(data.(varToPlot{vp}).(varToAlign{f})(:,k,5)==1));
        endLine(k,1) = max(find(data.(varToPlot{vp}).(varToAlign{f})(:,k,5)==1));
        if ~isempty(find(data.(varToPlot{vp}).(varToAlign{f})(:,k,6)==1,1))
            lickLine(k,1) = min(find(data.(varToPlot{vp}).(varToAlign{f})(:,k,6)==1));
        else
            lickLine(k,1) = nan;
        end
    end
    
    for k = 1:size(data.(varToPlot{vp}).(varToAlign{f}),2);
        line([TimeScaled(startLine(k)) TimeScaled(startLine(k))], [k-1 k],'Color','g', 'LineWidth',2);
        line([TimeScaled(endLine(k)) TimeScaled(endLine(k))], [k-1 k],'Color','g', 'LineWidth',2);
        if (isnan(lickLine(k,1))==0) %if it doesn't equal nan
            line([TimeScaled(lickLine(k,1)) TimeScaled(lickLine(k,1))], [k-1 k],'Color','r', 'LineWidth',3);
        end
    end
    set(gca,'YDir','reverse');
    set(gca,'FontSize',16);
    xlim([min(TimeScaled), max(TimeScaled)]);
    ylim([1 size(data.(varToPlot{vp}).(varToAlign{f}),2)])
    ylabel('Trial number', 'FontSize', fontSize,'fontWeight','bold');
    xlabel('Time (ms)', 'FontSize', fontSize,'fontWeight','bold');
    
    subplot(5,3,[10 11]) %Plot Average Velocity
    tempVelDATA = data.(varToPlot{vp}).(varToAlign{f})(:,:,1);
     AvgVelocity = nanmean(tempVelDATA,2);
     semVelocity = nanstd(tempVelDATA,0,2)/sqrt(size(tempVelDATA,2));
%     AvgVelocity = nanmean(data.(varToPlot{vp}).(varToAlign{f})(:,:,1),2);
%     semVelocity = std(data.(varToPlot{vp}).(varToAlign{f})(:,:,1),2)/sqrt(size(data.(varToPlot{vp}).(varToAlign{f})(:,:,1),2));
    
    [running_speed,time_bins] = JB_computeSpeed(AvgVelocity,sampleRateSpacing,TimeScaled);
    [running_speedSEM] = JB_computeSpeed(semVelocity,sampleRateSpacing,TimeScaled);
    
    shadedErrorBar(time_bins,running_speed,running_speedSEM,{'-k'})
%     plot(time_bins,running_speed,'k','LineWidth',3);
%     hold on
%     errorbar(running_speed,running_speedSEM)
    %     smoothedVelocity = smooth(AvgVelocity,smoothVal);
    %     plot(TimeScaled,smoothedVelocity,'k','LineWidth',3);
    set(gca,'FontSize',16);
    %     xlim([min(TimeScaled), max(TimeScaled)])
    ylabel('Velocity cm/s', 'FontSize', fontSize,'fontWeight','bold');
    xlabel('Time (ms)', 'FontSize', fontSize,'fontWeight','bold');
    %ylim([0 20]);
    
    subplot(5,3,[13 14]) %Average licking
    binSize = 200;
    t = TimeScaled;
    tpoints =linspace(min(t),max(t),binSize);
    
        tempLickDATA = data.(varToPlot{vp}).(varToAlign{f})(:,:,2);
     AvgLicking = nanmean(tempLickDATA,2);
     semLicking = nanstd(tempLickDATA,0,2)/sqrt(size(tempLickDATA,2));
    [licking_rate,time_bins] = JB_computeSpeed(AvgLicking,sampleRateSpacing,TimeScaled);
    [licking_rateSEM] = JB_computeSpeed(semLicking,sampleRateSpacing,TimeScaled);
        shadedErrorBar(time_bins,licking_rate,licking_rateSEM,{'-k'})
        
    set(gca,'FontSize',16);
    %     xlim([min(TimeScaled), max(TimeScaled)])
    ylabel('Licks per s', 'FontSize', fontSize,'fontWeight','bold');
    xlabel('Time (ms)', 'FontSize', fontSize,'fontWeight','bold');
    
    subplot(5,3,[3 6 9])
    for i=1:size(data.(varToPlot{vp}).(varToAlign{f}),2)
        if data.(varToPlot{vp}).(varToAlign{f})(1,i,3)==1 && data.(varToPlot{vp}).(varToAlign{f})(1,i,4)<=midPoint; % Correct Lick on GoTrials stimulus = Hit
            line([0 1], [i i], 'Color', 'g','LineWidth',2);
            hold on;
        elseif data.(varToPlot{vp}).(varToAlign{f})(1,i,3)==0 && data.(varToPlot{vp}).(varToAlign{f})(1,i,4)<=midPoint; % No Lick on GoTrials stimulus = Miss
            line([0 1], [i i], 'Color', 'r','LineWidth',2);
            hold on;
        elseif data.(varToPlot{vp}).(varToAlign{f})(1,i,3)==0 && data.(varToPlot{vp}).(varToAlign{f})(1,i,4)>=midPoint; % No Lick on Distractor stimulus = Correct Rejection
            line([2 3], [i i], 'Color', 'g','LineWidth',2);
            hold on;
        elseif data.(varToPlot{vp}).(varToAlign{f})(1,i,3)==1 && data.(varToPlot{vp}).(varToAlign{f})(1,i,4)>=midPoint; % Lick on Distractor stimulus = False Alarm
            line([2 3], [i i], 'Color', 'r','LineWidth',2);
            hold on;
        end
    end
    
    xlim([-0.5 3.5])
    ylim([1 size(data.(varToPlot{vp}).(varToAlign{f}),2)])
    set(gca,'YDir','reverse');
    
    if length(varToPlot)>1
        figure(99)
        set(gcf,'Position',positionGraph3)
        subplot(2,1,1) %Plot Average Velocity
        if strcmp(varToPlot{vp},'fullSession')
%             plot(time_bins,running_speed,'k','LineWidth',lineWidthPlots,'Color', 'w')
            %             plot(TimeScaled,smoothedVelocity,'LineWidth',lineWidthPlots,'Color', 'w');
        elseif strcmp(varToPlot{vp},'stimTrials') || strcmp(varToPlot{vp},'NoGoTrials')
                     shadedErrorBar(time_bins,running_speed,running_speedSEM,{'-r'})
%             plot(time_bins,running_speed,'k','LineWidth',lineWidthPlots,'Color', 'r')
            %             plot(TimeScaled,smoothedVelocity,'LineWidth',lineWidthPlots,'Color', 'r');
        elseif strcmp(varToPlot{vp},'cntrTrials') || strcmp(varToPlot{vp},'GoTrials')
                     shadedErrorBar(time_bins,running_speed,running_speedSEM,{'-k'})
%             plot(time_bins,running_speed,'k','LineWidth',lineWidthPlots,'Color', 'k')
            %             plot(TimeScaled,smoothedVelocity,'LineWidth',lineWidthPlots,'Color','k');
        end
        hold on
        set(gca,'FontSize',16);
        %         xlim([min(TimeScaled), max(TimeScaled)])
        ylabel('Velocity cm/s', 'FontSize', fontSize,'fontWeight','bold');
        xlabel('Time (ms)', 'FontSize', fontSize,'fontWeight','bold');
        legend(varToPlot)
        %ylim([0 20]);
        
        subplot(2,1,2) %Average licking
        if strcmp(varToPlot{vp},'fullSession')
            
              shadedErrorBar(time_bins,licking_rate,licking_rateSEM,{'-w'})
%             plot(time_bins,licking_rate,'k','LineWidth',lineWidthPlots,'Color', 'w')
           
        elseif strcmp(varToPlot{vp},'stimTrials') || strcmp(varToPlot{vp},'NoGoTrials')
             shadedErrorBar(time_bins,licking_rate,licking_rateSEM,{'-r'})
%             plot(time_bins,licking_rate,'k','LineWidth',lineWidthPlots,'Color', 'r')
            %              plot(TimeScaled,smoothedLicking,'LineWidth',lineWidthPlots,'Color', 'r')
            %   plot(tpoints,full(mean(sparse(1:length(t),bin,AvgLicking))),'LineWidth',lineWidthPlots,'Color', 'r')
        elseif strcmp(varToPlot{vp},'cntrTrials')  || strcmp(varToPlot{vp},'GoTrials')
            shadedErrorBar(time_bins,licking_rate,licking_rateSEM,{'-k'})
%             plot(time_bins,licking_rate,'k','LineWidth',lineWidthPlots,'Color', 'k')
            %              plot(TimeScaled,smoothedLicking,'LineWidth',lineWidthPlots,'Color', 'k')
            
            %      plot(tpoints,full(mean(sparse(1:length(t),bin,AvgLicking))),'LineWidth',lineWidthPlots,'Color', 'k')
        end
        
        hold on
        set(gca,'FontSize',16);
        %         xlim([min(TimeScaled), max(TimeScaled)])
        ylabel('Licks per s', 'FontSize', fontSize,'fontWeight','bold');
        xlabel('Time (ms)', 'FontSize', fontSize,'fontWeight','bold');
        legend(varToPlot)
        
        if (plotAngles==1)
            
            figure(100);
            set(gcf,'Position',positionGraph3)
            subplot(2,1,1) %Plot Average Velocity
            if strncmp(varToPlot{vp},'stimAngle',9)
                plot(time_bins,running_speed,'k','LineWidth',lineWidthPlots,'Color',  rgb(plotColor{colorTally}))
                %                 plot(TimeScaled,smoothedVelocity,'LineWidth',lineWidthPlots,'Color', rgb(plotColor{colorTally}));
                hold on
                %                  colorTally = colorTally+1;
            end
            hold on
            set(gca,'FontSize',16);
            %             xlim([min(TimeScaled), max(TimeScaled)])
            ylabel('Velocity cm/s', 'FontSize', fontSize,'fontWeight','bold');
            xlabel('Time (ms)', 'FontSize', fontSize,'fontWeight','bold');
            legend(orientationfieldname)
            %ylim([0 20]);
            
            subplot(2,1,2) %Average licking
            if strncmp(varToPlot{vp},'stimAngle',9) 
                plot(time_bins,licking_rate,'k','LineWidth',lineWidthPlots,'Color', rgb(plotColor{colorTally}))
                %                 plot(TimeScaled,smoothedLicking,'LineWidth',lineWidthPlots,'Color', rgb(plotColor{colorTally}));
                %                   plot(tpoints,full(mean(sparse(1:length(t),bin,AvgLicking))),'LineWidth',2,'Color', rgb(plotColor{colorTally}));
                hold on
                colorTally = colorTally+1;
            end
            hold on
            set(gca,'FontSize',16);
            %             xlim([min(TimeScaled), max(TimeScaled)])
            ylabel('Licks per s', 'FontSize', fontSize,'fontWeight','bold');
            xlabel('Time (ms)', 'FontSize', fontSize,'fontWeight','bold');
            legend(orientationfieldname)
        end
    end
    
end

figure(21)
set(gcf,'Position',positionGraph2)
subplot(2,1,1);
plot(cumulativePerformance,'*b');
ylabel('cumulative performance', 'FontSize', fontSize,'fontWeight','bold');
xlabel('Trial number', 'FontSize', fontSize,'fontWeight','bold');
ylim([0,1])

subplot(2,1,2);
plot(blockPerformancePlot,'-b','LineWidth',2);
ylabel('block performance', 'FontSize', fontSize,'fontWeight','bold');
xlabel('block (12) number', 'FontSize', fontSize,'fontWeight','bold');
ylim([0,1]);

end




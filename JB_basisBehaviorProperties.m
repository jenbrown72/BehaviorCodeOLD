
load('DATA.mat')
basisPropertiesID{1,1} = 'sessionType';
sessionType = 1;
basisPropertiesID{2,1} = 'date';
date = 2;
basisPropertiesID{3,1} = 'datenum';
datenum = 3;
basisPropertiesID{4,1} = 'sessionDuration';
sessionDuration = 4;
basisPropertiesID{5,1} = 'totalStepsPerMin';
totalStepsPerMin = 5;

basisProperties = {};


encoder0Pos = 1;
rawSessionTime = 7;



for i=1:length(DATA.allFiles);
    
    basisProperties{date,i} = DATA.allFiles{i}.date;
    basisProperties{datenum,i} = DATA.allFiles{i}.datenum;

    tempDATA = DATA.allFiles{i}.rawData;
    basisProperties{sessionDuration,i} = ((tempDATA(end,rawSessionTime)-tempDATA(1,rawSessionTime))/1000)/60; % Calculate duration of session in minutes
    threshold = 50000;
    aboveThreshold = find(tempDATA(:,encoder0Pos)>threshold);
    for j = 1:length(aboveThreshold)
        tempDATA(aboveThreshold(j),:)=nan;
    end
    
    totalSteps = nansum(diff(tempDATA(:,1))>0);
    basisProperties{totalStepsPerMin,i} = totalSteps/basisProperties{sessionDuration,i};
    
    %import as metaDATA cell array
    n=1;
    
    metaData = DATA.allFiles{i}.metaData;
    for j=1:length(metaData)
        
        match(j,1) = strcmp('Orientation Selected = ', metaData{j,1});
        
        if (match(j,1)==1) && (metaData{j,2}>0)
            tempAngles(n,1) = metaData{j,2};
            n=n+1;
        end
        
        %AutoReward
        if(strcmp('Auto reward = ', metaData{j,1}));
            autoSession = metaData{j,2};
        end
        
    end
    
    %Define session type
    tempsessionType = strcat('S',num2str(length(tempAngles)));
    
    if autoSession==1
        tempsessionType = strcat(tempsessionType,'auto');
    end
    
    basisProperties{sessionType,i} = tempsessionType;
    
end



numPoints = 1:1:length(basisProperties);
for j = 1:length(basisProperties);
    plot(numPoints(j),basisProperties{4,j},'*r')
    hold on
end

DatesX = basisProperties(2,:)

for k = 1:length(basisProperties)
    
    days{k} = basisProperties{2,k}(1:11)
    
end


[a n] = unique(days)

for k =1:length(a)
    
    for j = 1:length(days)
      strcmp(a{k},days{j})
      match(1,j) = strcmp(a{k},days{j});
    end
    
    idx = find(match>0);
    
    if length(idx)>1
        
        for h = 1:length(idx)
            time(idx(h)) = basisProperties{4,idx(h)}
        end
        
        [~,p] = max(time)
        
        for kk = 1:length(idx)
            
            if idx(kk)==p
                
                continue
                
            else
                basisProperties{:,idx(1)} %delect this column
                
            end
            
.

            
    end
        
        
end


    
    

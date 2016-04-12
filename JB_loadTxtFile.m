function [newTxtFileCount] = JB_loadTxtFile(MouseID)

% Initialize variables.
delimiter = ',';
endRow = 30;
startRow = 35;

% Read columns of data as strings:
% For more information, see the TEXTSCAN documentation.
formatSpecMetaDATA = '%s%s%[^\n\r]';
formatSpecTempDATA = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

startDir = 'C:\Users\adesniklab\Documents\BehaviorRawData\MouseID\';

if nargin==1; % if a filename was inputted into function
    
    currfile = MouseID;
    tempdataDir = char(strcat(startDir, currfile, '\txtFiles'));
    cd(tempdataDir);
    
    display(' ')
    display(' ')
    disp(['________________Analysing file: ' currfile ' __________________________________'])
    display(' ')
    display(' ')
end

%Check to see if files have already been loaded
DATA.loadedFiles = [];
DATA.allFiles = [];

%Initialise variables
% Angles used in behavior
stimAngles = [225, 241, 254, 263, 266, 268, 270, 272, 274, 277, 284, 299, 315];

currDir = cd;
dataDir = strcat(currDir,'\txtFiles\data');
cd(dataDir);
olderDir = dir('DATA.mat');

%if there is already a DATA.mat file - reload this
if ~isempty(olderDir)
    load('DATA.mat')
    loaded = 1;
else
    loaded = 0;
end

cd(currDir);
D = dir('*.txt'); %create a list of the txt files in the current directory

currDir = pwd; %identify the current folder
[~, deepestFolder,~] = fileparts(currDir); %returns the pathstr, name and extension. we care about the name
DATA.mouseID = deepestFolder;
newTxtFileCount = 0;

for i=1:size((D),1);
    filed = 0; %keep track of the txtfiles loaded
    if (loaded==1)
        if(any(strcmp(D(i).name, DATA.loadedFiles))) % if already loaded this data, move on
            disp(['DATA file: ' D(i).name ' Already Loaded'])
            filed = 1;
        end
    end
    
    if(filed==0)
        % Open the text file.
        fileID = fopen(D(i).name,'r');
        % Read columns of data according to format string.
        dataArray = textscan(fileID, formatSpecMetaDATA, endRow, 'Delimiter', delimiter, 'ReturnOnError', false);
        dataArrayTempDATA = textscan(fileID, formatSpecTempDATA, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
        % Close the text file.
        fclose(fileID);
        
        %For MetaDATA
        metaData = [dataArray{1:end-1}];
       % metaData = raw;
        
        %For tempDATA
        tempDATA = csvread(D(i).name,30,0); % Dont read first line - sometimes contains a string 'A', start from row 26 as the top is metaData
        tempDATA(length(tempDATA),:) = []; %Delete the last line incase not fully read
        
        disp(' ')
        disp(['Loading file: ' D(i).name])
        
        C = unique(tempDATA(:,11)); %Identify how many angles were presented
        stimNumber = size(find(ismember(stimAngles, C)),2);
        currFileNo = length(DATA.loadedFiles)+1;
        
        %Save into data format
        DATA.allFiles{currFileNo} = D(i);
        DATA.allFiles{currFileNo}.rawData = tempDATA;
        DATA.allFiles{currFileNo}.metaData = metaData;
        DATA.loadedFiles{1,currFileNo} = D(i).name;
        disp(['loaded: ' D(i).name ])
        disp(' ')
        filed = 1;
        newTxtFileCount = newTxtFileCount+1;
        
        %Extract date from filename
        tempName = DATA.allFiles{currFileNo}.name;
        startIdx = findstr('_201',tempName);
        endIdx = findstr('_Box',tempName);
        tempName = tempName(startIdx:endIdx);
        tempLocation = findstr(tempName,'_');
        
        %add a 0to the start of file date names so that sort with work later, e.g.
        %7 will be 07.
        for ii = 1:length(tempLocation)-1;
            if(tempLocation(ii+1)-tempLocation(ii))==2;
                startName = tempName(1:tempLocation(ii));
                endName = tempName(tempLocation(ii)+1:end);
                tempName = strcat(startName,'0',endName);
                tempLocation = tempLocation+1;
            end
        end
        toDelete = strfind(tempName,'_');
        tempName(toDelete)=[];
        DATA.allFiles{currFileNo}.dateFromFile = tempName;
    end
    
    if (filed==0)
        countNumber = max(tempDATA(:,6));
        if (countNumber<12)
            msgbox(strcat(D(i).name,' was not uploaded - count < 12'));
            disp(' ')
            disp(' ')
            disp(['WARNING..... NO FOLDER FOR: ' D(i).name ])
            disp(' ')
            disp(' ')
        else
            msgbox(strcat(D(i).name,' was not uploaded - stimNumber = ' , num2str(stimNumber)));
            disp(' ')
            disp(' ')
            disp(['WARNING..... NO FOLDER FOR: ' D(i).name ])
            disp(' ')
            disp(' ')
        end
    end
end

disp(['Upload complete'])
curr = cd;
cd(strcat(curr, '\txtFiles\data'))

if newTxtFileCount>0
    disp(['Ordering DATA files'])
    for i=1:length(DATA.allFiles);
        dataDate(i,1) = str2num(DATA.allFiles{i}.dateFromFile);
    end
    
    [~,S] = sort(dataDate);
    
    for j = 1:length(S)
        sortedFiles{j} = DATA.allFiles{S(j)};
        sortedLoadedFiles{j} = DATA.loadedFiles{S(j)};
    end
    
    %save files
    DATA.allFiles = sortedFiles;
    DATA.loadedFiles = sortedLoadedFiles;
    save('DATA.mat', 'DATA', '-v7.3')
    
else
    disp(['No new files to analyse'])
end

end
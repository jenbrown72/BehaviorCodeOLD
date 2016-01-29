function [DATA, loadedFiles] = JB_loadTxtFileTEST(MouseID)

% Initialize variables.
delimiter = ',';
endRow = 24;
startRow = 26;

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
% number of angles used in behavior
stimAngles = [225, 241, 254, 263, 266, 268, 270, 272, 274, 277, 284, 299, 315];

currDir = cd;
dataDir = strcat(currDir,'\txtFiles\data');
cd(dataDir);
olderDir = dir('DATA.mat');

if ~isempty(olderDir)
    load('DATA.mat')
    loaded = 1;
    
else
    loaded = 0;
end

cd(currDir);
D = dir('*.txt');

currDir = pwd;
[~, deepestFolder,~] = fileparts(currDir);
DATA.mouseID = deepestFolder;

for i=1:size((D),1);
    
    filed = 0;
    
    if (loaded==1)
        
        if(any(strcmp(D(i).name, DATA.loadedFiles))) % if already loaded this data, move on
            disp(['DATA file: ' D(i).name ' Already Loaded'])
            filed = 1;
            %   continue
            
        end
        
    end
    
    if(filed==0)
        
        %Load meta data
        % Open the text file.
        fileID = fopen(D(i).name,'r');
        % Read columns of data according to format string.
        dataArray = textscan(fileID, formatSpecMetaDATA, endRow, 'Delimiter', delimiter, 'ReturnOnError', false);
        dataArrayTempDATA = textscan(fileID, formatSpecTempDATA, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
        % Close the text file.
        fclose(fileID);
        
        %For MetaDATA
        % Convert the contents of columns containing numeric strings to numbers.
        % Replace non-numeric strings with NaN.
        raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
        for col=1:length(dataArray)-1
            raw(1:length(dataArray{col}),col) = dataArray{col};
        end
        numericData = NaN(size(dataArray{1},1),size(dataArray,2));
        
        % Converts strings in the input cell array to numbers. Replaced non-numeric
        % strings with NaN.
        rawData = dataArray{2};
        for row=1:size(rawData, 1);
            % Create a regular expression to detect and remove non-numeric prefixes and
            % suffixes.
            regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
            try
                result = regexp(rawData{row}, regexstr, 'names');
                numbers = result.numbers;
                
                % Detected commas in non-thousand locations.
                invalidThousandsSeparator = false;
                if any(numbers==',');
                    thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                    if isempty(regexp(thousandsRegExp, ',', 'once'));
                        numbers = NaN;
                        invalidThousandsSeparator = true;
                    end
                end
                % Convert numeric strings to numbers.
                if ~invalidThousandsSeparator;
                    numbers = textscan(strrep(numbers, ',', ''), '%f');
                    numericData(row, 2) = numbers{1};
                    raw{row, 2} = numbers{1};
                end
            catch me
            end
        end
        
        % Create output variable
        metaData = raw;
        
        %For tempDATA
        tempDATA = csvread(D(i).name,26,0); % Dont read first line - sometimes contains a string 'A'
        tempDATA(length(tempDATA),:) = []; %Delete the last line incase not fully read
        
        %  end
        
        disp(' ')
        disp(['Loading file: ' D(i).name])
        
        C = unique(tempDATA(:,11));
        stimNumber = size(find(ismember(stimAngles, C)),2);
        
        currFileNo = length(DATA.loadedFiles)+1;
        
        DATA.allFiles{currFileNo} = D(i);
        DATA.allFiles{currFileNo}.rawData = tempDATA;
        DATA.allFiles{currFileNo}.metaData = metaData;
        
        DATA.loadedFiles{1,currFileNo} = D(i).name;
        disp(['loaded: ' D(i).name ])
        disp(' ')
        filed = 1;
        
             
                
                
        
        %Extract date from filename
        
        tempName = DATA.allFiles{currFileNo}.name;
        startIdx = findstr('_201',tempName);  
endIdx = findstr('_Box',tempName);
tempName = tempName(startIdx:endIdx);
tempLocation = findstr(tempName,'_');


for i = 1:length(tempLocation)-1;
    
    if(tempLocation(i+1)-tempLocation(i))==2;
        
        startName = tempName(1:tempLocation(i));
        endName = tempName(tempLocation(i)+1:end);
        
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
%save('DATA.mat', 'DATA', '-v7.3')

disp(['Ordering DATA files'])

for i=1:length(DATA.allFiles);
    
    dataDate(i,1) = str2num(DATA.allFiles{i}.dateFromFile);
    
end

[~,S] = sort(dataDate);

for j = 1:length(S)
    sortedFiles{j} = DATA.allFiles{S(j)};
    sortedLoadedFiles{j} = DATA.loadedFiles{S(j)};
end

DATA.allFiles = sortedFiles;
DATA.loadedFiles = sortedLoadedFiles;
save('DATA.mat', 'DATA', '-v7.3')

end
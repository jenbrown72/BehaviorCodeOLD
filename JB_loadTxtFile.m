function [DATA, loadedFiles] = JB_loadTxtFile(MouseID)

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

%Initialise variables
% number of angles used in behavior
angleNumber = {'Auto','S1','S2','S6','S12'};
stimAngles = [225, 241, 254, 263, 266, 268, 270, 272, 274, 277, 284, 299, 315];

for i = 1:length(angleNumber)
    DATA.(char((angleNumber{i}))) = [];
end


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

for i=1:size((D),1);
    
    filed = 0;
    %    for i=1;
    
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
        % This call is based on the structure of the file used to generate this
        % code. If an error occurs for a different file, try regenerating the code
        % from the Import Tool.
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
        
%         %For TempDATA
%         
%         % Convert the contents of columns containing numeric strings to numbers.
% % Replace non-numeric strings with NaN.
% raw = repmat({''},length(dataArrayTempDATA{1}),length(dataArrayTempDATA)-1);
% for col=1:length(dataArrayTempDATA)-1
%     raw(1:length(dataArrayTempDATA{col}),col) = dataArrayTempDATA{col};
% end
% numericData = NaN(size(dataArrayTempDATA{1},1),size(dataArrayTempDATA,2));
% 
% for col=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
%     % Converts strings in the input cell array to numbers. Replaced non-numeric
%     % strings with NaN.
%     rawData = dataArrayTempDATA{col};
%     for row=1:size(rawData, 1);
%         % Create a regular expression to detect and remove non-numeric prefixes and
%         % suffixes.
%         regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
%         try
%             result = regexp(rawData{row}, regexstr, 'names');
%             numbers = result.numbers;
%             
%             % Detected commas in non-thousand locations.
%             invalidThousandsSeparator = false;
%             if any(numbers==',');
%                 thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
%                 if isempty(regexp(thousandsRegExp, ',', 'once'));
%                     numbers = NaN;
%                     invalidThousandsSeparator = true;
%                 end
%             end
%             % Convert numeric strings to numbers.
%             if ~invalidThousandsSeparator;
%                 numbers = textscan(strrep(numbers, ',', ''), '%f');
%                 numericData(row, col) = numbers{1};
%                 raw{row, col} = numbers{1};
%             end
%         catch me
%         end
%     end
% end
% 
% % Replace non-numeric cells with NaN
% R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
% raw(R) = {NaN}; % Replace non-numeric cells
% 
% % Create output variable
% tempDATA = cell2mat(raw);
% 
%     end
%     
% end


        
        tempDATA = csvread(D(i).name,26,0); % Dont read first line - sometimes contains a string 'A'
        tempDATA(length(tempDATA),:) = []; %Delete the last line incase not fully read
        
    end
        
        disp(' ')
        disp(['Loading file: ' D(i).name])
        
        C = unique(tempDATA(:,11));
        stimNumber = size(find(ismember(stimAngles, C)),2);
        
        if tempDATA(1,8) ==1; %% auto trial
            
            Key = 'Auto';
            Index = strfind(angleNumber,Key);
            
            for f=1:length(Index);
                
                if Index{f}==1
                    currFileNo = length(DATA.Auto)+1;
                    
                    DATA.(angleNumber{f}){currFileNo} = D(i);
                    DATA.(angleNumber{f}){currFileNo}.rawData = tempDATA;
                    DATA.(angleNumber{f}){currFileNo}.metaData = metaData;
                    
                    DATA.loadedFiles{i,1} = D(i).name;
                    DATA.loadedFiles{i,2} = D(i).datenum;
                    DATA.loadedFiles{i,3} = D(i).date;
                    
                    disp([(angleNumber{f}) ': ' D(i).name ])
                    disp(' ')
                    filed = 1;
                    
                else
                    continue
                end
                
            end
            
        else
            
            for k=1:length(angleNumber)
                
                Key = 'S';
                Index = strfind(angleNumber{k},Key);
                
                if ~isempty(Index)
                    
                    Value = sscanf(angleNumber{k}(Index(1) + length(Key):end), '%g', 1);
                    
                    if (Value==stimNumber)
                        
                        disp([angleNumber{k} ': ' D(i).name ])
                        disp(' ')
                        
                        currFileNo = length(DATA.(angleNumber{k}))+1;
                        
                        DATA.(angleNumber{k}){currFileNo} = D(i);
                        DATA.(angleNumber{k}){currFileNo}.rawData = tempDATA;
                        DATA.(angleNumber{k}){currFileNo}.metaData = metaData;
                        
                        DATA.loadedFiles{i,1} = D(i).name;
                        DATA.loadedFiles{i,2} = D(i).datenum;
                        DATA.loadedFiles{i,3} = D(i).date;
                        filed = 1;
                        
                    end
                    
                else
                    continue
                end
                
            end
            
        end
        
   % end
    
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
save('DATA.mat', 'DATA', '-v7.3')

end
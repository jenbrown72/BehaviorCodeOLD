function [data] = JB_LoadandAnalyseBehaviorfiles(currfile, fileSession)

% To analyse a particular filename: data = JB_Loadbehaviorfiles('JB00132_022415.mat')
% To analyse all files in current directory data = JB_Loadbehaviorfiles


load('DATA.mat');

if nargin==1; % if a filename was inputted into function
    
   % currfile = filename;
    display(' ')
    display(' ')
    disp(['________________Analysing file: ' currfile ' __________________________________'])
    display(' ')
    display(' ')
    
    data = JB_sessionAnalysisTEST(currfile, fileSession);
    
elseif nargin==2; % if a file session was inputted into function
        
    data = JB_sessionAnalysisTEST(currfile, fileSession);
    
elseif nargin==0; %if no filenames were inputted find out with a file called data.mat already exsists in the current directory
    
    dataFile = dir('data.mat');
    
    if ~isempty(dataFile); % if there is a data.mat files in the current directory - load
        
        load('data.mat')
        Directory = dir('*.mat');
        fileNames = struct2cell(Directory);
        
        for i=1:size((fileNames),2);
            match=0;
            currfile = fileNames{1,i};
            
            display(' ')
            disp(['________________Checking file: ' currfile ' __________________________________'])
            display(' ')
            
            if strcmp(currfile, 'data.mat'); % ignore this file
                continue;
                
            end
            
            for j = 1:size((data),2)
                
                if (strcmp(currfile(1,1:end-4),fieldnames(data(j).raw))==1)
                    match=1;
                    
                else
                    continue;
                    
                end
            end
            
            if match==1
                
                continue
                
            else
                
                display(' ')
                disp(['________________Analysing file: ' currfile ' __________________________________'])
                display(' ')
                
                data(size((data),2)+1) = JB_sessionAnalysis(currfile);  % if this is a new file = annalyse and store results in data.mat
                flag=1;
                continue
                
            end
            
        end
        
    else
        
        %find all csv files in this folder
        Directory = dir('*.mat');
        
        
        %load csv files
        for i=1:length(Directory)
            % prompt = 'Load Next Session? ';
            fileNames = struct2cell(Directory);
            
            currfile = fileNames{1,i};
            
            display(' ')
            display('Load  Next Session?')
            pause
            
            display(' ')
            disp(['________________Analysing file: ' currfile ' __________________________________'])
            display(' ')
            data(i) = JB_sessionAnalysis(currfile);
            
            
        end
        
    end
    
end

clear dataFile currfile i fileNames Directory;

save('data.mat', '-v7.3') ;
%save('data.mat')

end

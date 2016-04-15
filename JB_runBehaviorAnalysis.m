function [] = JB_runBehaviorAnalysis(plotON)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
if  nargin==1 
    plotON = plotON;
else
    plotON=0;
end


currDirectory = cd;
dirFolders = dir;
folders = dirFolders(arrayfun(@(x)x.name(1),dirFolders)~='.');

for h = 1:length(folders) %for each mouse ID folder
    
    tempMousedataDir = char(strcat(currDirectory,'\',folders(h).name)); 
    cd(tempMousedataDir); %go to the specific mouse folder
    
    [newTxtFileCount] = JB_loadTxtFile; %run load all txt files, will only upload new files
    
  if newTxtFileCount>0 %if new files are found, run basic BehaviorProperties
      
   JB_basicBehaviorProperties(0)
         
  end
  
cd(currDirectory);
end


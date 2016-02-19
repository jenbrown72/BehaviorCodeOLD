function [] = JB_runBehaviorAnalysis
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

currDirectory = cd;
dirFolders = dir;
folders = dirFolders(arrayfun(@(x)x.name(1),dirFolders)~='.');

for h = 1:length(folders) %for each mouse ID folder
    
    tempMousedataDir = char(strcat(currDirectory,'\',folders(h).name)); 
    cd(tempMousedataDir); %go to the specific mouse folder
    
    [newTxtFileCount] = JB_loadTxtFile; %run load all txt files, will only upload new files
    
  if newTxtFileCount>0 %if new files are found, run basic BehaviorProperties
    
    JB_basicBehaviorProperties(0); %0 dont plot figures, 1 plot figures
     
  end
    
currDirectory;
end


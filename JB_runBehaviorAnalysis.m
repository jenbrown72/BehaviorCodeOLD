function [] = JB_runBehaviorAnalysis
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

currDirectory = cd;
dirFolders = dir;
folders = dirFolders(arrayfun(@(x)x.name(1),dirFolders)~='.');

for h = 1:length(folders)
    
    tempMousedataDir = char(strcat(currDirectory,'\',folders(h).name));
    cd(tempMousedataDir);
    
    JB_loadTxtFileTEST;
    
    JB_basicBehaviorProperties;
     
currDirectory;
end


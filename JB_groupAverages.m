function [AllDATA] = JB_groupAverages(AllDATA,listToAnalyse,type)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%listToAnalyse = {'JB022_RL_';'SG020_RC_';'JB016_R_';'JB017_RL_';'JB017_RC_';'JB017_NT_'}
%listToAnalyse = {'JB022_R_';'JB022_NT_';'SG020_RL_';'SG020_NT_';'JB016_CE_';'JB017_R_';'JB017_L_'}
currDirectory = cd;

for h = 1:length(listToAnalyse) %for each mouse ID folder
    
    cd(char(strcat(currDirectory,'\',listToAnalyse{h},'\txtFiles\data')));  %go to the specific mouse folder
    
    [basicPropertiesToPlot] = JB_basicBehaviorPropertiesNeat(1); %run basic properties
    [DATAtrim] = JB_plotTrimming(basicPropertiesToPlot,1,1);
    
    AllDATA{type}.data{h} = DATAtrim;
    AllDATA{type}.data{h}.name = listToAnalyse{h};
    cd(currDirectory);
end

end


function [AllDATA] = JB_groupAverages(AllDATA,listToAnalyse,condition)
%JB_groupAverages Selects already analysed data from independent mouse IDs
%for group analysis

%funtion inputs

% AllDATAL data structure to collect, if concaternation to an already
% generated AllDATA at structure name to function input, or leave empty []:
% listToAnalyse: Mouse IDs to select from
% eg:
% listToAnalyse = {'JB022_RL_';'SG020_RC_';'JB016_R_';'JB017_RL_';'JB017_RC_';'JB017_NT_'}
% listToAnalyse = {'JB022_R_';'JB022_NT_';'SG020_RL_';'SG020_NT_';'JB016_CE_';'JB017_R_';'JB017_L_'}
% listToAnalyse = {'SG022_RC_';'SG022_RL_'};
% listToAnalyse = {'JB022_RL_';'JB022_R_';'JB022_NT_';'SG020_RL_';'SG020_RC_';'SG020_NT_'};
% listToAnalyse = {'SG022_RC_';'SG022_RL_';'SG023_L_'};
%concat if from different directories:  C = [AllDATA{1,1}.data AllDATA2{1,1}.data]

% listToAnalyse = {'JB022_RL_';'JB022_R_';'JB022_NT_';'SG020_RL_';'SG020_RC_';'SG020_NT_';'JB016_R_';'JB016_CE_';'JB016_NT_';'JB017_RL_';'JB017_RC_';'JB017_R_';'JB017_L_';'JB017_NT_'};

%Stim
%listToAnalyse = {'JB016_CE_';'JB016_L_';'JB016_NT_';'JB016_R_';'JB016_RL_'};

% condition: puts different conditions (e.g. stim on/ stim off) in different structures

% e.g. [AllDATA] = JB_groupAverages([],listToAnalyse,1)
% e.g. [AllDATA] = JB_groupAverages([AllDATA],listToAnalyse,2)
% e.g. [AllDATA] = JB_groupAverages([AllDATA],{'JB023_R_';'JB024_NT_',2)


currDirectory = cd;

for h = 1:length(listToAnalyse) %for each mouse ID folder
    
    cd(char(strcat(currDirectory,'\',listToAnalyse{h},'\txtFiles\data')));  %go to the specific mouse folder
    
    %Add here whish functions you want to run
    [basicPropertiesToPlot,possibleAngles] = JB_basicBehaviorProperties(0,0); %run basic properties
   % [DATAtrim] = JB_plotTrimmingTest(basicPropertiesToPlot,1,0);
%     [DATA] = JB_plotOptogenetics(basicPropertiesToPlot,possibleAngles,1);

    
    %save data into structure
%      AllDATA{condition}.data{h} = DATA;
    AllDATA{condition}.data{h} = basicPropertiesToPlot;
%     AllDATA{condition}.data{h}.name = listToAnalyse{h};
    cd(currDirectory);
end

end

